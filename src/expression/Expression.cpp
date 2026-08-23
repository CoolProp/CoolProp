#include "CoolProp/expression/Expression.h"
#include "CoolProp/expression/detail/Lexer.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/strings.h"
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace CoolProp {
namespace expression {
namespace detail {

// Function intrinsics supported by the grammar.
// NOLINTNEXTLINE(performance-enum-size) -- already uses std::uint8_t; false positive on some clang-tidy versions
enum class Func : std::uint8_t
{
    exp,
    ln,
    log10,
    sqrt,
    abs,
    pow,
    sinh,
    cosh,
    tanh,
    sin,
    cos,
    atan
};

// Evaluation context: scalars indexed by slot, arrays indexed by slot, and the
// single active summation index `i` (v1 forbids nested sums, so one is enough).
struct EvalState
{
    // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
    // non-owning views; valid for the bounded eval lifetime
    std::vector<double>& scalars;
    const std::vector<std::vector<double>>& arrays;
    // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
    long i = 0;
};

struct Node
{
    Node() = default;
    virtual ~Node() = default;
    Node(const Node&) = delete;
    Node& operator=(const Node&) = delete;
    Node(Node&&) = delete;
    Node& operator=(Node&&) = delete;
    virtual double eval(EvalState& st) const = 0;
};
using NodePtr = std::unique_ptr<Node>;

struct NumNode : Node
{
    double v = 0.0;
    explicit NumNode(double value) : v(value) {}
    double eval(EvalState&) const override {
        return v;
    }
};

struct ScalarNode : Node  // intrinsic, constant, derived, or let — all read a slot
{
    int slot = 0;
    explicit ScalarNode(int s) : slot(s) {}
    double eval(EvalState& st) const override {
        return st.scalars[slot];
    }
};

struct IndexNode : Node  // bare summation index used as a value
{
    double eval(EvalState& st) const override {
        return static_cast<double>(st.i);
    }
};

struct ArrayNode : Node  // arr[i]
{
    int slot = 0;
    explicit ArrayNode(int s) : slot(s) {}
    double eval(EvalState& st) const override {
        return st.arrays[slot][static_cast<std::size_t>(st.i)];
    }
};

struct NegNode : Node
{
    NodePtr child{};
    explicit NegNode(NodePtr c) : child(std::move(c)) {}
    double eval(EvalState& st) const override {
        return -child->eval(st);
    }
};

struct BinNode : Node
{
    char op = '+';  // '+','-','*','/','^'
    NodePtr l{}, r{};
    BinNode(char o, NodePtr a, NodePtr b) : op(o), l(std::move(a)), r(std::move(b)) {}
    double eval(EvalState& st) const override {
        double a = l->eval(st), b = r->eval(st);
        // NOLINTBEGIN(bugprone-branch-clone) -- distinct arithmetic ops; clang-tidy false-positive
        switch (op) {
            case '+':
                return a + b;
            case '-':
                return a - b;
            case '*':
                return a * b;
            case '/':
                return a / b;
            case '^':
                return std::pow(a, b);
            default:
                throw ValueError("internal: invalid binary operator");
        }
        // NOLINTEND(bugprone-branch-clone)
    }
};

struct CallNode : Node
{
    Func fn = Func::exp;
    std::vector<NodePtr> args{};
    CallNode(Func f, std::vector<NodePtr> a) : fn(f), args(std::move(a)) {}
    double eval(EvalState& st) const override {
        double x = args[0]->eval(st);
        // NOLINTBEGIN(bugprone-branch-clone) -- distinct std:: calls; clang-tidy false-positives them as identical
        switch (fn) {
            case Func::exp:
                return std::exp(x);
            case Func::ln:
                return std::log(x);
            case Func::log10:
                return std::log10(x);
            case Func::sqrt:
                return std::sqrt(x);
            case Func::abs:
                return std::fabs(x);
            case Func::pow:
                return std::pow(x, args[1]->eval(st));
            case Func::sinh:
                return std::sinh(x);
            case Func::cosh:
                return std::cosh(x);
            case Func::tanh:
                return std::tanh(x);
            case Func::sin:
                return std::sin(x);
            case Func::cos:
                return std::cos(x);
            case Func::atan:
                return std::atan(x);
            default:
                throw ValueError("internal: invalid function");
        }
        // NOLINTEND(bugprone-branch-clone)
    }
};

struct SumNode : Node
{
    long n = 0;  // fixed array length, resolved at bind time
    NodePtr body{};
    SumNode(long length, NodePtr b) : n(length), body(std::move(b)) {}
    double eval(EvalState& st) const override {
        double s = 0.0;
        long saved = st.i;
        for (long k = 0; k < n; ++k) {
            st.i = k;
            s += body->eval(st);
        }
        st.i = saved;
        return s;
    }
};

// ---------------------------------------------------------------------------
// Parse-time-only nodes: identifiers are unresolved (by name). The binder
// (Task 5) rewrites each into a runtime node once slots are known.
// Their eval() must never run.
// ---------------------------------------------------------------------------

struct NameNode : Node  // a bare identifier: intrinsic/constant/derived/let, or sum index
{
    std::string name{};
    explicit NameNode(std::string nm) : name(std::move(nm)) {}
    double eval(EvalState&) const override {
        throw ValueError("internal: unbound NameNode evaluated");
    }
};

struct ArrayRefNode : Node  // name[index]; index identifier captured for validation
{
    std::string arrayName{}, indexName{};
    ArrayRefNode(std::string a, std::string idx) : arrayName(std::move(a)), indexName(std::move(idx)) {}
    double eval(EvalState&) const override {
        throw ValueError("internal: unbound ArrayRefNode evaluated");
    }
};

struct ParseSumNode : Node  // sum(indexName: body); length resolved at bind
{
    std::string indexName{};
    NodePtr body{};
    double eval(EvalState&) const override {
        throw ValueError("internal: unbound ParseSumNode evaluated");
    }
};

struct NamedCallNode : Node  // funcName(args...); resolved to Func at bind
{
    std::string name{};
    std::vector<NodePtr> args{};
    NamedCallNode(std::string nm, std::vector<NodePtr> a) : name(std::move(nm)), args(std::move(a)) {}
    double eval(EvalState&) const override {
        throw ValueError("internal: unbound NamedCallNode evaluated");
    }
};

struct LetStmt
{
    std::string name{};
    NodePtr expr{};
    std::size_t col = 0;
};

struct ParseResult
{
    std::vector<LetStmt> lets;
    NodePtr result;
};

// ---------------------------------------------------------------------------
// Pratt (recursive-descent) parser
// ---------------------------------------------------------------------------

class Parser
{
   public:
    explicit Parser(std::vector<Token> toks) : t(std::move(toks)) {}

    ParseResult parse() {
        ParseResult pr;
        skipSeps();
        while (peek().type == TokenType::Keyword_let) {
            pr.lets.push_back(parseLet());
            skipSeps();
        }
        pr.result = parseExpr(0);
        skipSeps();
        expect(TokenType::End, "trailing tokens after final expression");
        return pr;
    }

   private:
    std::vector<Token> t{};
    std::size_t pos = 0;

    [[nodiscard]] const Token& peek() const {
        return t[pos];
    }
    const Token& advance() {
        return t[pos++];
    }
    bool match(TokenType tt) {
        if (t[pos].type == tt) {
            ++pos;
            return true;
        }
        return false;
    }
    void skipSeps() {
        while (t[pos].type == TokenType::StmtSep)
            ++pos;
    }
    const Token& expect(TokenType tt, const char* what) {
        if (t[pos].type != tt) throw ValueError(format("expected %s at col %d", what, (int)t[pos].col));
        return t[pos++];
    }

    LetStmt parseLet() {
        std::size_t col = advance().col;  // consume 'let'
        const Token& id = expect(TokenType::Ident, "identifier after 'let'");
        std::string name = id.text;
        expect(TokenType::Equals, "'=' in let binding");
        NodePtr e = parseExpr(0);
        if (peek().type != TokenType::StmtSep && peek().type != TokenType::End)
            throw ValueError(format("expected end of let statement at col %d", (int)peek().col));
        return LetStmt{name, std::move(e), col};
    }

    [[nodiscard]] int lbp(TokenType tt) const {
        switch (tt) {
            case TokenType::Plus:
            case TokenType::Minus:
                return 10;
            case TokenType::Star:
            case TokenType::Slash:
                return 20;
            case TokenType::Caret:
                return 30;
            default:
                return 0;
        }
    }

    [[nodiscard]] char opChar(TokenType tt) const {
        switch (tt) {
            case TokenType::Plus:
                return '+';
            case TokenType::Minus:
                return '-';
            case TokenType::Star:
                return '*';
            case TokenType::Slash:
                return '/';
            case TokenType::Caret:
                return '^';
            default:
                return '?';
        }
    }

    NodePtr parseExpr(int minbp) {
        NodePtr left = parsePrefix();
        for (;;) {
            TokenType tt = peek().type;
            int bp = lbp(tt);
            if (bp == 0 || bp < minbp) break;
            char op = opChar(tt);
            advance();
            int nextMin = (tt == TokenType::Caret) ? bp : bp + 1;  // ^ right-assoc
            NodePtr right = parseExpr(nextMin);
            left = std::make_unique<BinNode>(op, std::move(left), std::move(right));
        }
        return left;
    }

    NodePtr parsePrefix() {
        if (match(TokenType::Minus)) return std::make_unique<NegNode>(parseExpr(25));  // unary minus: looser than ^ (30), tighter than +/- (10)
        return parsePrimary();
    }

    NodePtr parsePrimary() {
        const Token& tk = peek();
        if (tk.type == TokenType::Number) {
            double v = advance().number;
            return std::make_unique<NumNode>(v);
        }
        if (tk.type == TokenType::LParen) {
            advance();
            NodePtr e = parseExpr(0);
            expect(TokenType::RParen, "')'");
            return e;
        }
        if (tk.type == TokenType::Keyword_sum) return parseSum();
        if (tk.type == TokenType::Ident) {
            std::string name = advance().text;
            if (peek().type == TokenType::LParen) return parseCall(name);
            if (match(TokenType::LBracket)) {
                const Token& idx = expect(TokenType::Ident, "index identifier inside '[]'");
                std::string idxname = idx.text;
                expect(TokenType::RBracket, "']'");
                return std::make_unique<ArrayRefNode>(name, idxname);
            }
            return std::make_unique<NameNode>(name);
        }
        throw ValueError(format("unexpected token at col %d", (int)tk.col));
    }

    NodePtr parseSum() {
        advance();  // 'sum'
        expect(TokenType::LParen, "'(' after sum");
        const Token& idx = expect(TokenType::Ident, "index identifier in sum");
        std::string idxname = idx.text;
        expect(TokenType::Colon, "':' in sum");
        NodePtr body = parseExpr(0);
        expect(TokenType::RParen, "')' to close sum");
        auto s = std::make_unique<ParseSumNode>();
        s->indexName = idxname;
        s->body = std::move(body);
        return s;
    }

    NodePtr parseCall(const std::string& name) {
        advance();  // '('
        std::vector<NodePtr> args;
        if (peek().type != TokenType::RParen) {
            args.push_back(parseExpr(0));
            while (match(TokenType::Comma))
                args.push_back(parseExpr(0));
        }
        expect(TokenType::RParen, "')' to close call");
        return std::make_unique<NamedCallNode>(name, std::move(args));
    }
};

// ---------------------------------------------------------------------------

std::vector<Token> lex(const std::string& s) {
    std::vector<Token> out;
    std::size_t i = 0, n = s.size();
    auto col = [&](std::size_t p) { return p + 1; };
    while (i < n) {
        char c = s[i];
        if (c == ' ' || c == '\t' || c == '\r') {
            ++i;
            continue;
        }
        if (c == '\n' || c == ';') {
            out.push_back({TokenType::StmtSep, 0.0, "", col(i)});
            ++i;
            continue;
        }
        if (std::isdigit(static_cast<unsigned char>(c)) != 0 || c == '.') {
            const char* start = s.c_str() + i;
            char* end = nullptr;
            double v = std::strtod(start, &end);
            if (end == start) throw ValueError(format("malformed number at col %d", (int)col(i)));
            Token tk{TokenType::Number, v, "", col(i)};
            i += static_cast<std::size_t>(end - start);
            out.push_back(tk);
            continue;
        }
        if (std::isalpha(static_cast<unsigned char>(c)) != 0 || c == '_') {
            std::size_t j = i;
            while (j < n && (std::isalnum(static_cast<unsigned char>(s[j])) != 0 || s[j] == '_'))
                ++j;
            std::string id = s.substr(i, j - i);
            TokenType tt = TokenType::Ident;
            if (id == "let")
                tt = TokenType::Keyword_let;
            else if (id == "sum")
                tt = TokenType::Keyword_sum;
            out.push_back({tt, 0.0, id, col(i)});
            i = j;
            continue;
        }
        TokenType tt;
        switch (c) {
            case '+':
                tt = TokenType::Plus;
                break;
            case '-':
                tt = TokenType::Minus;
                break;
            case '*':
                tt = TokenType::Star;
                break;
            case '/':
                tt = TokenType::Slash;
                break;
            case '^':
                tt = TokenType::Caret;
                break;
            case '(':
                tt = TokenType::LParen;
                break;
            case ')':
                tt = TokenType::RParen;
                break;
            case '[':
                tt = TokenType::LBracket;
                break;
            case ']':
                tt = TokenType::RBracket;
                break;
            case ':':
                tt = TokenType::Colon;
                break;
            case ',':
                tt = TokenType::Comma;
                break;
            case '=':
                tt = TokenType::Equals;
                break;
            default:
                throw ValueError(format("illegal character '%c' at col %d", c, (int)col(i)));
        }
        out.push_back({tt, 0.0, "", col(i)});
        ++i;
    }
    out.push_back({TokenType::End, 0.0, "", col(i)});
    return out;
}

// ---------------------------------------------------------------------------
// Name-resolution helpers
// ---------------------------------------------------------------------------

// The single bucket of thermodynamic inputs: DSL spelling -> CoolProp::parameters
// key.  See the rationale for the curated allowlist in Expression.h.
static const std::vector<std::pair<std::string, parameters>>& inputTableImpl() {
    // State variables and pure-fluid metadata.  See below for what is deliberately
    // NOT here, and why.
    static const std::vector<std::pair<std::string, parameters>> table = {
      {"T", iT}, {"rhomolar", iDmolar}, {"rhomass", iDmass}, {"molar_mass", imolar_mass}, {"p", iP}};
    return table;
}

// Deliberately ABSENT from the table above: the critical point (T_critical,
// rhomolar_critical, p_critical) and the EOS reducing state (T_reducing,
// rhomolar_reducing).  Both were tried and removed.
//
//  * The critical-point accessors are not the fluid file's STATES.critical.  With
//    superancillaries enabled -- the DEFAULT -- calc_T_critical and
//    calc_rhomolar_critical return the NUMERICAL critical point from the
//    superancillary (the state satisfying dp/drho|T = d2p/drho2|T = 0).  A formula
//    reading them would give different numbers depending on
//    ENABLE_SUPERANCILLARIES, and fluid data must not depend on a runtime config
//    flag.  Xenon is the worked example: its 2021 viscosity paper says it adopts
//    Tc and rho_c from the Lemmon-Span EOS, Xenon.json carries exactly those
//    (289.733 K, 8400 mol/m^3), and reading them through the inputs still moved a
//    verification point by 7e-5.
//
//  * The reducing state is a fitting convenience of the equation of state, while
//    a correlation writing `T_reducing` means its OWN fitted value.  Exposing the
//    name would advise the wrong thing and silently steal it from every formula
//    already using it as a constant -- six golden formulas in the test suite do.
//
//  * Underlying both: a reducing parameter is a defining constant of the fit that
//    produced the correlation, quoted to finite precision in the source paper.
//    Reproducing the correlation means using THAT number, not a better estimate of
//    the same physical quantity.  Reducing parameters belong in the block's
//    `constants`, frozen.  The only legitimate use for a live critical point is
//    physics that references the critical point itself -- a crossover critical
//    enhancement -- which no shipped viscosity correlation implements.
//
// Adding a row to this table is a BREAKING change to the language: it steals that
// name from every formula already using it as a constant.  Weigh it that way.

static bool inputForName(const std::string& nm, parameters& out) {
    for (const auto& e : inputTableImpl())
        if (nm == e.first) {
            out = e.second;
            return true;
        }
    return false;
}
static bool funcForName(const std::string& nm, Func& out, int& arity) {
    struct E
    {
        const char* n;
        Func f;
        int a;
    };
    static const std::array<E, 12> table = {{{"exp", Func::exp, 1},
                                             {"ln", Func::ln, 1},
                                             {"log10", Func::log10, 1},
                                             {"sqrt", Func::sqrt, 1},
                                             {"abs", Func::abs, 1},
                                             {"pow", Func::pow, 2},
                                             {"sinh", Func::sinh, 1},
                                             {"cosh", Func::cosh, 1},
                                             {"tanh", Func::tanh, 1},
                                             {"sin", Func::sin, 1},
                                             {"cos", Func::cos, 1},
                                             {"atan", Func::atan, 1}}};
    for (const auto& e : table)
        if (nm == e.n) {
            out = e.f;
            arity = e.a;
            return true;
        }
    return false;
}

// ---------------------------------------------------------------------------
// Real ProgramData
// ---------------------------------------------------------------------------

struct ProgramData
{
    int numScalars = 0;
    std::vector<std::pair<int, double>> constantInits{};   // (slot, value)
    std::vector<std::pair<parameters, int>> inputSlots{};  // (key, slot)
    std::vector<std::pair<int, NodePtr>> lets{};           // (slot, node) in order
    NodePtr result{};
    std::vector<std::vector<double>> arrays{};  // by array slot
    std::vector<parameters> inputOrder{};       // cached for requiredInputs()
};

// ---------------------------------------------------------------------------
// Binder: rewrites parse-time nodes into runtime nodes; assigns slots
// ---------------------------------------------------------------------------

class Binder
{
   public:
    Binder(ProgramData& pd, const std::map<std::string, double>& consts, const std::map<std::string, std::vector<double>>& arrays)
      : d(pd), constants(consts), arraysIn(arrays) {}

    void run(ParseResult& pr) {
        for (auto& L : pr.lets) {
            NodePtr bound = bind(L.expr, /*indexName*/ "");
            int slot = newScalarSlot();
            letNames[L.name] = slot;  // available to later statements + result
            d.lets.emplace_back(slot, std::move(bound));
        }
        d.result = bind(pr.result, "");
        d.numScalars = nextScalar;
    }

   private:
    // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
    // non-owning views; valid for the bounded bind lifetime
    ProgramData& d;
    const std::map<std::string, double>& constants;
    const std::map<std::string, std::vector<double>>& arraysIn;
    // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
    int nextScalar = 0;
    std::map<std::string, int> letNames{};
    std::map<std::string, int> constSlots{};
    std::map<std::string, int> inputSlots{};
    std::map<std::string, int> arraySlots{};

    int newScalarSlot() {
        return nextScalar++;
    }

    int resolveScalar(const std::string& nm) {
        // Resolution order: let -> thermodynamic input -> JSON constant.
        // Inputs are checked before constants so that a JSON constant can never
        // shadow a state variable (T/rhomolar/rhomass/molar_mass/p).
        auto itL = letNames.find(nm);
        if (itL != letNames.end()) return itL->second;
        parameters key;
        if (inputForName(nm, key)) {
            auto s = inputSlots.find(nm);
            if (s != inputSlots.end()) return s->second;
            int slot = newScalarSlot();
            inputSlots[nm] = slot;
            d.inputSlots.emplace_back(key, slot);
            d.inputOrder.push_back(key);
            return slot;
        }
        auto itC = constants.find(nm);
        if (itC != constants.end()) {
            auto s = constSlots.find(nm);
            if (s != constSlots.end()) return s->second;
            int slot = newScalarSlot();
            constSlots[nm] = slot;
            d.constantInits.emplace_back(slot, itC->second);
            return slot;
        }
        throw ValueError(format("unknown variable '%s'", nm.c_str()));
    }

    int resolveArray(const std::string& nm) {
        auto it = arraysIn.find(nm);
        if (it == arraysIn.end()) throw ValueError(format("unknown array '%s'", nm.c_str()));
        auto s = arraySlots.find(nm);
        if (s != arraySlots.end()) return s->second;
        int slot = static_cast<int>(d.arrays.size());
        d.arrays.push_back(it->second);
        arraySlots[nm] = slot;
        return slot;
    }

    NodePtr bind(NodePtr& node, const std::string& indexName) {
        if (auto* num = dynamic_cast<NumNode*>(node.get())) {
            return std::make_unique<NumNode>(num->v);
        }
        if (auto* nm = dynamic_cast<NameNode*>(node.get())) {
            if (!indexName.empty() && nm->name == indexName) return std::make_unique<IndexNode>();
            return std::make_unique<ScalarNode>(resolveScalar(nm->name));
        }
        if (auto* ar = dynamic_cast<ArrayRefNode*>(node.get())) {
            if (indexName.empty() || ar->indexName != indexName)
                throw ValueError(format("array '%s' subscript must be the enclosing sum index", ar->arrayName.c_str()));
            return std::make_unique<ArrayNode>(resolveArray(ar->arrayName));
        }
        if (auto* neg = dynamic_cast<NegNode*>(node.get())) {
            return std::make_unique<NegNode>(bind(neg->child, indexName));
        }
        if (auto* bn = dynamic_cast<BinNode*>(node.get())) {
            NodePtr l = bind(bn->l, indexName);
            NodePtr r = bind(bn->r, indexName);
            return std::make_unique<BinNode>(bn->op, std::move(l), std::move(r));
        }
        if (auto* call = dynamic_cast<NamedCallNode*>(node.get())) {
            Func f;
            int arity;
            if (!funcForName(call->name, f, arity)) throw ValueError(format("unknown function '%s'", call->name.c_str()));
            if (static_cast<int>(call->args.size()) != arity)
                throw ValueError(format("function '%s' expects %d argument(s)", call->name.c_str(), arity));
            std::vector<NodePtr> bound;
            bound.reserve(call->args.size());
            for (auto& a : call->args)
                bound.push_back(bind(a, indexName));
            return std::make_unique<CallNode>(f, std::move(bound));
        }
        if (auto* sum = dynamic_cast<ParseSumNode*>(node.get())) {
            if (!indexName.empty()) throw ValueError("nested summation is not supported in v1");
            long n = sumLength(*sum, sum->indexName);
            NodePtr body = bind(sum->body, sum->indexName);
            return std::make_unique<SumNode>(n, std::move(body));
        }
        throw ValueError("internal: unhandled parse node in binder");
    }

    long sumLength(ParseSumNode& sum, const std::string& idx) {
        long len = -1;
        std::function<void(Node*)> scan = [&](Node* nd) {
            if (auto* ar = dynamic_cast<ArrayRefNode*>(nd)) {
                if (ar->indexName == idx) {
                    auto it = arraysIn.find(ar->arrayName);
                    if (it == arraysIn.end()) throw ValueError(format("unknown array '%s'", ar->arrayName.c_str()));
                    long n = static_cast<long>(it->second.size());
                    if (len == -1)
                        len = n;
                    else if (len != n)
                        throw ValueError(format("array '%s' length %ld != %ld in sum", ar->arrayName.c_str(), n, len));
                }
            } else if (auto* neg = dynamic_cast<NegNode*>(nd)) {
                scan(neg->child.get());
            } else if (auto* bn = dynamic_cast<BinNode*>(nd)) {
                scan(bn->l.get());
                scan(bn->r.get());
            } else if (auto* call = dynamic_cast<NamedCallNode*>(nd)) {
                for (auto& a : call->args)
                    scan(a.get());
            }
        };
        scan(sum.body.get());
        if (len <= 0) throw ValueError("sum body references no array subscripted by the index");
        return len;
    }
};

}  // namespace detail

// ---------------------------------------------------------------------------
// Program methods + compile (namespace CoolProp::expression)
// ---------------------------------------------------------------------------

const std::vector<std::pair<std::string, parameters>>& inputTable() {
    return detail::inputTableImpl();
}

std::string inputName(parameters key) {
    for (const auto& e : detail::inputTableImpl())
        if (e.second == key) return e.first;
    throw ValueError(format("parameter key %d is not an expression input", static_cast<int>(key)));
}

const detail::ProgramData& Program::data() const {
    // Program is default-constructible (it has to be, to sit inside the transport
    // structs), so a caller can reach the accessors before compile() ever ran.
    if (!m_data) throw ValueError("expression Program has no compiled body");
    return *m_data;
}

double Program::evaluate(const std::vector<double>& inputVals) const {
    const detail::ProgramData& d = data();
    if (inputVals.size() != d.inputOrder.size()) {
        throw ValueError(
          format("expression needs %d input value(s), got %d", static_cast<int>(d.inputOrder.size()), static_cast<int>(inputVals.size())));
    }
    std::vector<double> scalars(static_cast<std::size_t>(d.numScalars), 0.0);
    for (const auto& c : d.constantInits)
        scalars[c.first] = c.second;
    for (std::size_t k = 0; k < d.inputSlots.size(); ++k)
        scalars[d.inputSlots[k].second] = inputVals[k];
    detail::EvalState st{scalars, d.arrays, 0};
    for (const auto& L : d.lets)
        scalars[L.first] = L.second->eval(st);
    return d.result->eval(st);
}

const std::vector<parameters>& Program::requiredInputs() const {
    return data().inputOrder;
}

Program compile(const std::string& source, const std::map<std::string, double>& constants, const std::map<std::string, std::vector<double>>& arrays) {
    using namespace detail;
    // Inputs resolve before constants, so a constant named `T` or `p` would be
    // dead data and the formula would quietly mean something other than what its
    // author wrote.  Reject the collision at compile time rather than silently
    // discarding the constant.
    for (const auto& c : constants) {
        parameters key;
        if (inputForName(c.first, key)) {
            throw ValueError(format("constant '%s' collides with the thermodynamic input of the same name", c.first.c_str()));
        }
    }
    std::vector<Token> toks = lex(source);
    Parser parser(std::move(toks));
    ParseResult pr = parser.parse();
    auto pd = std::make_shared<ProgramData>();
    Binder binder(*pd, constants, arrays);
    binder.run(pr);
    Program prog;
    prog.m_data = pd;
    return prog;
}

}  // namespace expression
}  // namespace CoolProp
