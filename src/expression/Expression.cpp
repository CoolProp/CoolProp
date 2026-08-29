#include "CoolProp/expression/Expression.h"
#include "CoolProp/expression/detail/Lexer.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/DataStructures.h"
#include "CoolProp/detail/strings.h"
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <algorithm>
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

// Thermodynamic names are NOT resolved from a list kept here.  They are resolved by
// CoolProp::is_valid_parameter(), so the DSL's vocabulary is exactly CoolProp's and a
// new quantity costs no code at all -- krypton's entropy-scaling correlation needed
// Smolar_residual, Bvirial and dBvirial_dT, and under the old curated table each was
// a row plus a recompile for something keyed_output() already knew how to produce.
//
// What replaces the table is not "resolve everything", but two narrower rules.
//
// FIRST, resolution is opt-in per formula.  A block declares `state_variables`, and
// only those names are state inside it.  The old design had one namespace shared by
// state and by the author's own constants and arrays, so every name the DSL knew was
// reserved everywhere -- and the DSL spelled pressure with an invented lowercase `p`,
// which is what viscosity papers universally call their exponent array.  Propylene
// glycol's dilute term had to rename p_i to np_i for a collision with a quantity it
// never referenced.  Under opt-in that block declares ["T"] and keeps `p`.
//
// SECOND, two classes stay refused even when declared, because opt-in expresses the
// author's intent and these are cases where the intent cannot be honoured:
// Spellings the DSL used to accept, mapped to what replaced them.  This exists ONLY
// to produce a good error: fluid JSON is data that third parties ship, and a block
// written against the previous format must fail at load time saying WHAT changed,
// not "unknown variable 'rhomolar'".  These names were the DSL's own inventions --
// CoolProp always called them Dmolar/Dmass/P -- and the invented lowercase `p` is
// what once collided with the exponent array every viscosity paper writes as p_i.
static const char* retiredSpelling(const std::string& nm) {
    if (nm == "rhomolar") return "Dmolar";
    if (nm == "rhomass") return "Dmass";
    if (nm == "p") return "P";
    return nullptr;
}

static bool deniedStateVariable(parameters key, std::string& why) {
    switch (key) {
        // keyed_output() for a transport output re-enters the correlation being
        // defined.  A viscosity formula that reads `V` is asking for itself.
        case iviscosity:
        case iconductivity:
        case iPrandtl:
        case isurface_tension:
            why = "it is a transport output, so reading it re-enters the correlation being defined";
            return true;
        // calc_T_critical()/calc_rhomolar_critical() return the NUMERICAL critical
        // point when ENABLE_SUPERANCILLARIES is set (the default), not STATES.critical.
        // A correlation reducing on them changes answer with configuration, which is
        // how xenon once missed its reference values by 7e-5.  The fix is to freeze
        // the number the paper's authors actually regressed against, as a constant.
        case iT_critical:
        case iP_critical:
        case irhomolar_critical:
        case irhomass_critical:
            why = "the critical point is configuration-dependent (ENABLE_SUPERANCILLARIES); "
                  "freeze the paper's value as a constant instead";
            return true;
        // The EOS reducing state is refused for a DIFFERENT reason, and saying so
        // matters: superancillaries do not touch get_reducing_state().  A correlation
        // that writes `T_reducing` means its OWN fitted reducing parameter, which is
        // in general not the EOS's and moves whenever the EOS section is revised.
        //
        // Tau and Delta are the same quantities wearing a different name --
        // keyed_output() computes them as _reducing.T/_T and _rhomolar/_reducing.rhomolar
        // -- so refusing the reducing state while admitting these would be a guard
        // with a door next to it.  They are also composition-dependent for mixtures.
        case iT_reducing:
        case iP_reducing:
        case irhomolar_reducing:
        case irhomass_reducing:
        case iTau:
        case iDelta:
            why = "the EOS reducing state is not the correlation's own; freeze the paper's "
                  "reducing parameters as constants instead (Tau and Delta are that state too)";
            return true;
        // Not real numbers.  Phase is an enum ordinal widened to double (and depends
        // on any imposed phase); Q is -1 outside the dome, a sentinel that would enter
        // arithmetic as data.  Both would compile into a plausible-looking formula.
        case iPhase:
        case iQ:
        case iQmass:
            why = "it is a phase index or a sentinel-valued quality, not a continuous state function";
            return true;
        // Fluid metadata and incompressible-only accessors: not functions of the
        // current state at all.  Several are unimplemented for HEOS and throw from
        // inside the evaluation, where fluid-file data has no business failing.
        case ifraction_min:
        case ifraction_max:
        case iT_freeze:
        case iGWP20:
        case iGWP100:
        case iGWP500:
        case iFH:
        case iHH:
        case iPH:
        case iODP:
            why = "it is fluid metadata, not a function of the current state";
            return true;
        default:
            return false;
    }
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
    std::vector<std::string> inputNames{};      // same order, as the author spelled them
};

// ---------------------------------------------------------------------------
// Binder: rewrites parse-time nodes into runtime nodes; assigns slots
// ---------------------------------------------------------------------------

class Binder
{
   public:
    Binder(ProgramData& pd, const std::map<std::string, double>& consts, const std::map<std::string, std::vector<double>>& arrays,
           const std::map<std::string, parameters>& declared)
      : d(pd), constants(consts), arraysIn(arrays), declaredState(declared) {}

    void run(ParseResult& pr) {
        for (auto& L : pr.lets) {
            // A `let` resolves before a declared state variable, so `let T = 500`
            // would quietly detach the formula from the state -- and only partly, in
            // a formula that also read `T` on an earlier line.  Declaring a name and
            // then rebinding it is never what an author means.
            if (declaredState.count(L.name))
                throw ValueError(format("let '%s' shadows the declared state variable of the same name", L.name.c_str()));
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
    const std::map<std::string, parameters>& declaredState;
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
        // Resolution order: let -> DECLARED state variable -> JSON constant.  A name
        // the block did not declare is never state, however well CoolProp knows it.
        auto itL = letNames.find(nm);
        if (itL != letNames.end()) return itL->second;
        auto itS = declaredState.find(nm);
        if (itS != declaredState.end()) {
            const parameters key = itS->second;
            auto s = inputSlots.find(nm);
            if (s != inputSlots.end()) return s->second;
            int slot = newScalarSlot();
            inputSlots[nm] = slot;
            d.inputSlots.emplace_back(key, slot);
            d.inputOrder.push_back(key);
            d.inputNames.push_back(nm);
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
        // A name CoolProp knows, but that this block did not declare, is the most
        // likely authoring mistake -- say so instead of "unknown variable".
        if (const char* now = retiredSpelling(nm)) {
            throw ValueError(format("unknown variable '%s' -- it was this DSL's own spelling for a thermodynamic "
                                    "quantity before blocks declared their inputs; CoolProp calls it '%s', and it must "
                                    "be listed in this block's \"state_variables\"",
                                    nm.c_str(), now));
        }
        parameters probe;
        std::string ignored;
        if (resolveStateVariable(nm, probe, ignored)) {
            throw ValueError(format("unknown variable '%s' -- it is a thermodynamic quantity; add it to this block's "
                                    "\"state_variables\" to read it from the state",
                                    nm.c_str()));
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

bool resolveStateVariable(const std::string& name, parameters& key, std::string& reason) {
    key = iundefined_parameter;  // both out-params written on every path
    reason.clear();
    parameters k;
    if (!is_valid_parameter(name, k)) {
        if (const char* now = detail::retiredSpelling(name)) {
            reason = format("'%s' was this DSL's own spelling; CoolProp calls it '%s'", name.c_str(), now);
        } else {
            reason = format("'%s' is not a CoolProp parameter name", name.c_str());
        }
        return false;
    }
    std::string why;
    if (detail::deniedStateVariable(k, why)) {
        reason = format("'%s' may not be used as a state variable: %s", name.c_str(), why.c_str());
        return false;
    }

    // is_valid_parameter() also accepts CoolProp's back-compat aliases and an
    // upper-cased spelling of every name: `A` is the speed of sound, `D` is Dmass,
    // `M` is molar_mass, `TAU` is Tau.  Those are a trap here in a way they are not
    // in PropsSI, because a DSL formula is FULL OF single-letter coefficient names.
    // An author who forgets to declare a paper coefficient `A` in `constants` gets
    // told to add it to state_variables, and doing so silently binds the speed of
    // sound.  Only the canonical spelling is accepted, so the DSL's vocabulary is
    // CoolProp's with exactly one name per quantity -- and still no list kept here.
    const std::string canonical = get_parameter_information(static_cast<int>(k), "short");
    if (canonical != name) {
        reason = format("'%s' is an alias; use CoolProp's canonical name '%s'", name.c_str(), canonical.c_str());
        return false;
    }
    key = k;
    return true;
}

std::string inputName(parameters key) {
    return get_parameter_information(static_cast<int>(key), "short");
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

Program compile(const std::string& source, const std::map<std::string, double>& constants, const std::map<std::string, std::vector<double>>& arrays,
                const std::vector<std::string>& state_variables) {
    using namespace detail;
    // Resolve the declaration first, so that a bad `state_variables` entry is
    // reported as such rather than surfacing later as "unknown variable".
    std::map<std::string, parameters> declared;
    for (const auto& nm : state_variables) {
        parameters key;
        std::string reason;
        if (!resolveStateVariable(nm, key, reason)) throw ValueError(format("state_variables: %s", reason.c_str()));
        if (!declared.emplace(nm, key).second) throw ValueError(format("state_variables lists '%s' twice", nm.c_str()));
        // Dedupe on the resolved key, not the spelling: two names for one quantity
        // would take two slots and cost two keyed_output() calls per evaluation.
        for (const auto& prev : declared) {
            if (prev.first != nm && prev.second == key)
                throw ValueError(format("state_variables lists '%s' and '%s', which are the same quantity", prev.first.c_str(), nm.c_str()));
        }
        // A declared name shadowing the author's own data would make that data dead
        // and quietly change what the formula means.  This is the collision the old
        // single-namespace design hit from the other direction, where the language
        // reserved the name whether the block wanted it or not.
        if (constants.count(nm)) throw ValueError(format("state variable '%s' collides with the constant of the same name", nm.c_str()));
        if (arrays.count(nm)) throw ValueError(format("state variable '%s' collides with the array of the same name", nm.c_str()));
    }
    std::vector<Token> toks = lex(source);
    Parser parser(std::move(toks));
    ParseResult pr = parser.parse();
    auto pd = std::make_shared<ProgramData>();
    Binder binder(*pd, constants, arrays, declared);
    binder.run(pr);
    // Declaring a state variable the formula never reads costs a keyed_output() call
    // per evaluation and misstates the block's dependencies, so it is an error rather
    // than dead weight.
    for (const auto& nm : state_variables) {
        if (std::find(pd->inputNames.begin(), pd->inputNames.end(), nm) == pd->inputNames.end())
            throw ValueError(format("state variable '%s' is declared but never used in the formula", nm.c_str()));
    }
    Program prog;
    prog.m_data = pd;
    return prog;
}

}  // namespace expression
}  // namespace CoolProp
