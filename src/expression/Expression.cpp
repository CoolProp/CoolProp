#include "CoolProp/expression/Expression.h"
#include "CoolProp/expression/detail/Lexer.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/strings.h"
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <memory>
#include <vector>

namespace CoolProp {
namespace expression {
namespace detail {

// Function intrinsics supported by the grammar.
enum class Func : std::uint8_t { exp, ln, log10, sqrt, abs, pow, sinh, cosh, tanh, sin, cos, atan };

// Evaluation context: scalars indexed by slot, arrays indexed by slot, and the
// single active summation index `i` (v1 forbids nested sums, so one is enough).
struct EvalState
{
    std::vector<double>& scalars;
    const std::vector<std::vector<double>>& arrays;
    long i = 0;
};

struct Node
{
    virtual ~Node() = default;
    virtual double eval(EvalState& st) const = 0;
};
using NodePtr = std::unique_ptr<Node>;

struct NumNode : Node
{
    double v;
    explicit NumNode(double value) : v(value) {}
    double eval(EvalState&) const override { return v; }
};

struct ScalarNode : Node  // intrinsic, constant, derived, or let — all read a slot
{
    int slot;
    explicit ScalarNode(int s) : slot(s) {}
    double eval(EvalState& st) const override { return st.scalars[slot]; }
};

struct IndexNode : Node  // bare summation index used as a value
{
    double eval(EvalState& st) const override { return static_cast<double>(st.i); }
};

struct ArrayNode : Node  // arr[i]
{
    int slot;
    explicit ArrayNode(int s) : slot(s) {}
    double eval(EvalState& st) const override {
        return st.arrays[slot][static_cast<std::size_t>(st.i)];
    }
};

struct NegNode : Node
{
    NodePtr child;
    explicit NegNode(NodePtr c) : child(std::move(c)) {}
    double eval(EvalState& st) const override { return -child->eval(st); }
};

struct BinNode : Node
{
    char op;  // '+','-','*','/','^'
    NodePtr l, r;
    BinNode(char o, NodePtr a, NodePtr b) : op(o), l(std::move(a)), r(std::move(b)) {}
    double eval(EvalState& st) const override {
        double a = l->eval(st), b = r->eval(st);
        switch (op) {
            case '+': return a + b;
            case '-': return a - b;
            case '*': return a * b;
            case '/': return a / b;
            case '^': return std::pow(a, b);
        }
        return 0.0;  // unreachable; op validated at parse time
    }
};

struct CallNode : Node
{
    Func fn;
    std::vector<NodePtr> args;
    CallNode(Func f, std::vector<NodePtr> a) : fn(f), args(std::move(a)) {}
    double eval(EvalState& st) const override {
        double x = args[0]->eval(st);
        switch (fn) {
            case Func::exp: return std::exp(x);
            case Func::ln: return std::log(x);
            case Func::log10: return std::log10(x);
            case Func::sqrt: return std::sqrt(x);
            case Func::abs: return std::fabs(x);
            case Func::pow: return std::pow(x, args[1]->eval(st));
            case Func::sinh: return std::sinh(x);
            case Func::cosh: return std::cosh(x);
            case Func::tanh: return std::tanh(x);
            case Func::sin: return std::sin(x);
            case Func::cos: return std::cos(x);
            case Func::atan: return std::atan(x);
        }
        return 0.0;
    }
};

struct SumNode : Node
{
    long n;  // fixed array length, resolved at bind time
    NodePtr body;
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
    std::string name;
    explicit NameNode(std::string nm) : name(std::move(nm)) {}
    double eval(EvalState&) const override {
        throw ValueError("internal: unbound NameNode evaluated");
    }
};

struct ArrayRefNode : Node  // name[index]; index identifier captured for validation
{
    std::string arrayName, indexName;
    ArrayRefNode(std::string a, std::string idx) : arrayName(std::move(a)), indexName(std::move(idx)) {}
    double eval(EvalState&) const override {
        throw ValueError("internal: unbound ArrayRefNode evaluated");
    }
};

struct ParseSumNode : Node  // sum(indexName: body); length resolved at bind
{
    std::string indexName;
    NodePtr body;
    double eval(EvalState&) const override {
        throw ValueError("internal: unbound ParseSumNode evaluated");
    }
};

struct NamedCallNode : Node  // funcName(args...); resolved to Func at bind
{
    std::string name;
    std::vector<NodePtr> args;
    NamedCallNode(std::string nm, std::vector<NodePtr> a) : name(std::move(nm)), args(std::move(a)) {}
    double eval(EvalState&) const override {
        throw ValueError("internal: unbound NamedCallNode evaluated");
    }
};

struct LetStmt
{
    std::string name;
    NodePtr expr;
    std::size_t col;
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
    std::vector<Token> t;
    std::size_t pos = 0;

    const Token& peek() const { return t[pos]; }
    const Token& advance() { return t[pos++]; }
    bool match(TokenType tt) {
        if (t[pos].type == tt) {
            ++pos;
            return true;
        }
        return false;
    }
    void skipSeps() {
        while (t[pos].type == TokenType::StmtSep) ++pos;
    }
    const Token& expect(TokenType tt, const char* what) {
        if (t[pos].type != tt)
            throw ValueError(format("expected %s at col %d", what, (int)t[pos].col));
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

    int lbp(TokenType tt) const {
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

    char opChar(TokenType tt) const {
        switch (tt) {
            case TokenType::Plus: return '+';
            case TokenType::Minus: return '-';
            case TokenType::Star: return '*';
            case TokenType::Slash: return '/';
            case TokenType::Caret: return '^';
            default: return '?';
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
        if (match(TokenType::Minus)) return std::make_unique<NegNode>(parseExpr(40));  // unary binds tight
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
            while (match(TokenType::Comma)) args.push_back(parseExpr(0));
        }
        expect(TokenType::RParen, "')' to close call");
        return std::make_unique<NamedCallNode>(name, std::move(args));
    }
};

// ---------------------------------------------------------------------------

struct ProgramData
{
    std::vector<Intrinsic> intrinsics;
    std::vector<Derived> deriveds;
};

std::vector<Token> lex(const std::string& s) {
    std::vector<Token> out;
    std::size_t i = 0, n = s.size();
    auto col = [&](std::size_t p) { return p + 1; };
    while (i < n) {
        char c = s[i];
        if (c == ' ' || c == '\t' || c == '\r') { ++i; continue; }
        if (c == '\n' || c == ';') {
            out.push_back({TokenType::StmtSep, 0.0, "", col(i)});
            ++i;
            continue;
        }
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            const char* start = s.c_str() + i;
            char* end = nullptr;
            double v = std::strtod(start, &end);
            if (end == start) throw ValueError(format("malformed number at col %d", (int)col(i)));
            Token tk{TokenType::Number, v, "", col(i)};
            i += static_cast<std::size_t>(end - start);
            out.push_back(tk);
            continue;
        }
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            std::size_t j = i;
            while (j < n && (std::isalnum(static_cast<unsigned char>(s[j])) || s[j] == '_')) ++j;
            std::string id = s.substr(i, j - i);
            TokenType tt = TokenType::Ident;
            if (id == "let") tt = TokenType::Keyword_let;
            else if (id == "sum") tt = TokenType::Keyword_sum;
            out.push_back({tt, 0.0, id, col(i)});
            i = j;
            continue;
        }
        TokenType tt;
        switch (c) {
            case '+': tt = TokenType::Plus; break;
            case '-': tt = TokenType::Minus; break;
            case '*': tt = TokenType::Star; break;
            case '/': tt = TokenType::Slash; break;
            case '^': tt = TokenType::Caret; break;
            case '(': tt = TokenType::LParen; break;
            case ')': tt = TokenType::RParen; break;
            case '[': tt = TokenType::LBracket; break;
            case ']': tt = TokenType::RBracket; break;
            case ':': tt = TokenType::Colon; break;
            case ',': tt = TokenType::Comma; break;
            case '=': tt = TokenType::Equals; break;
            default: throw ValueError(format("illegal character '%c' at col %d", c, (int)col(i)));
        }
        out.push_back({tt, 0.0, "", col(i)});
        ++i;
    }
    out.push_back({TokenType::End, 0.0, "", col(i)});
    return out;
}

}  // namespace detail

double Program::evaluate(const double*, const double*) const {
    // Stub: real implementation lands in Task 5.
    return 7.0;
}
const std::vector<Intrinsic>& Program::requiredIntrinsics() const {
    return m_data->intrinsics;
}
const std::vector<Derived>& Program::requiredDerived() const {
    return m_data->deriveds;
}

Program compile(const std::string&, const std::map<std::string, double>&,
                const std::map<std::string, std::vector<double>>&) {
    Program p;
    p.m_data = std::make_shared<detail::ProgramData>();
    return p;
}

}  // namespace expression
}  // namespace CoolProp
