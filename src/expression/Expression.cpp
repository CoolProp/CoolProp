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
