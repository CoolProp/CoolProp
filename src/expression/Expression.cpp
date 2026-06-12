#include "CoolProp/expression/Expression.h"
#include "CoolProp/expression/detail/Lexer.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/strings.h"
#include <cctype>
#include <cstdlib>

namespace CoolProp {
namespace expression {
namespace detail {
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
