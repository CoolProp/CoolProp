#ifndef COOLPROP_EXPRESSION_LEXER_H
#define COOLPROP_EXPRESSION_LEXER_H

#include <cstdint>
#include <string>
#include <vector>

namespace CoolProp {
namespace expression {
namespace detail {

enum class TokenType : std::uint8_t {
    Number, Ident, Keyword_let, Keyword_sum,
    Plus, Minus, Star, Slash, Caret,
    LParen, RParen, LBracket, RBracket,
    Colon, Comma, Equals, StmtSep, End
};

struct Token
{
    TokenType type;
    double number = 0.0;     // valid when type == Number
    std::string text;        // valid when type == Ident
    std::size_t col = 0;     // 1-based column for error messages
};

/// Tokenize. Throws CoolProp::ValueError on an illegal character.
std::vector<Token> lex(const std::string& source);

}  // namespace detail
}  // namespace expression
}  // namespace CoolProp

#endif
