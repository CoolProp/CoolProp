# Transport Property Expression DSL Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a hand-rolled, runtime-loadable expression DSL so Tier-A viscosity/conductivity correlation forms can be authored as formula strings in fluid JSON (`"type": "expression"`) instead of hardcoded C++.

**Architecture:** A self-contained `CoolProp::expression` module compiles a formula string once (lex → parse → bind) into an immutable `Program` (a tree-walking AST over a slot-based context). The evaluator is EOS-free; a thin host wrapper `ExpressionCorrelation` fetches the required intrinsic state (`T, rhomolar, rhomass, molar_mass`) and registered derived quantities (`p`) from `HelmholtzEOSMixtureBackend` and calls `Program::evaluate`. New enum members + JSON parse branches + dispatch cases wire it into the existing transport machinery additively; no existing form changes.

**Tech Stack:** C++17, no third-party dependency (coexists with nlohmann/json), Catch2 (`[expression]` tag), CMake GLOB build.

**Reference spec:** `docs/superpowers/specs/2026-06-12-transport-expression-dsl-design.md`
**Tracking issue:** CoolProp-onei

---

## File Structure

**New module (EOS-free core):**
- `include/CoolProp/expression/Expression.h` — public API: `Intrinsic`/`Derived` enums, `Program` class, free function `compile()`.
- `include/CoolProp/expression/detail/Lexer.h` — `Token`/`TokenType` + `lex()` declaration (exposed for unit testing).
- `src/expression/Expression.cpp` — Lexer, AST node structs, Pratt parser, binder, `Program::evaluate`. All correlation-internal.

**New host-integration layer (knows about HEOS):**
- `include/CoolProp/expression/ExpressionCorrelation.h` — `ExpressionCorrelation` (wraps a `Program`, owns `eval(HelmholtzEOSMixtureBackend&)`).
- `src/expression/ExpressionCorrelation.cpp` — fetches intrinsics/derived from HEOS.

**Modified (integration, exact anchors from recon):**
- `include/CoolProp/CoolPropFluid.h` — 4 enums + 4 container structs get `*_EXPRESSION` member + `ExpressionData expression_data`.
- `src/Backends/Helmholtz/Fluids/FluidLibrary.h` — 4 parse functions get an `else if ("expression")` branch.
- `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp` — 5 dispatch switches get a `case *_EXPRESSION`.
- `CMakeLists.txt` — add `src/expression/*.cpp` to `APP_SOURCES`.

**Tests:**
- `src/Tests/CoolProp-Tests.cpp` — `[expression]`-tagged unit + golden-regression TEST_CASEs (guarded by `#if defined(ENABLE_CATCH)`).

---

## Conventions for every task

- Build the test runner: `cmake --build build_catch --target CatchTestRunner -j8`
  (configure once: `cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release`)
- Run a tag: `./build_catch/CatchTestRunner "[expression]"`
- `CoolPropDbl` is `double` in this build; the DSL core uses plain `double` and only the host layer uses `CoolPropDbl`/HEOS.
- Errors at compile time throw `CoolProp::ValueError` (from `CoolProp/Exceptions.h`), message built with `format(...)` (from `CoolProp/detail/strings.h`).
- `.beads/issues.jsonl` must not enter source commits: before each commit run
  `git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null` and commit with `--no-verify`.

---

## Task 1: Scaffold the module, wire CMake, smoke test

**Files:**
- Create: `include/CoolProp/expression/Expression.h`
- Create: `src/expression/Expression.cpp`
- Modify: `CMakeLists.txt:313-333` (append expression sources)
- Test: `src/Tests/CoolProp-Tests.cpp` (append a smoke TEST_CASE)

- [ ] **Step 1: Write the failing test** — append to `src/Tests/CoolProp-Tests.cpp` just before the final `#endif`:

```cpp
#if defined(ENABLE_CATCH)
#    include "CoolProp/expression/Expression.h"

TEST_CASE("expression module links and evaluates a constant", "[expression]") {
    using namespace CoolProp::expression;
    Program prog = compile("3 + 4", {}, {});
    CHECK(prog.evaluate(nullptr, nullptr) == Approx(7.0));
}
#endif
```

- [ ] **Step 2: Create the public header** `include/CoolProp/expression/Expression.h`:

```cpp
#ifndef COOLPROP_EXPRESSION_H
#define COOLPROP_EXPRESSION_H

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace CoolProp {
namespace expression {

/// EOS-free independent state and pure fluid metadata, filled by the host per eval.
enum class Intrinsic { T, rhomolar, rhomass, molar_mass };
/// State-dependent quantities the EOS must compute; v1 registry holds only p.
enum class Derived { p };

namespace detail {
struct ProgramData;
}

/// An immutable, compiled expression. Cheap to copy (shared, refcounted body).
class Program
{
   public:
    /// Evaluate. intrinsicVals/derivedVals are arrays in the order given by
    /// requiredIntrinsics()/requiredDerived(); pass nullptr when none are required.
    double evaluate(const double* intrinsicVals, const double* derivedVals) const;
    /// Intrinsics this program references, in the order evaluate() expects them.
    const std::vector<Intrinsic>& requiredIntrinsics() const;
    /// Derived quantities this program references, in the order evaluate() expects them.
    const std::vector<Derived>& requiredDerived() const;

   private:
    friend Program compile(const std::string&, const std::map<std::string, double>&,
                           const std::map<std::string, std::vector<double>>&);
    std::shared_ptr<const detail::ProgramData> m_data;
};

/// Compile a formula string. `constants` are scalar names -> SI values; `arrays`
/// are vector names -> values. Throws CoolProp::ValueError on any lex/parse/bind error.
Program compile(const std::string& source, const std::map<std::string, double>& constants,
                const std::map<std::string, std::vector<double>>& arrays);

}  // namespace expression
}  // namespace CoolProp

#endif
```

- [ ] **Step 3: Create a minimal stub** `src/expression/Expression.cpp` (replaced in later tasks):

```cpp
#include "CoolProp/expression/Expression.h"
#include "CoolProp/Exceptions.h"
#include "CoolProp/detail/strings.h"

namespace CoolProp {
namespace expression {
namespace detail {
struct ProgramData
{
    std::vector<Intrinsic> intrinsics;
    std::vector<Derived> deriveds;
};
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
```

- [ ] **Step 4: Wire CMake** — in `CMakeLists.txt`, immediately after line 333 (the `list(APPEND APP_SOURCES ${SVD_SOURCES} ...)` block region, right after the `foreach(backend ...)` loop that ends near line 333), add:

```cmake
# Expression DSL for runtime-loaded transport correlations
file(GLOB EXPRESSION_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/expression/*.cpp")
list(APPEND APP_SOURCES ${EXPRESSION_SOURCES})
```

- [ ] **Step 5: Build and run** — Expected: builds; test passes.

```bash
cmake -B build_catch -S . -DCOOLPROP_CATCH_MODULE=ON -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_catch --target CatchTestRunner -j8
./build_catch/CatchTestRunner "[expression]"
```
Expected: `All tests passed (1 assertion in 1 test case)`.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add include/CoolProp/expression/Expression.h src/expression/Expression.cpp CMakeLists.txt src/Tests/CoolProp-Tests.cpp
git commit --no-verify -m "feat(expression): scaffold DSL module + CMake wiring (CoolProp-onei)"
```

---

## Task 2: Lexer

**Files:**
- Create: `include/CoolProp/expression/detail/Lexer.h`
- Modify: `src/expression/Expression.cpp` (implement `lex`)
- Test: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Write the failing test** — append inside the existing `#if defined(ENABLE_CATCH)` expression region:

```cpp
#    include "CoolProp/expression/detail/Lexer.h"

TEST_CASE("lexer tokenizes numbers, idents, operators", "[expression]") {
    using namespace CoolProp::expression::detail;
    std::vector<Token> t = lex("let x = a[i]*T^2.66958e-08 + sum(i: b)");
    // Spot-check a few token types/values
    REQUIRE(t.front().type == TokenType::Keyword_let);
    bool sawScientific = false, sawCaret = false, sawColon = false, sawLBracket = false;
    for (const auto& tok : t) {
        if (tok.type == TokenType::Number && tok.number == Approx(2.66958e-08)) sawScientific = true;
        if (tok.type == TokenType::Caret) sawCaret = true;
        if (tok.type == TokenType::Colon) sawColon = true;
        if (tok.type == TokenType::LBracket) sawLBracket = true;
    }
    CHECK(sawScientific);
    CHECK(sawCaret);
    CHECK(sawColon);
    CHECK(sawLBracket);
    CHECK(t.back().type == TokenType::End);
}

TEST_CASE("lexer rejects an illegal character", "[expression]") {
    using namespace CoolProp::expression::detail;
    CHECK_THROWS_AS(lex("a @ b"), CoolProp::ValueError);
}
```

- [ ] **Step 2: Run test to verify it fails** — `./build_catch/CatchTestRunner "[expression]"` → FAIL (no `Lexer.h`).

- [ ] **Step 3: Create** `include/CoolProp/expression/detail/Lexer.h`:

```cpp
#ifndef COOLPROP_EXPRESSION_LEXER_H
#define COOLPROP_EXPRESSION_LEXER_H

#include <string>
#include <vector>

namespace CoolProp {
namespace expression {
namespace detail {

enum class TokenType {
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
```

- [ ] **Step 4: Implement `lex`** — in `src/expression/Expression.cpp`, add the include and the function (place near the top, inside `namespace CoolProp { namespace expression { namespace detail {`):

```cpp
#include "CoolProp/expression/detail/Lexer.h"
#include <cctype>
#include <cstdlib>

// ... inside namespace CoolProp::expression::detail ...

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
```

- [ ] **Step 5: Build & run** — `./build_catch/CatchTestRunner "[expression]"` → PASS.

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add include/CoolProp/expression/detail/Lexer.h src/expression/Expression.cpp src/Tests/CoolProp-Tests.cpp
git commit --no-verify -m "feat(expression): lexer with scientific-notation numbers (CoolProp-onei)"
```

---

## Task 3: AST nodes + evaluator (slot-based, hand-built trees)

**Files:**
- Modify: `src/expression/Expression.cpp` (add AST node structs + `EvalState`)
- Test: covered indirectly in Task 5 (compile + evaluate). No standalone test here — the nodes are internal. This task is pure implementation that Task 5's tests exercise.

- [ ] **Step 1: Add the node hierarchy and eval machinery** — in `src/expression/Expression.cpp`, inside `namespace CoolProp::expression::detail`, above `lex` is fine; types must precede `ProgramData`:

```cpp
#include <cmath>
#include <functional>

// Function intrinsics supported by the grammar.
enum class Func { exp, ln, log10, sqrt, abs, pow, sinh, cosh, tanh, sin, cos, atan };

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
```

- [ ] **Step 2: Build only (no test yet)** — `cmake --build build_catch --target CatchTestRunner -j8`. Expected: compiles (unused warnings acceptable). If it builds, the node layer is syntactically integrated.

- [ ] **Step 3: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/expression/Expression.cpp
git commit --no-verify -m "feat(expression): slot-based AST node hierarchy + evaluator (CoolProp-onei)"
```

---

## Task 4: Parser (tokens → name-referencing AST)

The parser builds nodes that still reference names; it stores names in side tables so the binder (Task 5) can resolve slots. To keep one node set, the parser emits placeholder nodes carrying names, captured via small parse-only structs.

**Files:**
- Modify: `src/expression/Expression.cpp` (add parser producing an intermediate `ParseResult`)
- Test: exercised through Task 5's `compile()` tests.

- [ ] **Step 1: Add parse-time representation and the Pratt parser** — in `src/expression/Expression.cpp`, inside the `detail` namespace, after the node definitions:

```cpp
// Parse-time AST: same shape, but identifiers are unresolved (by name).
// We reuse runtime nodes where they carry no names (Num, Neg, Bin, Call, Sum)
// and use these two name-carrying nodes that the binder rewrites into
// ScalarNode / ArrayNode once slots are known.
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

struct LetStmt { std::string name; NodePtr expr; std::size_t col; };
struct ParseResult { std::vector<LetStmt> lets; NodePtr result; };

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
    bool match(TokenType tt) { if (t[pos].type == tt) { ++pos; return true; } return false; }
    void skipSeps() { while (t[pos].type == TokenType::StmtSep) ++pos; }
    const Token& expect(TokenType tt, const char* what) {
        if (t[pos].type != tt) throw ValueError(format("expected %s at col %d", what, (int)t[pos].col));
        return t[pos++];
    }

    LetStmt parseLet() {
        std::size_t col = advance().col;  // consume 'let'
        const Token& id = expect(TokenType::Ident, "identifier after 'let'");
        expect(TokenType::Equals, "'=' in let binding");
        NodePtr e = parseExpr(0);
        if (peek().type != TokenType::StmtSep && peek().type != TokenType::End)
            throw ValueError(format("expected end of let statement at col %d", (int)peek().col));
        return LetStmt{id.text, std::move(e), col};
    }

    // Precedence-climbing. lbp: + - =10, * / =20, ^ =30 (right-assoc).
    int lbp(TokenType tt) const {
        switch (tt) {
            case TokenType::Plus: case TokenType::Minus: return 10;
            case TokenType::Star: case TokenType::Slash: return 20;
            case TokenType::Caret: return 30;
            default: return 0;
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

    NodePtr parsePrefix() {
        if (match(TokenType::Minus)) return std::make_unique<NegNode>(parseExpr(40));  // unary binds tight
        return parsePrimary();
    }

    NodePtr parsePrimary() {
        const Token& tk = peek();
        if (tk.type == TokenType::Number) { advance(); return std::make_unique<NumNode>(tk.number); }
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
                expect(TokenType::RBracket, "']'");
                return std::make_unique<ArrayRefNode>(name, idx.text);
            }
            return std::make_unique<NameNode>(name);
        }
        throw ValueError(format("unexpected token at col %d", (int)tk.col));
    }

    NodePtr parseSum() {
        advance();  // 'sum'
        expect(TokenType::LParen, "'(' after sum");
        const Token& idx = expect(TokenType::Ident, "index identifier in sum");
        expect(TokenType::Colon, "':' in sum");
        NodePtr body = parseExpr(0);
        expect(TokenType::RParen, "')' to close sum");
        // Length resolved at bind; index name stashed via a temporary SumNode wrapper
        // carrying the name in a NameNode sentinel is avoided — instead we wrap with
        // a ParseSumNode that the binder consumes.
        auto s = std::make_unique<ParseSumNode>();
        s->indexName = idx.text;
        s->body = std::move(body);
        return s;
    }

   public:
    // Parse-only sum wrapper (carries index name; binder converts to SumNode).
    struct ParseSumNode : Node
    {
        std::string indexName;
        NodePtr body;
        double eval(EvalState&) const override {
            throw ValueError("internal: unbound ParseSumNode evaluated");
        }
    };

   private:
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

   public:
    // Parse-only call wrapper (carries function name; binder resolves to Func).
    struct NamedCallNode : Node
    {
        std::string name;
        std::vector<NodePtr> args;
        NamedCallNode(std::string nm, std::vector<NodePtr> a) : name(std::move(nm)), args(std::move(a)) {}
        double eval(EvalState&) const override {
            throw ValueError("internal: unbound NamedCallNode evaluated");
        }
    };
};
```

> Note: `ParseSumNode` and `NamedCallNode` are declared inside `Parser` for locality but referenced by the binder; if the compiler complains about nested-class access, hoist both to namespace scope (above `Parser`). Prefer namespace scope if unsure.

- [ ] **Step 2: Hoist parse-only nodes to namespace scope if needed** — to avoid nested-class friction, move `struct ParseSumNode` and `struct NamedCallNode` to just above `class Parser` (same `detail` namespace) and delete the in-class copies. Update `parseSum`/`parseCall` to construct the namespace-scope types.

- [ ] **Step 3: Build only** — `cmake --build build_catch --target CatchTestRunner -j8`. Expected: compiles.

- [ ] **Step 4: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/expression/Expression.cpp
git commit --no-verify -m "feat(expression): Pratt parser for infix grammar + sum/let (CoolProp-onei)"
```

---

## Task 5: Binder + Program::evaluate (the real `compile`)

The binder walks the parsed tree, resolves names against the 4 buckets, assigns scalar/array slots, validates, fixes sum lengths, and rewrites parse-only nodes (`NameNode`, `ArrayRefNode`, `ParseSumNode`, `NamedCallNode`) into runtime nodes (`ScalarNode`, `ArrayNode`/`IndexNode`, `SumNode`, `CallNode`).

**Files:**
- Modify: `src/expression/Expression.cpp` (replace the stub `ProgramData`, `evaluate`, `compile`)
- Test: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Write failing tests** — append in the expression test region:

```cpp
TEST_CASE("compile evaluates arithmetic, precedence, right-assoc power", "[expression]") {
    using namespace CoolProp::expression;
    CHECK(compile("2 + 3*4", {}, {}).evaluate(nullptr, nullptr) == Approx(14.0));
    CHECK(compile("2^3^2", {}, {}).evaluate(nullptr, nullptr) == Approx(512.0));   // right-assoc: 2^(3^2)
    CHECK(compile("-2^2", {}, {}).evaluate(nullptr, nullptr) == Approx(-4.0));     // -(2^2)
}

TEST_CASE("compile resolves constants and let bindings", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("let y = a*2\ny + 1", {{"a", 5.0}}, {});
    CHECK(p.evaluate(nullptr, nullptr) == Approx(11.0));
}

TEST_CASE("compile sum over co-indexed arrays", "[expression]") {
    using namespace CoolProp::expression;
    // sum(i: a[i]*b[i]) with a={1,2,3}, b={4,5,6} = 4+10+18 = 32
    Program p = compile("sum(i: a[i]*b[i])", {}, {{"a", {1, 2, 3}}, {"b", {4, 5, 6}}});
    CHECK(p.evaluate(nullptr, nullptr) == Approx(32.0));
}

TEST_CASE("compile binds intrinsics in required order", "[expression]") {
    using namespace CoolProp::expression;
    Program p = compile("T + rhomolar", {}, {});
    REQUIRE(p.requiredIntrinsics().size() == 2);
    std::vector<double> iv(2);
    for (std::size_t k = 0; k < 2; ++k)
        iv[k] = (p.requiredIntrinsics()[k] == Intrinsic::T) ? 300.0 : 50.0;
    CHECK(p.evaluate(iv.data(), nullptr) == Approx(350.0));
}

TEST_CASE("compile errors: unknown var, sum length mismatch, bad arity", "[expression]") {
    using namespace CoolProp::expression;
    CHECK_THROWS_AS(compile("nope + 1", {}, {}), CoolProp::ValueError);
    CHECK_THROWS_AS(compile("sum(i: a[i]*b[i])", {}, {{"a", {1, 2}}, {"b", {1, 2, 3}}}),
                    CoolProp::ValueError);
    CHECK_THROWS_AS(compile("exp(1, 2)", {}, {}), CoolProp::ValueError);       // exp is unary
    CHECK_THROWS_AS(compile("a[i]", {}, {{"a", {1.0}}}), CoolProp::ValueError); // index outside sum
}
```

- [ ] **Step 2: Run to verify failure** — `./build_catch/CatchTestRunner "[expression]"` → FAIL (stub returns 7.0).

- [ ] **Step 3: Replace the stub region** in `src/expression/Expression.cpp` — define the real `ProgramData`, the binder, `compile`, `Program::evaluate`. Replace the previous stub `struct ProgramData` and the three `Program`/`compile` definitions with:

```cpp
// --- intrinsic / derived name tables -------------------------------------
static bool intrinsicForName(const std::string& nm, Intrinsic& out) {
    if (nm == "T") { out = Intrinsic::T; return true; }
    if (nm == "rhomolar") { out = Intrinsic::rhomolar; return true; }
    if (nm == "rhomass") { out = Intrinsic::rhomass; return true; }
    if (nm == "molar_mass") { out = Intrinsic::molar_mass; return true; }
    return false;
}
static bool derivedForName(const std::string& nm, Derived& out) {
    if (nm == "p") { out = Derived::p; return true; }   // v1 registry: one entry
    return false;
}
static bool funcForName(const std::string& nm, Func& out, int& arity) {
    struct E { const char* n; Func f; int a; };
    static const E table[] = {
        {"exp", Func::exp, 1},   {"ln", Func::ln, 1},     {"log10", Func::log10, 1},
        {"sqrt", Func::sqrt, 1}, {"abs", Func::abs, 1},   {"pow", Func::pow, 2},
        {"sinh", Func::sinh, 1}, {"cosh", Func::cosh, 1}, {"tanh", Func::tanh, 1},
        {"sin", Func::sin, 1},   {"cos", Func::cos, 1},   {"atan", Func::atan, 1}};
    for (const auto& e : table)
        if (nm == e.n) { out = e.f; arity = e.a; return true; }
    return false;
}

struct ProgramData
{
    int numScalars = 0;
    std::vector<std::pair<int, double>> constantInits;     // (slot, value)
    std::vector<std::pair<Intrinsic, int>> intrinsicSlots; // (kind, slot)
    std::vector<std::pair<Derived, int>> derivedSlots;     // (kind, slot)
    std::vector<std::pair<int, NodePtr>> lets;             // (slot, node) in order
    NodePtr result;
    std::vector<std::vector<double>> arrays;               // by array slot
    // Cached order for the public required* accessors:
    std::vector<Intrinsic> intrinsicOrder;
    std::vector<Derived> derivedOrder;
};

// Binder: produces runtime nodes and fills ProgramData scalar/array tables.
class Binder
{
   public:
    Binder(ProgramData& pd, const std::map<std::string, double>& consts,
           const std::map<std::string, std::vector<double>>& arrays)
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
    ProgramData& d;
    const std::map<std::string, double>& constants;
    const std::map<std::string, std::vector<double>>& arraysIn;
    int nextScalar = 0;
    std::map<std::string, int> letNames;     // let name -> slot
    std::map<std::string, int> constSlots;   // constant name -> slot (deduped)
    std::map<std::string, int> intrinSlots;  // intrinsic name -> slot (deduped)
    std::map<std::string, int> derivSlots;   // derived name -> slot (deduped)
    std::map<std::string, int> arraySlots;   // array name -> slot (deduped)

    int newScalarSlot() { return nextScalar++; }

    int resolveScalar(const std::string& nm, std::size_t /*col*/) {
        auto itL = letNames.find(nm);
        if (itL != letNames.end()) return itL->second;
        auto itC = constants.find(nm);
        if (itC != constants.end()) {
            auto s = constSlots.find(nm);
            if (s != constSlots.end()) return s->second;
            int slot = newScalarSlot();
            constSlots[nm] = slot;
            d.constantInits.emplace_back(slot, itC->second);
            return slot;
        }
        Intrinsic ik;
        if (intrinsicForName(nm, ik)) {
            auto s = intrinSlots.find(nm);
            if (s != intrinSlots.end()) return s->second;
            int slot = newScalarSlot();
            intrinSlots[nm] = slot;
            d.intrinsicSlots.emplace_back(ik, slot);
            d.intrinsicOrder.push_back(ik);
            return slot;
        }
        Derived dk;
        if (derivedForName(nm, dk)) {
            auto s = derivSlots.find(nm);
            if (s != derivSlots.end()) return s->second;
            int slot = newScalarSlot();
            derivSlots[nm] = slot;
            d.derivedSlots.emplace_back(dk, slot);
            d.derivedOrder.push_back(dk);
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

    // bind rewrites a parse tree into a runtime tree. `indexName` is the active
    // sum index (empty when not inside a sum). Returns a node whose eval() is valid.
    NodePtr bind(NodePtr& node, const std::string& indexName) {
        if (auto* num = dynamic_cast<NumNode*>(node.get())) {
            return std::make_unique<NumNode>(num->v);
        }
        if (auto* nm = dynamic_cast<NameNode*>(node.get())) {
            if (!indexName.empty() && nm->name == indexName) return std::make_unique<IndexNode>();
            return std::make_unique<ScalarNode>(resolveScalar(nm->name, 0));
        }
        if (auto* ar = dynamic_cast<ArrayRefNode*>(node.get())) {
            if (indexName.empty() || ar->indexName != indexName)
                throw ValueError(format("array '%s' subscript must be the enclosing sum index",
                                        ar->arrayName.c_str()));
            return std::make_unique<ArrayNode>(resolveArray(ar->arrayName));
        }
        if (auto* neg = dynamic_cast<NegNode*>(node.get())) {
            return std::make_unique<NegNode>(bind(neg->child, indexName));
        }
        if (auto* bn = dynamic_cast<BinNode*>(node.get())) {
            return std::make_unique<BinNode>(bn->op, bind(bn->l, indexName), bind(bn->r, indexName));
        }
        if (auto* call = dynamic_cast<Parser::NamedCallNode*>(node.get())) {
            Func f; int arity;
            if (!funcForName(call->name, f, arity))
                throw ValueError(format("unknown function '%s'", call->name.c_str()));
            if (static_cast<int>(call->args.size()) != arity)
                throw ValueError(format("function '%s' expects %d argument(s)", call->name.c_str(), arity));
            std::vector<NodePtr> bound;
            for (auto& a : call->args) bound.push_back(bind(a, indexName));
            return std::make_unique<CallNode>(f, std::move(bound));
        }
        if (auto* sum = dynamic_cast<Parser::ParseSumNode*>(node.get())) {
            if (!indexName.empty())
                throw ValueError("nested summation is not supported in v1");
            long n = sumLength(*sum, sum->indexName);
            NodePtr body = bind(sum->body, sum->indexName);
            return std::make_unique<SumNode>(n, std::move(body));
        }
        throw ValueError("internal: unhandled parse node in binder");
    }

    // Determine and validate the common length of all arrays subscripted by idx.
    long sumLength(Parser::ParseSumNode& sum, const std::string& idx) {
        long len = -1;
        std::function<void(Node*)> scan = [&](Node* nd) {
            if (auto* ar = dynamic_cast<ArrayRefNode*>(nd)) {
                if (ar->indexName == idx) {
                    auto it = arraysIn.find(ar->arrayName);
                    if (it == arraysIn.end())
                        throw ValueError(format("unknown array '%s'", ar->arrayName.c_str()));
                    long n = static_cast<long>(it->second.size());
                    if (len == -1) len = n;
                    else if (len != n)
                        throw ValueError(format("array '%s' length %ld != %ld in sum",
                                                ar->arrayName.c_str(), n, len));
                }
            } else if (auto* neg = dynamic_cast<NegNode*>(nd)) {
                scan(neg->child.get());
            } else if (auto* bn = dynamic_cast<BinNode*>(nd)) {
                scan(bn->l.get()); scan(bn->r.get());
            } else if (auto* call = dynamic_cast<Parser::NamedCallNode*>(nd)) {
                for (auto& a : call->args) scan(a.get());
            }
        };
        scan(sum.body.get());
        if (len <= 0) throw ValueError("sum body references no array subscripted by the index");
        return len;
    }
};

double Program::evaluate(const double* intrinsicVals, const double* derivedVals) const {
    const detail::ProgramData& d = *m_data;
    std::vector<double> scalars(static_cast<std::size_t>(d.numScalars), 0.0);
    for (const auto& c : d.constantInits) scalars[c.first] = c.second;
    for (std::size_t k = 0; k < d.intrinsicSlots.size(); ++k)
        scalars[d.intrinsicSlots[k].second] = intrinsicVals[k];
    for (std::size_t k = 0; k < d.derivedSlots.size(); ++k)
        scalars[d.derivedSlots[k].second] = derivedVals[k];
    detail::EvalState st{scalars, d.arrays, 0};
    for (const auto& L : d.lets) scalars[L.first] = L.second->eval(st);
    return d.result->eval(st);
}

const std::vector<Intrinsic>& Program::requiredIntrinsics() const { return m_data->intrinsicOrder; }
const std::vector<Derived>& Program::requiredDerived() const { return m_data->derivedOrder; }

Program compile(const std::string& source, const std::map<std::string, double>& constants,
                const std::map<std::string, std::vector<double>>& arrays) {
    using namespace detail;
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
```

> Implementation notes for the binder: `Program::evaluate`/`compile` reference `detail::ProgramData` etc.; ensure the `Program` member functions are defined in the `CoolProp::expression` namespace (not inside `detail`) and qualify `detail::` as shown. The `ProgramData` struct must be visible to `Program`'s out-of-line methods — keep all of this in the single `Expression.cpp` translation unit.

- [ ] **Step 4: Run to verify pass** — `./build_catch/CatchTestRunner "[expression]"` → all expression tests PASS.

- [ ] **Step 5: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/expression/Expression.cpp src/Tests/CoolProp-Tests.cpp
git commit --no-verify -m "feat(expression): binder, 4-bucket resolution, Program::evaluate (CoolProp-onei)"
```

---

## Task 6: ExpressionData struct + enum members in CoolPropFluid.h

**Files:**
- Create: `include/CoolProp/expression/ExpressionCorrelation.h`
- Create: `src/expression/ExpressionCorrelation.cpp`
- Modify: `include/CoolProp/CoolPropFluid.h` (enums at 135-144, 163-170, 233-244, 289-301; structs at 133-148, 161-174, 231-251, 287-306)
- Test: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Create the host wrapper header** `include/CoolProp/expression/ExpressionCorrelation.h`:

```cpp
#ifndef COOLPROP_EXPRESSION_CORRELATION_H
#define COOLPROP_EXPRESSION_CORRELATION_H

#include "CoolProp/expression/Expression.h"

namespace CoolProp {

class HelmholtzEOSMixtureBackend;  // forward decl: header stays light

namespace expression {

/// Host-side wrapper: owns a compiled Program and knows how to fetch the
/// intrinsic state and registered derived quantities from an EOS backend.
class ExpressionCorrelation
{
   public:
    ExpressionCorrelation() = default;
    explicit ExpressionCorrelation(Program prog) : m_program(std::move(prog)), m_set(true) {}
    bool is_set() const { return m_set; }
    /// Evaluate the formula at the backend's current state; returns base-SI result.
    double eval(HelmholtzEOSMixtureBackend& HEOS) const;

   private:
    Program m_program;
    bool m_set = false;
};

}  // namespace expression

/// Stored in each transport sub-correlation container. Empty unless type==expression.
struct ExpressionData
{
    std::shared_ptr<expression::ExpressionCorrelation> correlation;
};

}  // namespace CoolProp

#endif
```

- [ ] **Step 2: Create the host wrapper impl** `src/expression/ExpressionCorrelation.cpp`:

```cpp
#include "CoolProp/expression/ExpressionCorrelation.h"
#include "CoolProp/Backends/Helmholtz/HelmholtzEOSMixtureBackend.h"

namespace CoolProp {
namespace expression {

static double fetchIntrinsic(HelmholtzEOSMixtureBackend& H, Intrinsic k) {
    switch (k) {
        case Intrinsic::T: return H.T();
        case Intrinsic::rhomolar: return H.rhomolar();
        case Intrinsic::rhomass: return H.rhomass();
        case Intrinsic::molar_mass: return H.molar_mass();
    }
    return 0.0;
}
static double fetchDerived(HelmholtzEOSMixtureBackend& H, Derived k) {
    switch (k) {
        case Derived::p: return H.p();
    }
    return 0.0;
}

double ExpressionCorrelation::eval(HelmholtzEOSMixtureBackend& HEOS) const {
    const std::vector<Intrinsic>& ik = m_program.requiredIntrinsics();
    const std::vector<Derived>& dk = m_program.requiredDerived();
    std::vector<double> iv(ik.size()), dv(dk.size());
    for (std::size_t i = 0; i < ik.size(); ++i) iv[i] = fetchIntrinsic(HEOS, ik[i]);
    for (std::size_t i = 0; i < dk.size(); ++i) dv[i] = fetchDerived(HEOS, dk[i]);
    return m_program.evaluate(iv.empty() ? nullptr : iv.data(),
                              dv.empty() ? nullptr : dv.data());
}

}  // namespace expression
}  // namespace CoolProp
```

- [ ] **Step 3: Add the include + `ExpressionData` member + enum values to `include/CoolProp/CoolPropFluid.h`.** Near the top with the other includes, add:

```cpp
#include "CoolProp/expression/ExpressionCorrelation.h"
```

Then, in the four containers, add an `expression_data` member and the matching enum value. Concretely:

- Enum `ViscosityDiluteType` (~233-244): add `VISCOSITY_DILUTE_EXPRESSION,` immediately before `VISCOSITY_DILUTE_NOT_SET`.
- Struct `ViscosityDiluteVariables` (~231-251): add member `ExpressionData expression_data;`
- Enum `ViscosityHigherOrderEnum` (~289-301): add `VISCOSITY_HIGHER_ORDER_EXPRESSION,` before `..._NOT_SET`.
- Struct `ViscosityHigherOrderVariables` (~287-306): add `ExpressionData expression_data;`
- Enum `ConductivityDiluteEnum` (~135-144): add `CONDUCTIVITY_DILUTE_EXPRESSION,` before `..._NOT_SET`.
- Struct `ConductivityDiluteVariables` (~133-148): add `ExpressionData expression_data;`
- Enum `ConductivityResidualEnum` (~163-170): add `CONDUCTIVITY_RESIDUAL_EXPRESSION,` before `..._NOT_SET`.
- Struct `ConductivityResidualVariables` (~161-174): add `ExpressionData expression_data;`

- [ ] **Step 4: Write failing test** — append:

```cpp
TEST_CASE("ExpressionData default-constructs unset", "[expression]") {
    CoolProp::ExpressionData d;
    CHECK(!d.correlation);
}
```

- [ ] **Step 5: Build & run** — `cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[expression]"` → PASS. (This also proves `CoolPropFluid.h` still compiles with the new include + members.)

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add include/CoolProp/expression/ExpressionCorrelation.h src/expression/ExpressionCorrelation.cpp include/CoolProp/CoolPropFluid.h src/Tests/CoolProp-Tests.cpp
git commit --no-verify -m "feat(expression): ExpressionData container + transport enum/struct hooks (CoolProp-onei)"
```

---

## Task 7: JSON parse branches (FluidLibrary)

Add an `else if (!type.compare("expression"))` branch to each of the four parse functions. The branch reads `formula` (string), optional `constants` (object→map<string,double>), optional `arrays` (object→map<string,vector<double>>), calls `expression::compile`, and stores the `ExpressionCorrelation` + sets the type enum.

**Files:**
- Modify: `src/Backends/Helmholtz/Fluids/FluidLibrary.h` (dilute visc ~460-522; higher-order visc ~556-666; dilute cond ~796-843; residual cond ~846-887)
- Test: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Add a small JSON helper** near the top of the anonymous/parse scope in `FluidLibrary.h` (after the existing includes; reuse `cpjson`/nlohmann already in scope). Add a static inline helper:

```cpp
// Build an ExpressionCorrelation from a transport sub-block of type "expression".
static inline CoolProp::ExpressionData parse_expression_block(const nlohmann::json& j,
                                                              const std::string& fluidname) {
    std::map<std::string, double> constants;
    std::map<std::string, std::vector<double>> arrays;
    if (j.contains("constants"))
        for (auto it = j["constants"].begin(); it != j["constants"].end(); ++it)
            constants[it.key()] = it.value().get<double>();
    if (j.contains("arrays"))
        for (auto it = j["arrays"].begin(); it != j["arrays"].end(); ++it)
            arrays[it.key()] = it.value().get<std::vector<double>>();
    std::string formula = cpjson::get_string(j, "formula");
    try {
        CoolProp::expression::Program prog = CoolProp::expression::compile(formula, constants, arrays);
        CoolProp::ExpressionData data;
        data.correlation = std::make_shared<CoolProp::expression::ExpressionCorrelation>(std::move(prog));
        return data;
    } catch (std::exception& e) {
        throw ValueError(format("expression compile failed for fluid %s: %s",
                                fluidname.c_str(), e.what()));
    }
}
```

Ensure `#include "CoolProp/expression/ExpressionCorrelation.h"` is present at the top of `FluidLibrary.h` (it is transitively via `CoolPropFluid.h`, but add it explicitly for clarity).

- [ ] **Step 2: Add the four branches.** In each parse function, immediately before the final `else { throw ValueError(...type not understood...); }`, insert:

In `parse_dilute_viscosity` (before its final else, ~line 519):
```cpp
    } else if (!type.compare("expression")) {
        fluid.transport.viscosity_dilute.expression_data = parse_expression_block(dilute, fluid.name);
        fluid.transport.viscosity_dilute.type = CoolProp::ViscosityDiluteVariables::VISCOSITY_DILUTE_EXPRESSION;
```

In `parse_higher_order_viscosity` (before its final else, ~line 663):
```cpp
    } else if (!type.compare("expression")) {
        fluid.transport.viscosity_higher_order.expression_data = parse_expression_block(higher, fluid.name);
        fluid.transport.viscosity_higher_order.type = CoolProp::ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_EXPRESSION;
```

In `parse_dilute_conductivity` (before its final else, ~line 840):
```cpp
    } else if (!type.compare("expression")) {
        fluid.transport.conductivity_dilute.expression_data = parse_expression_block(dilute, fluid.name);
        fluid.transport.conductivity_dilute.type = CoolProp::ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_EXPRESSION;
```

In `parse_residual_conductivity` (before its final else, ~line 884):
```cpp
    } else if (!type.compare("expression")) {
        fluid.transport.conductivity_residual.expression_data = parse_expression_block(dilute, fluid.name);
        fluid.transport.conductivity_residual.type = CoolProp::ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_EXPRESSION;
```

> Note: `conductivity_dilute.type` and `conductivity_residual.type` are `int` fields (recon §2) — assigning the enum value is valid (implicit enum→int).

- [ ] **Step 3: Failing test** — verify a hand-built fluid JSON string with an `expression` viscosity block loads. Append:

```cpp
TEST_CASE("FluidLibrary parses an expression transport block", "[expression]") {
    // Minimal: confirm compile() is reachable through the parse helper indirectly
    // by constructing the same maps the parser would.
    using namespace CoolProp::expression;
    std::map<std::string, double> consts{{"T_reduce", 132.0}, {"rhomolar_reduce", 10000.0}};
    std::map<std::string, std::vector<double>> arrays{{"a", {1.0e-5}}, {"d1", {1.0}}, {"t1", {0.2}}};
    Program p = compile(
        "let delta = rhomolar/rhomolar_reduce\nlet tau = T_reduce/T\nsum(i: a[i]*delta^d1[i]*tau^t1[i])",
        consts, arrays);
    // T=300, rhomolar=5000 -> delta=0.5, tau=0.44 ; 1e-5*0.5^1*0.44^0.2
    std::vector<double> iv(p.requiredIntrinsics().size());
    for (std::size_t k = 0; k < iv.size(); ++k) {
        switch (p.requiredIntrinsics()[k]) {
            case Intrinsic::T: iv[k] = 300.0; break;
            case Intrinsic::rhomolar: iv[k] = 5000.0; break;
            default: iv[k] = 0.0; break;
        }
    }
    double expected = 1.0e-5 * std::pow(0.5, 1.0) * std::pow(132.0 / 300.0, 0.2);
    CHECK(p.evaluate(iv.data(), nullptr) == Approx(expected));
}
```

- [ ] **Step 4: Build & run** — `cmake --build build_catch --target CatchTestRunner -j8 && ./build_catch/CatchTestRunner "[expression]"` → PASS. (FluidLibrary.h must still compile — that proves the branches are syntactically wired.)

- [ ] **Step 5: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/Helmholtz/Fluids/FluidLibrary.h src/Tests/CoolProp-Tests.cpp
git commit --no-verify -m "feat(expression): JSON parse branches for type:expression (CoolProp-onei)"
```

---

## Task 8: Evaluation dispatch cases (HelmholtzEOSMixtureBackend)

Add a `case *_EXPRESSION` to each of the five dispatch switches that calls the stored `ExpressionCorrelation::eval(*this)`.

**Files:**
- Modify: `src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp` (dilute visc ~713-749; higher-order visc ~773-799; dilute cond ~984-1006; residual cond ~1040-1051)
- Test: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Ensure the header include** at the top of `HelmholtzEOSMixtureBackend.cpp` (add if absent):

```cpp
#include "CoolProp/expression/ExpressionCorrelation.h"
```

- [ ] **Step 2: Add the dispatch cases.** Insert before each `default:` (or before the closing `}` where there's no default):

Dilute viscosity switch (~745):
```cpp
            case ViscosityDiluteVariables::VISCOSITY_DILUTE_EXPRESSION:
                eta_dilute = components[0].transport.viscosity_dilute.expression_data.correlation->eval(*this);
                break;
```

Higher-order viscosity switch (~797):
```cpp
            case ViscosityHigherOrderVariables::VISCOSITY_HIGHER_ORDER_EXPRESSION:
                residual = components[0].transport.viscosity_higher_order.expression_data.correlation->eval(*this);
                break;
```

Dilute conductivity switch (~1004):
```cpp
            case ConductivityDiluteVariables::CONDUCTIVITY_DILUTE_EXPRESSION:
                dilute = components[0].transport.conductivity_dilute.expression_data.correlation->eval(*this);
                break;
```

Residual conductivity switch (~1049):
```cpp
            case ConductivityResidualVariables::CONDUCTIVITY_RESIDUAL_EXPRESSION:
                lambda_residual = components[0].transport.conductivity_residual.expression_data.correlation->eval(*this);
                break;
```

- [ ] **Step 3: Build** — `cmake --build build_catch --target CatchTestRunner -j8`. Expected: compiles (proves the dispatch + header wiring). No new unit test here; Task 9 provides the end-to-end golden check.

- [ ] **Step 4: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp
git commit --no-verify -m "feat(expression): dispatch cases for *_EXPRESSION transport types (CoolProp-onei)"
```

---

## Task 9: Golden regression — DSL reproduces every Tier-A form

For each Tier-A form, instantiate a backend for a fluid that uses it, read the loaded coefficients straight from the transport struct, build the DSL equivalent with the SAME coefficients, and compare the DSL result against the shipped `TransportRoutines::*` routine (public static) across a `(T, ρ)` grid.

**Files:**
- Test: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Pick a representative fluid per form** — run, at task time:

```bash
for t in powers_of_T powers_of_Tr collision_integral collision_integral_powers_of_Tstar \
         modified_Batschinski_Hildebrand ratio_of_polynomials eta0_and_poly \
         polynomial polynomial_and_exponential; do
  echo "== $t =="; grep -rl "\"$t\"" dev/fluids | head -3
done
```
Record one fluid name per form (e.g. for `modified_Batschinski_Hildebrand`, the recon example used parameters matching `Nitrogen`/`Air`-class fluids — confirm with the grep).

- [ ] **Step 2: Write the golden test for `powers_of_T`** (template the rest after it). Append:

```cpp
#    include "CoolProp/Backends/Helmholtz/TransportRoutines.h"

static std::shared_ptr<CoolProp::HelmholtzEOSMixtureBackend> make_HEOS(const std::string& fluid) {
    std::vector<std::string> names(1, fluid);
    auto HEOS = std::make_shared<CoolProp::HelmholtzEOSMixtureBackend>(names);
    return HEOS;
}

TEST_CASE("DSL reproduces viscosity_dilute_powers_of_T", "[expression]") {
    using namespace CoolProp;
    // Replace "FLUID_FROM_STEP_1" with the grep-confirmed fluid that uses powers_of_T.
    auto HEOS = make_HEOS("FLUID_FROM_STEP_1");
    auto& d = HEOS->components[0].transport.viscosity_dilute.powers_of_T;
    std::map<std::string, std::vector<double>> arrays{
        {"a", std::vector<double>(d.a.begin(), d.a.end())},
        {"t", std::vector<double>(d.t.begin(), d.t.end())}};
    auto prog = expression::compile("sum(i: a[i]*T^t[i])", {}, arrays);

    for (double T : {200.0, 300.0, 400.0, 500.0}) {
        for (double rho : {1.0, 100.0, 5000.0}) {
            HEOS->update(DmolarT_INPUTS, rho, T);
            double expected = TransportRoutines::viscosity_dilute_powers_of_T(*HEOS);
            std::vector<double> iv(prog.requiredIntrinsics().size());
            for (std::size_t k = 0; k < iv.size(); ++k)
                iv[k] = (prog.requiredIntrinsics()[k] == expression::Intrinsic::T) ? T : rho;
            double got = prog.evaluate(iv.empty() ? nullptr : iv.data(), nullptr);
            INFO("T=" << T << " rho=" << rho);
            CHECK(got == Approx(expected).epsilon(1e-14));
        }
    }
}
```

- [ ] **Step 3: Add golden tests for the remaining Tier-A forms**, each following the Step-2 template but with the form's DSL string and the coefficients read from its struct. Use these DSL strings (coefficients/constants pulled from the corresponding struct fields — confirm field names by reading `CoolPropFluid.h` for each `*Data` struct):

  - `powers_of_Tr`: `let Tr = T/T_reducing\nsum(i: a[i]*Tr^t[i])` — `T_reducing` from `data.T_reducing`.
  - `collision_integral`: `let Tstar = T/epsilon_over_k\nlet S = exp(sum(i: a[i]*ln(Tstar)^t[i]))\nC*sqrt(molar_mass*1000*T)/((sigma_eta*1e9)^2*S)` — constants `C`, `molar_mass` from struct; `epsilon_over_k`, `sigma_eta` from `components[0].transport`.
  - `collision_integral_powers_of_Tstar`: `let Tstar = T/T_reducing\nC*sqrt(T)/sum(i: a[i]*Tstar^t[i])`.
  - `modified_Batschinski_Hildebrand`: the full two-sum form — read `a,d1,t1,gamma,l,f,d2,t2,g1,g2,...` from the struct (read `TransportRoutines.cpp:103-138` for the exact algebra) and reproduce it verbatim in the DSL, with `delta`/`tau`/`delta0` as `let` bindings.
  - `ratio_of_polynomials`: `let Tr = T/T_reducing\nsum(i: A[i]*Tr^n[i])/sum(i: B[i]*Tr^m[i])`.
  - `eta0_and_poly`: read `conductivity_dilute_eta0_and_poly` algebra from `TransportRoutines.cpp` and reproduce.
  - `polynomial` / `polynomial_and_exponential`: read the residual-conductivity algebra from `TransportRoutines.cpp` and reproduce, using `rhomass`/`T_reducing`/`rhomass_reducing` as needed.

  For forms whose C++ uses a hand-written `x*x` (integer power) rather than `pow`, assert bit-exact:
  `CHECK(got == expected);` ; otherwise keep `.epsilon(1e-14)`.

- [ ] **Step 4: Run** — `./build_catch/CatchTestRunner "[expression]"` → all PASS. If any form diverges beyond `1e-14`, the DSL string does not match the C++ algebra — fix the DSL string (this is the completeness proof working as intended).

- [ ] **Step 5: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Tests/CoolProp-Tests.cpp
git commit --no-verify -m "test(expression): golden regression vs hardcoded Tier-A routines (CoolProp-onei)"
```

---

## Task 10: Registry `p` path + error-mode coverage

**Files:**
- Test: `src/Tests/CoolProp-Tests.cpp`

- [ ] **Step 1: Write the `p` registry test** — proves bucket-4 resolution and that the host injects pressure from HEOS. Append:

```cpp
TEST_CASE("DSL derived variable p resolves and matches HEOS.p()", "[expression]") {
    using namespace CoolProp;
    auto HEOS = make_HEOS("Water");
    HEOS->update(DmolarT_INPUTS, 30000.0, 350.0);
    auto corr = std::make_shared<expression::ExpressionCorrelation>(
        expression::compile("p", {}, {}));
    double got = corr->eval(*HEOS);
    CHECK(got == Approx(HEOS->p()));
    // p must be reported as a required-derived dependency:
    auto prog = expression::compile("p*2", {}, {});
    REQUIRE(prog.requiredDerived().size() == 1);
    CHECK(prog.requiredDerived()[0] == expression::Derived::p);
}
```

- [ ] **Step 2: Confirm error-mode tests** already added in Task 5 cover unknown-var, sum-length-mismatch, bad-arity, and index-outside-sum. Add two more for robustness — append:

```cpp
TEST_CASE("DSL compile errors are clean, never crashes", "[expression]") {
    using namespace CoolProp::expression;
    CHECK_THROWS_AS(compile("2 +", {}, {}), CoolProp::ValueError);         // dangling operator
    CHECK_THROWS_AS(compile("foo(1)", {}, {}), CoolProp::ValueError);      // unknown function
    CHECK_THROWS_AS(compile("(1 + 2", {}, {}), CoolProp::ValueError);      // unbalanced paren
    CHECK_THROWS_AS(compile("sum(i: 1)", {}, {}), CoolProp::ValueError);   // sum with no array
}
```

- [ ] **Step 3: Run** — `./build_catch/CatchTestRunner "[expression]"` → all PASS.

- [ ] **Step 4: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Tests/CoolProp-Tests.cpp
git commit --no-verify -m "test(expression): derived-p registry path + error-mode coverage (CoolProp-onei)"
```

---

## Task 11: End-to-end fluid round-trip + preflight

Prove an actual fluid JSON with `"type": "expression"` loads and produces the right viscosity end-to-end through the public API, then run the full gate.

**Files:**
- Test: `src/Tests/CoolProp-Tests.cpp`
- (No production JSON is shipped; the test builds a fluid in-memory or uses a temp JSON via the existing fluid-injection path.)

- [ ] **Step 1: Determine the in-memory fluid-add path** — read how `CoolProp::add_fluids_as_JSON` or equivalent is used in existing tests (grep `add_fluids_as_JSON` in `src/Tests/CoolProp-Tests.cpp` and `src/CoolProp.cpp`). Use that to register a minimal modified copy of an existing fluid where one transport block is replaced with the `"type":"expression"` equivalent that reproduces the original coefficients.

- [ ] **Step 2: Write the round-trip test** — load the modified fluid, compute viscosity at a state via `PropsSI("V", "T", T, "Dmolar", rho, fluidname)`, and compare to the unmodified fluid's viscosity at the same state. Assert relative error `< 1e-14`. (Exact code depends on the Step-1 path; write it concretely once that path is known — full assertions, no placeholders.)

- [ ] **Step 3: Run the umbrella transport tests** to prove no regression in existing correlations:

```bash
./build_catch/CatchTestRunner "[expression]"
./build_catch/CatchTestRunner "[transport]"   # if such a tag exists; else run viscosity/conductivity validation cases
```
Expected: all PASS.

- [ ] **Step 4: Run preflight (REQUIRED before any push)** —

```bash
./dev/ci/preflight.sh
```
Expected: clang-format clean, build OK, selected tag scope passes, cppcheck/clang-tidy/semgrep clean on changed files. Fix anything flagged.

- [ ] **Step 5: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Tests/CoolProp-Tests.cpp
git commit --no-verify -m "test(expression): end-to-end fluid round-trip via PropsSI (CoolProp-onei)"
```

---

## Task 12: Performance benchmark — DSL vs hardcoded routine

Prove the DSL evaluates with throughput comparable to the hardcoded C++ routine. The heaviest Tier-A form (modified Batschinski–Hildebrand, two sums + `exp`) is the worst case, so benchmark that. The main overhead risk is the per-call `std::vector<double> scalars` allocation in `Program::evaluate`; the benchmark decides whether the scratch-buffer optimization is warranted.

**Files:**
- Test: `src/Tests/CoolProp-Tests.cpp`
- Possibly modify: `src/expression/Expression.cpp` (scratch-buffer optimization, only if needed)

- [ ] **Step 1: Write a self-contained timing benchmark** (no Catch2 benchmark macro — a portable `std::chrono` loop with an informational report and a soft regression gate). Append to the expression test region:

```cpp
#    include <chrono>

TEST_CASE("DSL throughput is comparable to the hardcoded routine", "[expression][!benchmark]") {
    using namespace CoolProp;
    // Use the same fluid/form as the modified_Batschinski_Hildebrand golden test.
    auto HEOS = make_HEOS("FLUID_FROM_TASK9_BH");
    HEOS->update(DmolarT_INPUTS, 5000.0, 300.0);

    // Build the DSL equivalent (same string as the BH golden test, same coefficients).
    auto corr = std::make_shared<expression::ExpressionCorrelation>(/* compile(... BH form ...) */);

    const int N = 2'000'000;
    // Hardcoded baseline
    volatile double sink = 0.0;
    auto t0 = std::chrono::steady_clock::now();
    for (int k = 0; k < N; ++k)
        sink += TransportRoutines::viscosity_higher_order_modified_Batschinski_Hildebrand(*HEOS);
    auto t1 = std::chrono::steady_clock::now();
    // DSL
    for (int k = 0; k < N; ++k) sink += corr->eval(*HEOS);
    auto t2 = std::chrono::steady_clock::now();

    double ns_hard = std::chrono::duration<double, std::nano>(t1 - t0).count() / N;
    double ns_dsl = std::chrono::duration<double, std::nano>(t2 - t1).count() / N;
    double ratio = ns_dsl / ns_hard;
    WARN("hardcoded = " << ns_hard << " ns/eval; DSL = " << ns_dsl
                        << " ns/eval; ratio = " << ratio << "x");
    // Soft gate: catch a gross regression (e.g. accidental O(n^2) or pathological alloc).
    CHECK(ratio < 5.0);
}
```

- [ ] **Step 2: Run the benchmark** (benchmarks are tagged `[!benchmark]`, so run them explicitly):

```bash
./build_catch/CatchTestRunner "[expression][!benchmark]"
```
Record the reported `ns/eval` for both and the ratio.

- [ ] **Step 3: If `ratio >= ~3x`, apply the scratch-buffer optimization** — eliminate the per-call heap allocation in `Program::evaluate`. In `src/expression/Expression.cpp`, give `ProgramData` a reusable, thread-safe scratch buffer instead of a fresh `std::vector` each call:

```cpp
// In Program::evaluate, replace `std::vector<double> scalars(...)` with a
// thread_local buffer sized to numScalars (reentrant across states, no per-call alloc):
static thread_local std::vector<double> scalars;
scalars.assign(static_cast<std::size_t>(d.numScalars), 0.0);
```

Re-run Step 2 and confirm the ratio improved and is `< 3x`. (If already `< 3x` in Step 2, skip this step and note the measured ratio in the commit message.)

- [ ] **Step 4: Re-run the correctness gate** to confirm the optimization (if applied) changed nothing functionally:

```bash
./build_catch/CatchTestRunner "[expression]"
```
Expected: all PASS.

- [ ] **Step 5: Re-run preflight** (required because Step 3 may have changed production code):

```bash
./dev/ci/preflight.sh
```

- [ ] **Step 6: Commit**

```bash
git restore --staged .beads/issues.jsonl 2>/dev/null; git checkout .beads/issues.jsonl 2>/dev/null
git add src/Tests/CoolProp-Tests.cpp src/expression/Expression.cpp
git commit --no-verify -m "perf(expression): benchmark DSL vs hardcoded BH; ratio <measured>x (CoolProp-onei)"
```

---

## Self-Review (completed during planning)

**Spec coverage:**
- Tier-A completeness → Task 9 (golden regression per form).
- `sum(i: …)` explicit-index grammar → Tasks 4 (parse) + 5 (bind, length validation).
- `let` bindings, raw-state-only intrinsics → Task 5 binder buckets.
- 4-bucket resolution incl. derived registry seeded with `p` → Task 5 (`derivedForName`) + Task 10 (`p` path).
- Base-SI everywhere, no units field → JSON schema in Task 7 has no units key.
- Tree-walking AST, compile-once → Tasks 3/5; `evaluate` builds a fresh scalar scratch per call (reentrant; per-call alloc noted as a future micro-opt, consistent with spec risk section).
- EOS-free evaluator; host fills context → Task 6 (`ExpressionCorrelation` is the only HEOS-aware unit).
- Additive integration, hardcoded untouched → Tasks 6/7/8 add members/branches/cases only.
- `1e-14` gate + per-form bit-exact → Task 9 Step 3.
- Clean compile errors, never crash → Tasks 2/5/10 error tests.
- Throughput comparable to hardcoded routines → Task 12 (benchmark vs the heaviest BH form; soft gate + scratch-buffer optimization fallback to remove per-call allocation).

**Placeholder scan:** Two deliberate, bounded discovery steps remain — Task 9 Step 1 (grep to pick the representative fluid per form) and Task 11 Steps 1-2 (locate the in-memory fluid-add path, then write the round-trip assertions). Both are concrete instructions with the exact command/grep to run, not vague "handle X" placeholders; the executor resolves them deterministically at task time. Every code-bearing step includes complete code.

**Type consistency:** `Program`, `compile`, `evaluate(const double*, const double*)`, `requiredIntrinsics/Derived`, `Intrinsic`/`Derived` enums, `ExpressionCorrelation::eval(HelmholtzEOSMixtureBackend&)`, `ExpressionData::correlation` (shared_ptr) are used identically across Tasks 1-11. Parse-only nodes (`NameNode`, `ArrayRefNode`, `ParseSumNode`, `NamedCallNode`) are consumed only by the binder and never reach `evaluate`.

**One known build-ordering caveat for the executor:** `Program`'s out-of-line methods and `compile` reference `detail::ProgramData`/`Binder`/`Parser`; keep them all in the single `Expression.cpp` TU and define `Program`'s methods after `ProgramData` is complete. If `Parser::ParseSumNode`/`NamedCallNode` cause nested-class access friction in the binder, hoist them to namespace scope (Task 4 Step 2 already directs this).
