# CoolProp CI scripts and contributor tooling

This directory holds CI helper scripts and contributor-facing tooling docs.
The full code-quality workflow doc lives elsewhere (CoolProp-2uw.4); this file
covers the pieces that CI directly depends on.

## compile_commands.json

`clang-tidy`, `include-what-you-use`, and most editor LSP integrations need a
JSON compilation database. CoolProp configures `CMAKE_EXPORT_COMPILE_COMMANDS`
ON unconditionally (CMakeLists.txt), so any cmake configure produces it:

```bash
cmake -G Ninja -B build -S .            # Ninja or Makefile generator
ls build/compile_commands.json          # 90+ entries covering src/
```

### macOS contributors

The Homebrew `llvm@18` and `llvm` packages ship `clang-tidy`, but they don't
know about the Xcode SDK that Apple's `/usr/bin/c++` links against. If
`clang-tidy` reports `'iterator' file not found` or `__builtin_clzg`-style
errors against system headers, point it at the Xcode sysroot:

```bash
SDK=$(xcrun --show-sdk-path)
clang-tidy -p build --extra-arg=--sysroot=$SDK src/CPstrings.cpp
```

clang-tidy 19+ is recommended on macOS — earlier versions don't recognize new
libc++ builtins (`__builtin_clzg`, `__builtin_ctzg`) that Apple's libc++
headers use.

CI runs on Ubuntu where this issue doesn't apply.

## clang-format

`clang-format.sh` runs `clang-format -i` over C/C++ files changed between two
git refs. Used by `.github/workflows/dev_clangformat.yml` and runnable
locally:

```bash
./dev/ci/clang-format.sh HEAD master   # format files staged + committed vs master
```

Pinned version: clang-format 18.1.x (matches the Ubuntu image in CI).
