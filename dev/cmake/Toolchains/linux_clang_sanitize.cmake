# run with cmake -DCMAKE_TOOLCHAIN_FILE=...

find_program(CMAKE_C_COMPILER clang)
find_program(CMAKE_CXX_COMPILER clang++)

if(NOT CMAKE_C_COMPILER)
  message(FATAL_ERROR "clang not found")
endif()

if(NOT CMAKE_CXX_COMPILER)
  message(FATAL_ERROR "clang++ not found")
endif()

set(
    CMAKE_C_COMPILER
    "${CMAKE_C_COMPILER}"
    CACHE
    STRING
    "C compiler"
    FORCE
)

set(
    CMAKE_CXX_COMPILER
    "${CMAKE_CXX_COMPILER}"
    CACHE
    STRING
    "C++ compiler"
    FORCE
)

polly_add_cache_flag(CMAKE_CXX_FLAGS "-fsanitize=memory")
polly_add_cache_flag(CMAKE_CXX_FLAGS "-fsanitize-memory-track-origins")
polly_add_cache_flag(CMAKE_CXX_FLAGS "-g")


include(polly_add_cache_flag)
include(polly_fatal_error)

polly_add_cache_flag(CMAKE_CXX_FLAGS "-fsanitize=address")
polly_add_cache_flag(CMAKE_CXX_FLAGS "-g")

set(
    CMAKE_CXX_FLAGS_RELEASE
    "-O1 -DNDEBUG"
    CACHE
    STRING
    "C++ compiler flags"
    FORCE
)

polly_add_cache_flag(CMAKE_C_FLAGS "-fsanitize=address")
polly_add_cache_flag(CMAKE_C_FLAGS "-g")

set(
    CMAKE_C_FLAGS_RELEASE
    "-O1 -DNDEBUG"
    CACHE
    STRING
    "C compiler flags"
    FORCE
)


polly_add_cache_flag(CMAKE_CXX_FLAGS "-fsanitize=leak")
polly_add_cache_flag(CMAKE_CXX_FLAGS "-g")

polly_add_cache_flag(CMAKE_C_FLAGS "-fsanitize=leak")
polly_add_cache_flag(CMAKE_C_FLAGS "-g")

list(APPEND HUNTER_TOOLCHAIN_UNDETECTABLE_ID "sanitize-leak")



polly_add_cache_flag(CMAKE_CXX_FLAGS "-fsanitize=thread")
polly_add_cache_flag(CMAKE_CXX_FLAGS "-g")

polly_add_cache_flag(CMAKE_C_FLAGS "-fsanitize=thread")
polly_add_cache_flag(CMAKE_C_FLAGS "-g")