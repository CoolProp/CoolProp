# run with cmake -DCMAKE_TOOLCHAIN_FILE=...


polly_add_cache_flag(CMAKE_CXX_FLAGS "-fsanitize=memory")
polly_add_cache_flag(CMAKE_CXX_FLAGS "-fsanitize-memory-track-origins")
polly_add_cache_flag(CMAKE_CXX_FLAGS "-g")


