#ifndef CROSSPLATFORM_SHARED_PTR
#define CROSSPLATFORM_SHARED_PTR

#include "PlatformDetermination.h"

// Based on the platform and compiler, include the necessary header to give access to std::tr1::shared_ptr directly as shared_ptr

#if defined(__ANDROID__)
        #include <memory>
        using std::shared_ptr;
#elif defined(__ISLINUX__) && (defined(__llvm__) || defined(__clang__)) // CLANG
    #if __has_include(<tr1/memory>)
        // CLANG and -stdlib=libstdc++
        // See also http://stackoverflow.com/questions/13445742/apple-and-shared-ptr
        #include <tr1/memory>
        using std::tr1::shared_ptr;
    #else
        // CLANG and -stdlib=libc++
        #include <memory>
        using std::shared_ptr;
    #endif
#elif defined(__ISLINUX__) // GCC
    #include <memory>
    using std::shared_ptr;
#elif defined(__ISAPPLE__) && (defined(__llvm__) || defined(__clang__)) // CLANG
    // See docs for clang: http://clang.llvm.org/docs/LanguageExtensions.html#include-file-checking-macros
    #if __has_include(<tr1/memory>)
        // CLANG and -stdlib=libstdc++
        // See also http://stackoverflow.com/questions/13445742/apple-and-shared-ptr
        #include <tr1/memory>
        using std::tr1::shared_ptr;
    #else
        // CLANG and -stdlib=libc++
        #include <memory>
        using std::shared_ptr;
    #endif
#elif defined(__GNUC__)
    #include <tr1/memory>
    using std::tr1::shared_ptr;
#elif defined(__ISWINDOWS__) && defined(__MINGW32__)
    #include <tr1/memory>
    using std::tr1;
#elif defined(__ISWINDOWS__) && !defined(__MINGW32__)
    #include <memory>
    // VS2008 has std::shared_ptr from C++11
    #if defined(_MSC_VER) && _MSC_VER >= 1600
        using std::shared_ptr;
    #else
        using std::tr1::shared_ptr;
    #endif
#else
    #pragma error
#endif



#endif
