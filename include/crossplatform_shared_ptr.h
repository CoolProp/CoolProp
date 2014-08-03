#ifndef CROSSPLATFORM_SHARED_PTR
#define CROSSPLATFORM_SHARED_PTR

#include "PlatformDetermination.h"

// Based on the platform and compiler, include the necessary header to give access to std::tr1::shared_ptr directly as shared_ptr

#if defined(__ISLINUX__) && (defined(__llvm__) || defined(__clang__))
    #include <memory>
    using std::shared_ptr;
#elif defined(__ISLINUX__)
    #include <tr1/memory>
    using namespace std::tr1;
#elif defined(__ISAPPLE__)
    // See also http://stackoverflow.com/questions/13445742/apple-and-shared-ptr
    #include <tr1/memory>
    using std::tr1::shared_ptr;
#elif defined(__GNUC__)
    #include <tr1/memory>
    using namespace std::tr1;
#elif defined(__ISWINDOWS__) && defined(__MINGW32__)
    #include <tr1/memory>
    using namespace std::tr1;
#elif defined(__ISWINDOWS__) && !defined(__MINGW32__)
    #include <memory>
    using namespace std::tr1;
#else
    #pragma error
#endif



#endif
