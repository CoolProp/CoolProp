#ifndef CROSSPLATFORM_SHARED_PTR
#define CROSSPLATFORM_SHARED_PTR

// By default, we use shared_ptr from the std namespace, and include the memory header,
// but some compilers need different treatment. Cmake provides the tools to
// ensure that the correct header is identified as a compile-time check, and we use
// that capability to change the include and/or the namespace

#if defined(SHARED_PTR_TR1_MEMORY_HEADER)
#    include <tr1/memory>
#else
#    include <memory>
#endif

#if defined(SHARED_PTR_TR1_NAMESPACE)
using std::tr1::shared_ptr;
#else
using std::shared_ptr;
#endif

#endif