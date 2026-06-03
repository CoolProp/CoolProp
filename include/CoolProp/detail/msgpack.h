// Workaround MSVC warnings
#ifdef _MSC_VER
#    pragma warning(push)
#    pragma warning(disable : 4267)
#endif

//// Workaround MSVC endiannes issues
//#if defined(_MSC_VER) && ( defined(_M_ARM) || defined(_M_ARM64) )
//#    define MSGPACK_ENDIAN_LITTLE_BYTE
//#endif

#include "msgpack.hpp"

/* #if defined (MSGPACK_ENDIAN_LITTLE_BYTE)
#    undef MSGPACK_ENDIAN_LITTLE_BYTE
#endif*/

#ifdef _MSC_VER
#    pragma warning(pop)
#endif
