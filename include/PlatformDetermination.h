#ifndef PLATFORMDETERMINATION_H
#define PLATFORMDETERMINATION_H

// See also http://stackoverflow.com/questions/5919996/how-to-detect-reliably-mac-os-x-ios-linux-windows-in-c-preprocessor
#if _WIN64
#    define __ISWINDOWS__
#elif _WIN32
#    define __ISWINDOWS__
#elif __APPLE__
#    define __ISAPPLE__
#elif __linux || __unix || __posix
#    define __ISLINUX__
#elif __powerpc__
#    define __ISPOWERPC__
#else
#    pragma error
#endif

#endif
