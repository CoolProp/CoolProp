#ifndef PLATFORMDETERMINATION_H
#define PLATFORMDETERMINATION_H

#if defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WIN64__)
#  define __ISWINDOWS__
#elif __APPLE__
#  define __ISAPPLE__
#elif __linux
#  define __ISLINUX__
#endif

#endif