#ifndef CROSSPLATFORM_SHARED_PTR
#define CROSSPLATFORM_SHARED_PTR

#include "PlatformDetermination.h"

// Based on the platform and compiler, include the necessary header to give access to std::tr1::shared_ptr

#if defined(__ISLINUX__)
#include <tr1/memory>
#elif defined(__ISAPPLE__)
#include <tr1/memory>
#elif defined(__ISWINDOWS__) && defined(__MINGW32__)
#include <memory>
#else
#pragma error
#endif

using namespace std::tr1;

#endif
