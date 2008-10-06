#ifndef CROSSPLATFORM_H
#define CROSSPLATFORM_H

//#define WIN_PLATFORM

#if defined(WIN_PLATFORM)
#include <direct.h>
#define GETCWD _getcwd
#else
#define GETCWD getcwd
#endif

#endif
