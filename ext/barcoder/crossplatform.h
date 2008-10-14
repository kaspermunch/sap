#ifndef CROSSPLATFORM_H
#define CROSSPLATFORM_H

//#define WIN_PLATFORM

#if defined(WIN_PLATFORM)

#include <direct.h>
#define GETCWD _getcwd

#define PATH_SEPERATOR "\\"

#else

#define GETCWD getcwd
#define PATH_SEPERATOR "/"

#endif

#endif
