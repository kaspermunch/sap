#ifndef CROSSPLATFORM_H
#define CROSSPLATFORM_H

#ifdef _WIN32

#define _MSC_VER

#include <direct.h>
#define GETCWD _getcwd
#define PATH_SEPERATOR "\\"

#else

#define GETCWD getcwd
#define PATH_SEPERATOR "/"

#endif

#endif
