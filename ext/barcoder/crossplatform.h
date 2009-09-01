#ifndef CROSSPLATFORM_H
#define CROSSPLATFORM_H

#ifdef _WIN32

#include <direct.h>
#define GETCWD _getcwd
#define PATH_SEPERATOR "\\"

#else

#define GETCWD getcwd
#define PATH_SEPERATOR "/"

#endif

#endif
