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



/* In your header files, wherever you want an interface or API made public outside the current DSO, place __attribute__ ((visibility ("default"))) in struct, class and function declarations you wish to make public (it's easier if you define a macro as this). You don't need to specify it in the definition. Then, alter your make system to pass -fvisibility=hidden to each call of GCC compiling a source file. If you are throwing exceptions across shared object boundaries see the section "Problems with C++ exceptions" below. Use nm -C -D on the outputted DSO to compare before and after to see the difference it makes.
 */


defined _WIN32 || defined __CYGWIN__
  #ifdef __GNUC__
     #define DLL_PUBLIC __attribute__ ((dllexport))
  #else
     #define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
  #endif
  #define DLL_LOCAL
#else
  #if __GNUC__ >= 4
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define DLL_PUBLIC
    #define DLL_LOCAL
  #endif
#endif

/* #if defined _WIN32 || defined __CYGWIN__
 *   #ifdef BUILDING_DLL
 *     #ifdef __GNUC__
 *       #define DLL_PUBLIC __attribute__ ((dllexport))
 *     #else
 *       #define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
 *     #endif
 *   #else
 *     #ifdef __GNUC__
 *       #define DLL_PUBLIC __attribute__ ((dllimport))
 *     #else
 *       #define DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
 *     #endif
 *   #endif
 *   #define DLL_LOCAL
 * #else
 *   #if __GNUC__ >= 4
 *     #define DLL_PUBLIC __attribute__ ((visibility ("default")))
 *     #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
 *   #else
 *     #define DLL_PUBLIC
 *     #define DLL_LOCAL
 *   #endif
 * #endif
 */

/* extern "C" DLL_PUBLIC void function(int a);
 * class DLL_PUBLIC SomeClass
 * {
 *    int c;
 *    DLL_LOCAL void privateMethod();  // Only for use within this DSO
 * public:
 *    Person(int _c) : c(_c) { }
 *    static void foo(int a);
 * };
 */
