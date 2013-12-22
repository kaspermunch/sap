#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#ifdef _WIN32
#include <io.h>
#endif


extern void initwrapIMa (void) { /* should be init + name of module as given in setup.py */
  char *dummy = NULL; // dummy body
}

int run_main(int argc, char **argv);

int runprogram (int argc, char **argv, char *outputBaseName) {

  int backupStdout, backupStderr, fdOut, fdErr;
  mode_t mode;

  char *stdoutFile = NULL;
  char *stderrFile = NULL;

  stdoutFile = (char *) calloc(strlen(outputBaseName) + 5, sizeof(char));
  sprintf(stdoutFile,"%s.out", outputBaseName);  
  stderrFile = (char *) calloc(strlen(outputBaseName) + 5, sizeof(char));
  sprintf(stderrFile,"%s.err", outputBaseName);  
  
  
#ifdef _WIN32
  /* duplicate file stdout and stderr descriptors */
  backupStdout = _dup(1);
  backupStderr = _dup(2);
  
  /* open files to write to */
  mode = _S_IREAD | _S_IWRITE;
  fdOut = _open(stdoutFile, _O_WRONLY | _O_CREAT, _S_IREAD | _S_IWRITE);
  fdErr = _open(stderrFile, _O_WRONLY | _O_CREAT, _S_IREAD | _S_IWRITE);
  
  /* duplicate it into 1 and 2 (the place of stdout and stderr) */
  _dup2(fdOut, 1);
  _dup2(fdErr, 2);
#else
  /* same thing for possix */
  backupStdout = dup(1);
  backupStderr = dup(2);
  
  mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
  fdOut = open(stdoutFile, O_WRONLY | O_CREAT | O_TRUNC, mode);
  fdErr = open(stderrFile, O_WRONLY | O_CREAT | O_TRUNC, mode);
  
  dup2(fdOut, 1);
  dup2(fdErr, 2);
#endif
  /* stdout using printf now goes to the file */
  
  // fflush(NULL); // flushes all buffers
  
  /* run function that includes calls to printf */
  run_main(argc, argv);
  
#ifdef _WIN32
  /* duplicate the original stdout and stderr descriptor into 1 and 2 */   
  _dup2(backupStdout, 1);
  _dup2(backupStderr, 2);
#else
  /* same thing for possix */
  dup2(backupStdout, 1);
  dup2(backupStderr, 2);
#endif
  
  /* printf and python print now goes to stdout again */
  
  return 0;
}
