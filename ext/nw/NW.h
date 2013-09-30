#ifndef NJ_H
#define NJ_H

int  align( const char *seq_1, const char *seq_2,  double gap_init, double gap_ext, double gap_flank, char **seq_1_al_ptr, char **seq_2_al_ptr);

int   max             ( int, int, int, char * );

void rev_string(char *s, int n);


#endif
