

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "NW.h"

int  align( const char *seq_1, const char *seq_2, double gap_init, double gap_ext, double gap_flank, char **seq_1_al_ptr, char **seq_2_al_ptr) {

  int             L1, L2, j ;
  int             L1_al, L2_al ;
  int     i, k = 0, x, y ;
  int     fU, fD, fL ;
        char    ptr, nuc ;

        L1 = strlen( seq_1 ) ;
        L2 = strlen( seq_2 ) ;

	int **F = (int **) calloc(sizeof(int *), L2+1);
	for (i=0; i<L2+1; i++)
	  F[i] = (int *) calloc(sizeof(int), L1+1);

	int **traceback = (int **) calloc(sizeof(int *), L2+1);
	for (i=0; i<L2+1; i++)
	  traceback[i] = (int *) calloc(sizeof(int), L1+1);

	char *seq_1_al = (char *) calloc(L1+L2, sizeof(char));
	char *seq_2_al = (char *) calloc(L1+L2, sizeof(char));

	//        int     gap_init = 5 ;                   /* gap penalty */
        //int     gap_ext = 2 ;                    /* gap penalty */
        //int     gap_flank = 1 ;                    /* gap penalty */

/*         const int  s[ 4 ][ 4 ] = { {  5, -4, -4, -4 },
 *                                    { -4,  5, -4, -4 },
 *                                    { -4, -4,  5, -4 },
 *                                    { -4, -4, -4,  5 } } ;
 */
	/* identity matrix as used in clustalw2 */
        const int  s[ 5 ][ 5 ] = { { 10,  0,  0,  0, 0 },
                                   {  0, 10,  0,  0, 0 },
                                   {  0,  0, 10,  0, 0 },
                                   {  0,  0,  0, 10, 0 },
                                   {  0,  0,  0,  0, 0 } } ;


	x = 'N'; y = 'N'; /* just to avoid warnings about uninitialization */

        F[ 0 ][ 0 ] =  0 ;
        traceback[ 0 ][ 0 ] =  'n' ;

        for( j = 1; j <= L1; j++ )
        {
	        F[ 0 ][ j ] =   -j * gap_flank ;
                traceback[ 0 ][ j ] =  '-' ;
        }

        for( i = 1; i <= L2; i++ )
        {
                F[ i ][ 0 ] =   -i * gap_flank ;
                traceback[ i ][ 0 ] =  '|' ;
        }

        for( i = 1; i <= L2; i++ )
        {
                for( j = 1; j <= L1; j++ )
                {
                        nuc = seq_1[ j-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  x = 0 ;  break ;
                                case 'C':  x = 1 ;  break ;
                                case 'G':  x = 2 ;  break ;
            			case 'T':  x = 3 ;  break ;
                                default:  x = 4 ;
                        }

                        nuc = seq_2[ i-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  y = 0 ;  break ;
                                case 'C':  y = 1 ;  break ;
                                case 'G':  y = 2 ;  break ;
                                case 'T':  y = 3 ;  break ;
                                default:  x = 4 ;
                        }

			if (j == L1) {
			  fU = F[ i-1 ][ j ] - gap_flank ;
			} else if (traceback[ i-1 ][ j ] == '|') {
			  fU = F[ i-1 ][ j ] - gap_ext ;
			} else {
			  fU = F[ i-1 ][ j ] - gap_init ;
			}

                        fD = F[ i-1 ][ j-1 ] + s[ x ][ y ] ;

			if (i == L2) {
			  fL = F[ i ][ j-1 ] - gap_flank ;
			} else if (traceback[ i ][ j-1 ] == '-') {
			  fL = F[ i ][ j-1 ] - gap_ext ;
			} else {
			  fL = F[ i ][ j-1 ] - gap_init ;
			}
			
                        F[ i ][ j ] = max( fU, fD, fL, &ptr ) ;

                        traceback[ i ][ j ] =  ptr ;
                }
        }
        i-- ; j-- ;

/* 	printf("%d\n", F[L2][L1]);
 */

        while( i > 0 || j > 0 )
        {
                switch( traceback[ i ][ j ] )
                {
                        case '|' :  seq_1_al[ k ] = '-' ; 
                                    seq_2_al[ k ] = seq_2[ i-1 ] ; 
                                    i-- ;
                                    break ;

                        case '\\':  seq_1_al[ k ] = seq_1[ j-1 ] ; 
                                    seq_2_al[ k ] = seq_2[ i-1 ] ; 
                                    i-- ;  j-- ;
                                    break ;

                        case '-' :  seq_1_al[ k ] = seq_1[ j-1 ] ; 
                                    seq_2_al[ k ] = '-' ; 
                                    j-- ;
                }
                k++ ;
        }
        
        seq_1_al[ k ] = '\0' ; 
        seq_2_al[ k ] = '\0' ; 

        L1_al = strlen( seq_1_al ) ;
        L2_al = strlen( seq_2_al ) ;


/*          printf("\n\n") ;
 * 
 * 
 * 
 *         for( i = 0; i <= L2; i++ )
 *         {
 *                 for( j = 0; j <= L1; j++ )
 *                 {
 *                         printf("%-4d ", F[ i ][ j ] ) ;
 *                 }
 *                 printf("\n") ;
 *         }
 * 
 * 
 *         fprintf( stdout, "    " ) ;
 *         
 *         for( j = 0; j < L1; j++ )
 *         {
 *                 fprintf( stdout, "%c ", seq_1[ j ] ) ;
 *         }
 *         
 *         fprintf( stdout, "\n" ) ;
 *         fprintf( stdout, "  " ) ;
 * 
 *         for( i = 0; i <= L2; i++ )
 *         {
 *                 if( i > 0 )
 *                 {
 *                         fprintf( stdout, "%c ", seq_2[ i-1 ] ) ;
 *                 }
 * 
 *                 for( j = 0; j <= L1; j++ )
 *                 {
 *                         printf("%c ", traceback[ i ][ j ] ) ;
 *                 }
 *                 printf("\n") ;
 *         }
 *         printf("\n") ;
 * 
 * 
 * 
 *         for( j = L1_al - 1; j >= 0; j-- )
 *         {
 *                 fprintf( stdout, "%c ", seq_1_al[ j ] ) ;
 *         }
 *         fprintf( stdout, "\n" ) ;
 * 
 *         for( j = L2_al - 1; j >= 0; j-- )
 *         {
 *                 fprintf( stdout, "%c ", seq_2_al[ j ] ) ;
 *         }
 *         fprintf( stdout, "\n\n" ) ;
 */


	//printf("rev\n"); fflush(NULL);
	
	// reverse the strings:
 	rev_string(seq_1_al, L1_al);
 	rev_string(seq_2_al, L2_al);

	*seq_1_al_ptr = seq_1_al;
	*seq_2_al_ptr = seq_2_al;


	for (i=0; i<L2+1; i++) {
	  free(F[i]);
	  free(traceback[i]);
	}
	free(F);
	free(traceback);
	// Try and print these output strings from the python code

        return  1 ;
}


void rev_string(char *s, int n)
{
  int i, m;
  char t;
  m = floor(n/2.0);
  for (i=0; i<m; i++)
    {
      t = s[i];
      s[i] = s[n-1-i];
      s[n-1-i] = t;
    }
}


int  max(int f1, int f2, int f3, char *ptr )
/*         int    f1 ;
 *         int    f2 ;
 *         int    f3 ;
 *         char * ptr ;
 */
{
        int  max = 0 ;

        if( f1 >= f2 && f1 >= f3 )  
        {
                max = f1 ;
                *ptr = '|' ;
        }
        else if( f2 > f3 )              
        {
                max = f2 ;
                *ptr = '\\' ;
        }
        else
        {
                max = f3 ;
                *ptr = '-' ;
        }
        
        return  max ;   
}


int main(void)
{
  int i;

  for (i=0; i<1000; i++) 
    {
      const char str1[] = "AGTCCGGATCGAGAAGTCCGGATCGAGACCAGTAGTCCGGATCGAGACCAGTAGTCCGGATCGAGACCAGTAGTCCGGATCGAG";
      const char str2[] = "AGTCCGGATCGAGACCAGTAGTCCGGATCGAGACCAGTAGTCCGGATCCAGTAGTCCGGATCGAGACCAGTAGTCCGGATCGAG";
      
      
      char *out1;
      char *out2;
      
      align(str1, str2, 5, 2, 1, &out1, &out2);
      

      printf("out1: %s\n", out1);
      printf("out2: %s\n", out2);
    }

  return 1;
}
