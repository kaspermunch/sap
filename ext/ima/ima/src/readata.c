/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/


#undef GLOBVARS
#include <ctype.h>
#include "imamp.h"
#include "updateassignment.h"


/* read in the data */

/*********** LOCAL STUFF **********/
#define DATAFILEMAXLINELENGTH  300
static int infilelines;
double pi[MAXLOCI][4];          // used here and in initialize.c


/* prototypes of local functions*/
static int findsegsites (FILE * infile, int li, int numbases, int initseg[],
                         int MODEL, int **numsitesIS);
static void elimfrom (int li, int site);
static int seqid (int li, int i, int j);
static void sortseq (int li);
static void eliminategaps (int li);
static void readseqHKY (FILE * infile, int li);
static void readseqIS (FILE * infile, int li, int MODEL, int **numsitesIS);
static void readseqSW (FILE * infile, int li);
static void readseqSWP (FILE * infile, int li);
static void parse_locus_info (int li, int *uinext, char *cc, int *fpstri,
                              char fpstr[]);
static void skip_datafile_toplines (FILE * infile);

/* For IS model, identify and count segregating sites */
int
findsegsites (FILE * infile, int li, int numbases, int initseg[],
              int MODEL, int **numsitesIS)
{
  char c, *zeroc, *altc, **zerocpop, **altcpop;
  int **initsegpop;
  int i, pop, j, v, totseg = 0, A;
  int b, e, firstpopline;
  int np;

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    np = npops + 1;
  }
  else
  {
    np = npops;
  }

  zeroc = malloc (numbases * (sizeof (char)));
  altc = malloc (numbases * (sizeof (char)));

  zerocpop = malloc (np * (sizeof (char *)));
  altcpop = malloc (np * (sizeof (char *)));
  initsegpop = malloc (np * (sizeof (int *)));
  for (i = 0; i < np; i++)
  {
    zerocpop[i] = malloc (numbases * (sizeof (char)));
    altcpop[i] = malloc (numbases * (sizeof (char)));
    initsegpop[i] = calloc ((size_t) numbases, (sizeof (int)));
  }
  b = 0;
  e = L[li].samppop[0] - 1;
  pop = 0;
  firstpopline = 1;
  for (i = 0; i < L[li].numgenes; i++)
  {
    while (i > e)
    {
      pop++;
      b = e + 1;
      if (pop == npops)
      {
        assert (assignmentoptions[POPULATIONASSIGNMENT] == 1);
        e = b + L[li].numgenesunknown - 1;
      }
      else
      {
        e = b + L[li].samppop[pop] - 1;
      }
      firstpopline = 1;
    }
    /* assumes that the first sequence of 10 characters that does not 
     * contain a carriage return is the species name */
    for (v = 0; v < LENGENENAME; v++)
    {
      if ((char) fgetc (infile) == '\n')
        v = 0;
    }
    if (MODEL == JOINT_IS_SW)   // just read through theses 
    {
      for (j = 0; j < L[li].nAlinked; j++)
        fscanf (infile, "%d", &A);
    }
    j = 0;
    while ((c = (char) tolower ((fgetc (infile)))) != '\n'
           && j < L[li].numbases)

    {
      if (c != ' ')
      {
        if (i == 0)
        {
          zeroc[j] = c;
          altc[j] = ' ';
          if ((c != 'a' && c != 'c') && (c != 't' && c != 'g'))
            L[li].badsite[j] = 1;
        }
        else
        {
          if ((c != 'a' && c != 'c') && (c != 't' && c != 'g'))
            L[li].badsite[j] = 1;
          if (L[li].badsite[j] == 0 && c != zeroc[j])
          {
            if (altc[j] == ' ')
            {
              altc[j] = c;
            }
            else
            {
              if (c != altc[j])
                L[li].badsite[j] = 1;
            }
            initseg[j] = 1;
          }
        }
        if (firstpopline)
        {
          zerocpop[pop][j] = c;
          altcpop[pop][j] = ' ';
        }
        else
        {
          if (L[li].badsite[j] == 0 && c != zerocpop[pop][j])
          {
            if (altcpop[pop][j] == ' ')
              altcpop[pop][j] = c;
            initsegpop[pop][j] = 1;
          }
        }
        j++;
      }
    }
    firstpopline = 0;
  }

  printf ("Locus %d seg. sites: ", li);
  for (i = 0; i < numbases; i++)
  {
    if (L[li].badsite[i] == 1)
    {
      initseg[i] = 0;
      for (pop = 0; pop < np; pop++)
        initsegpop[pop][i] = 0;
    }
    if (initseg[i] == 1)
    {
      printf ("%d ", i);
      totseg++;
    }
    for (pop = 0; pop < np; pop++)
      if (initsegpop[pop][i] == 1)
        numsitesIS[li][pop]++;
  }
  printf ("\n");
  XFREE (zeroc);
  XFREE (altc);
  for (i = 0; i < np; i++)

  {
    XFREE (zerocpop[i]);
    XFREE (altcpop[i]);
    XFREE (initsegpop[i]);
  }
  XFREE (zerocpop);
  XFREE (altcpop);
  XFREE (initsegpop);
  return totseg;
}                               /* findsegsites */

void
elimfrom (int li, int site)     // called when reading in HKY data
{
  int i, j;
  for (i = 0; i < L[li].numgenes; i++)
  {
    for (j = site; j < L[li].numsites - 1; j++)
      L[li].seq[i][j] = L[li].seq[i][j + 1];
  }
}
int
seqid (int li, int i, int j)
{
  int n;
  for (n = 0; n < L[li].numgenes; n++)
  {
    if (L[li].seq[n][i] != L[li].seq[n][j])
      return 0;
  }
  return 1;
}

void
sortseq (int li)                // called when reading in HKY data
{
  int i, j;
  L[li].mult = malloc ((L[li].numsites) * (sizeof (int)));
  for (i = 0; i < L[li].numsites; i++)
    L[li].mult[i] = 1;
  for (i = 0; i < L[li].numsites; i++)
  {
    for (j = i + 1; j < L[li].numsites; j++)
    {
      /*printf("compares codon %i and %i\n",i+1,j+1); */
      if (seqid (li, i, j) == 1)
      {
        /*      printf("eliminates codon %i\n",j+1); */
        elimfrom (li, j);
        j--;
        L[li].numsites--;
        L[li].mult[i]++;
      }
    }
  }
}

void
eliminategaps (int li)          /* called when reading in HKY data */
{
  int i, j;
  for (i = 0; i < L[li].numsites; i++)
  {
    for (j = 0; j < L[li].numgenes; j++)
    {
      if (L[li].seq[j][i] == -1)
      {
        elimfrom (li, i);
        i--;
        L[li].numsites--;
      }
    }
  }
  printf ("Locus %i: %s, HKY  model %i sites after elimination of gaps \n",
          li, L[li].name, L[li].numsites);
}

void
readseqHKY (FILE * infile, int li)
{
  char gName[LENGENENAME + 1];
  int i, j, v, k = 0;
  char c;
  double PIstandard;

  L[li].umodel[0] = HKY;
  L[li].numsites = L[li].numbases;
  for (i = 0; i < 4; i++)
    pi[li][i] = 0.0;
  L[li].seq = malloc (L[li].numgenes * (sizeof (int *)));
  for (i = 0; i < L[li].numgenes; i++)
    L[li].seq[i] = malloc (L[li].numsites * (sizeof (int)));
  do
  {
    for (i = 0; i < L[li].numgenes; i++)
    {
      /* assumes that the first sequence of 10 characters that does not 
       * contain a carriage return is the species name */
      for (v = 0; v < LENGENENAME; v++)
      {
// is this trapping of eol really necessary ??? 
        if ((c = (char) fgetc (infile)) == '\n')
          v = -1;
        if (isprint (c) != 0)
        {
          gName[v] = c;
        }
        else
        {
          gName[v] = '\0';
          IM_err (IMERR_GENENAME,
                  "locus (%d) %d-th gene name, partial name %s", li, i,
                  gName);
        }
      }
      gName[LENGENENAME] = '\0';
      strcpy (L[li].gNames[i], gName);

      j = k;
      while ((c = (char) fgetc (infile)) != '\n')
      {
        if ((c != ' ') && (c != '\t'))
        {
          if (i >= L[li].numgenes || j >= L[li].numsites)
          {
            IM_err(IMERR_DATAREADOVERRUN,"HKY data locus %d gene# %d site# %d",li,i,j);
          }
          if (c == 'a' || c == 'A')
          {
            L[li].seq[i][j] = 0;
            pi[li][0]++;
          }
          else if (c == 'c' || c == 'C')
          {
            L[li].seq[i][j] = 1;
            pi[li][1]++;
          }
          else if (c == 'g' || c == 'G')
          {
            L[li].seq[i][j] = 2;
            pi[li][2]++;;
          }
          else if (c == 't' || c == 'u' || c == 'T' || c == 'U')
          {
            L[li].seq[i][j] = 3;
            pi[li][3]++;
          }
          else if (c == 'n' || c == '-' || c == 'N' || c == '.')
          {
            L[li].seq[i][j] = -1;
          }
          else
          {
            IM_err(IMERR_DATAERROR,"BAD BASE in locus %d species %i base %i: %c",li, i + 1, j + 1, c);
          }
          j++;
          if (i == (L[li].numgenes - 1))
            k++;
        }
      }
    }
  }
  while (k < L[li].numsites);
  eliminategaps (li);
  L[li].totsites=L[li].numsites;
  sortseq (li);
  PIstandard = 0.0;
  for (i = 0; i < 4; i++)
    PIstandard += pi[li][i];
  for (i = 0; i < 4; i++)
  {
    pi[li][i] = pi[li][i] / PIstandard;

    /*pi[li][i]=0.25; JC only */
  }
  infilelines += L[li].numgenes;
  return;
}                               /* readseqHKY */

void
readseqIS (FILE * infile, int li, int MODEL, int **numsitesIS)
{
  char gName[LENGENENAME + 1];
  int i, j, ai, v, k, sitek, sitej;
  int *initseg;
  char c, *refseq;
  int numinvar;

  L[li].umodel[0] = INFINITESITES;
  initseg = calloc ((size_t) L[li].numbases, (sizeof (int)));
  L[li].badsite = calloc ((size_t) L[li].numbases, (sizeof (int)));
  if (L[li].numbases > 0)
  {
    L[li].numsites =
      findsegsites (infile, li, L[li].numbases, initseg, MODEL, numsitesIS);
    if (L[li].model == INFINITESITES)
      printf ("Locus %i : %s, Infinite Sites model, %i sites are segregating\n", li, L[li].name, L[li].numsites);
    if (L[li].model == JOINT_IS_SW)
      printf ("Locus %i : %s, Joint Infinite Sites/Stepwise model, %i sites are segregating\n", li, L[li].name, L[li].numsites);
    rewind (infile);
    for (i = 0; i < infilelines; i++)
      while ((c = ((char) fgetc (infile))) != '\n');
  }
  else
  {
    L[li].numsites = 0;
  }
  numinvar = L[li].numbases - L[li].numsites;
  refseq = malloc (L[li].numbases * (sizeof (char)));
  L[li].seq = malloc (L[li].numgenes * (sizeof (int *)));
  for (i = 0; i < L[li].numgenes; i++)
    L[li].seq[i] = malloc ((L[li].numsites) * (sizeof (int)));

  assert (L[li].numbases > 0);
  if (MODEL == JOINT_IS_SW)
  {
    L[li].A = malloc (L[li].nlinked * sizeof (int *));
    for (ai = 1; ai < L[li].nlinked; ai++)
    {
      L[li].umodel[ai] = STEPWISE;
      L[li].A[ai] = malloc (L[li].numgenes * sizeof (int));
      L[li].maxA[ai] = -1;
      L[li].minA[ai] = 10000;   // something large
    }
  }
  if (L[li].numbases > 0)
  {
    k = 0;
    sitek = 0;

    do
    {
      for (i = 0; i < L[li].numgenes; i++)
      {
        /*assumes that the first sequence of 10 characters that
         * does not contain a carriage return is the species name
         * */
        for (v = 0; v < LENGENENAME; v++)
        {
          if ((c = (char) fgetc (infile)) == '\n')
            v = -1;
          if (isprint (c) != 0)
          {
            gName[v] = c;
          }
          else
          {
            gName[v] = '\0';
            IM_err (IMERR_GENENAME,
                    "locus (%d) %d-th gene name, partial name %s", li, i,
                    gName);
          }
        }
        gName[LENGENENAME] = '\0';
        strcpy (L[li].gNames[i], gName);

        if (MODEL == JOINT_IS_SW)
        {
          for (ai = 1; ai < L[li].nlinked; ai++)
          {
            fscanf (infile, "%d", &L[li].A[ai][i]);
            if (L[li].A[ai][i] > L[li].maxA[ai])
              L[li].maxA[ai] = L[li].A[ai][i];
            if (L[li].A[ai][i] < L[li].minA[ai])
              L[li].minA[ai] = L[li].A[ai][i];
          }
        }
        j = k;
        sitej = sitek;
        while ((c = (char) tolower ((fgetc (infile)))) != '\n'
               && j < L[li].numbases)
        {
          if (isdigit (c))
            IM_err(IMERR_DATAERROR,"locus %d, formatting of input file causes wrong lines to be read as data",li);
          if (c != ' ')
          {
            if (initseg[j] == 1)
            {
              if (i == 0)
              {
                refseq[j] = c;
                L[li].seq[0][sitej] = 0;
              }
              else if (c == refseq[j])
              {
                L[li].seq[i][sitej] = 0;
              }
              else
              {
                L[li].seq[i][sitej] = 1;
              }
              sitej++;
              if (i == (L[li].numgenes - 1))
                sitek++;
            }
            j++;
            if (i == (L[li].numgenes - 1))
              k++;
          }
        }
      }
    } while (k < L[li].numbases);
  }
  infilelines += L[li].numgenes;
  XFREE (initseg);
  XFREE (refseq);

#ifdef COMMENT
  assert (_CrtCheckMemory ());
#endif /* COMMENT */
}                               /* readseqIS */


/* read in the allele sizes for a data set under the STEPWISE model */
void
readseqSW (FILE * infile, int li)
{
  int i, j, ai, tempA, numA;
  char tempname[10], *c;
  char ch;
  char textline[301];
  /* int maxlinelength = DATAFILEMAXLINELENGTH; */
  char gName[11];
  int v;

  /* This for-loop may be placed here out of the main for-loop below */
  L[li].A = malloc (L[li].nlinked * sizeof (int *));
  for (ai = 0; ai < L[li].nlinked; ai++)
  {
    L[li].umodel[ai] = STEPWISE;
    L[li].A[ai] = malloc (L[li].numgenes * sizeof (int));
    L[li].maxA[ai] = -1;
    L[li].minA[ai] = 10000;     // something large
  }
  /* the main for-loop */
  for (i = 0; i < L[li].numgenes; i++)
  {
     //fgets (c, maxlinelength, infile);
    if (jheyoptions[SWINPUTOPTION])     // alernate input option  only works if there is only one SW portion
    {
      assert (0);
      c = &textline[0];
      sscanf (c, "%s ", &tempname[0]);
      strncpy (gName, textline, 10);
      gName[10] = '\0';
      strcpy (L[li].gNames[i], gName);

      c = nextwhite (c);
      sscanf (c, "%d %d", &tempA, &numA);
      if (tempA > L[li].maxA[ai])
        L[li].maxA[ai] = tempA;
      if (tempA < L[li].minA[ai])
        L[li].minA[ai] = tempA;

      for (j = 0; j < numA; j++)
      {
        L[li].A[0][i + j] = tempA;
      }
      i += j - 1;
    }
    else
    {
      /* We need this for handling 10-character gene names that are followed
       * allele numbers without any space */
      for (v = 0; v < LENGENENAME; v++)
      {
        if ((ch = (char) fgetc (infile)) == '\n')
          v = -1;
        if (isprint (ch) != 0)
        {
          gName[v] = ch;
        }
        else
        {
          gName[v] = '\0';
          IM_err (IMERR_GENENAME,
                  "locus (%d) %d-th gene name, partial name %s", li, i,
                  gName);
        }
      }
      gName[LENGENENAME] = '\0';

      //sscanf (textline, "%s ", &tempname[0]);
      //strncpy (gName, textline, 10);
      //gName[10] = '\0';
      strcpy (L[li].gNames[i], gName);
      //ch = nextwhite (ch);
      for (ai = 0; ai < L[li].nAlinked; ai++)
      {
        fscanf (infile, "%d", &(L[li].A[ai][i]));
        if (L[li].A[ai][i] > L[li].maxA[ai])
          L[li].maxA[ai] = L[li].A[ai][i];
        if (L[li].A[ai][i] < L[li].minA[ai])
          L[li].minA[ai] = L[li].A[ai][i];
        //ch = nextwhite (ch);
      }
      skip_a_line (infile);
      infilelines++;;
    }
  }
/*
  for (ai = 0; ai < L[li].nlinked; ai++)
  {
    if (L[li].maxA[ai] - L[li].minA[ai])
  }
*/
  printf ("Locus %i : %s, Stepwise Mutation Model\n", li, L[li].name);
}                               /* readseqSW */

/* read in the allele sizes for a data set under the STEPWISE model */
void
readseqSWP (FILE * infile, int li)
{
  int i, j, ai, tempA, numA;
  char tempname[10], *c;
  char ch;
  char textline[301];
  /* int maxlinelength = DATAFILEMAXLINELENGTH; */
  char gName[11];
  int v;

  /* This for-loop may be placed here out of the main for-loop below */
  L[li].A = malloc (L[li].nlinked * sizeof (int *));
  for (ai = 0; ai < L[li].nlinked; ai++)
  {
    L[li].umodel[ai] = STEPWISE;
    L[li].A[ai] = malloc (L[li].numgenes * sizeof (int));
    L[li].maxA[ai] = -1;
    L[li].minA[ai] = 10000;     // something large
  }
  /* the main for-loop */
  for (i = 0; i < L[li].numgenes; i++)
  {
     //fgets (c, maxlinelength, infile);
    if (jheyoptions[SWINPUTOPTION])     // alernate input option  only works if there is only one SW portion
    {
      assert (0);
      c = &textline[0];
      sscanf (c, "%s ", &tempname[0]);
      strncpy (gName, textline, 10);
      gName[10] = '\0';
      strcpy (L[li].gNames[i], gName);

      c = nextwhite (c);
      sscanf (c, "%d %d", &tempA, &numA);
      if (tempA > L[li].maxA[ai])
        L[li].maxA[ai] = tempA;
      if (tempA < L[li].minA[ai])
        L[li].minA[ai] = tempA;

      for (j = 0; j < numA; j++)
      {
        L[li].A[0][i + j] = tempA;
      }
      i += j - 1;
    }
    else
    {
      /* We need this for handling 10-character gene names that are followed
       * allele numbers without any space */
      for (v = 0; v < LENGENENAME; v++)
      {
        if ((ch = (char) fgetc (infile)) == '\n')
          v = -1;
        if (isprint (ch) != 0)
        {
          gName[v] = ch;
        }
        else
        {
          gName[v] = '\0';
          IM_err (IMERR_GENENAME,
                  "locus (%d) %d-th gene name, partial name %s", li, i,
                  gName);
        }
      }
      gName[LENGENENAME] = '\0';

      //sscanf (textline, "%s ", &tempname[0]);
      //strncpy (gName, textline, 10);
      //gName[10] = '\0';
      strcpy (L[li].gNames[i], gName);
      //ch = nextwhite (ch);
      for (ai = 0; ai < L[li].nAlinked; ai++)
      {
        fscanf (infile, "%d", &(L[li].A[ai][i]));
        if (L[li].A[ai][i] > L[li].maxA[ai])
          L[li].maxA[ai] = L[li].A[ai][i];
        if (L[li].A[ai][i] < L[li].minA[ai])
          L[li].minA[ai] = L[li].A[ai][i];
        //ch = nextwhite (ch);
      }
      skip_a_line (infile);
      infilelines++;;
    }
  }
/*
  for (ai = 0; ai < L[li].nlinked; ai++)
  {
    if (L[li].maxA[ai] - L[li].minA[ai])
  }
*/
  printf ("Locus %i : %s, Stepwise Mutation Model\n", li, L[li].name);
}                               /* readseqSW */



void
parse_locus_info (int li, int *uinext, char *cc, int *fpstri, char fpstr[])
{
  int ui, i;

  ui = *uinext;
  sscanf (cc, "%s", L[li].name);
  cc = nextnonspace (cc);
  L[li].numgenes = 0;
  L[li].numgenesknown = 0;
  for (i = 0; i < npops; i++)
  {
    sscanf (cc, "%d", &L[li].samppop[i]);
    cc = nextnonspace (cc);
    L[li].numgenes += L[li].samppop[i];
    L[li].numgenesknown += L[li].samppop[i];
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 0)
  {
    L[li].numgenesunknown = 0;
  }
  // set ghost pop sample add 0
  /*
  if (modeloptions[ADDGHOSTPOP])
  {
    assert (assignmentoptions[POPULATIONASSIGNMENT] == 0);
    if (npops + 1 > MAXPOPS)
    {
      IM_err (IMERR_INPUTFILEINVALID, 
              "ghost population makes number of population (%d) greater than MAXPOPS [%d]", 
              npops, MAXPOPS);
    }
    L[li].samppop[npops] = 0;
    npops++;
    numtreepops += 2;
  }*/

  SP "%d\t%s", li, L[li].name);
  for (i = 0; i < npops; i++)
  {
    SP "\t%3d", L[li].samppop[i]);
  }
  sscanf (cc, "%d", &L[li].numbases);
  cc = nextnonspace (cc);
  L[li].nlinked = 0;
  L[li].nAlinked = 0;
  switch (toupper (cc[0]))
  {
  case 'I':
    L[li].model = INFINITESITES;
    SP "\tIS");
    *uinext = ui + 1;
    L[li].nlinked = 1;
    break;
  case 'H':
    L[li].model = HKY;
    SP "\tHKY");
    *uinext = ui + 1;
    L[li].nlinked = 1;
    break;
  case 'S':
    L[li].model = STEPWISE;
    if (isdigit (cc[1]))
    {
      L[li].nAlinked = atoi (&cc[1]);
      L[li].nlinked = L[li].nAlinked;
      SP "\tSW_M");
    }
    else
    {
      /* These two lines have been absent. */
      L[li].nAlinked = 1; 
      L[li].nlinked = 1;
      SP "\tSW");
    }
    *uinext = ui + L[li].nlinked;
    break;
  case 'J':
    L[li].model = JOINT_IS_SW;
    if (isdigit (cc[1]))
    {
      L[li].nAlinked = atoi (&cc[1]);
      SP "\tIS+SW_M");
    }
    else
    {
      L[li].nAlinked = 1;
      SP "\tIS+SW");
    }
    L[li].nlinked = L[li].nAlinked + 1;
    *uinext = ui + L[li].nlinked;
    break;
  default:
    L[li].model = INFINITESITES;
    L[li].nlinked = 1;
    SP "\tIS");
  }
  if (L[li].nlinked < 1 || L[li].nlinked > MAXLINKED)
  {
    IM_err (IMERR_INPUTFILEINVALID,
            "The number of linked is less than 1 or greater than %d",
            MAXLINKED);
  }
  for (ui = 0; ui < L[li].nlinked; ui++)
    L[li].uperyear_prior[ui].min = L[li].uperyear_prior[ui].max = 0;
  cc = nextnonspace (cc);
  ui = 0;

  /* For ASSIGNMENT */
  if (cc != NULL && cc[0] == 'A')
  {
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      L[li].numgenesunknown = atoi (&cc[1]);
      assert (L[li].numgenesunknown > 0);
    }
    else
    {
      IM_err (IMERR_INPUTFILEINVALID,
              "There must be at least one gene with unknown origin: %s",
              cc);
    }
    cc = nextnonspace (cc);
  }

  /* get inheritance scalar info, mutation rate info */
  if (cc != NULL) 
  {
    sscanf (cc, "%lf", &(L[li].hval));
    SP "\t%5.3lf", L[li].hval);
    cc = nextnonspace (cc);

    /* For ASSIGNMENT */
    if (cc != NULL && cc[0] == 'A')
    {
      if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
      {
        L[li].numgenesunknown = atoi (&cc[1]);
        assert (L[li].numgenesunknown > 0);
      }
      else
      {
        IM_err (IMERR_INPUTFILEINVALID,
                "There must be at least one gene with unknown origin: %s",
                cc);
      }
      cc = nextnonspace (cc);
    }

    i = 0;
    while (cc != NULL && i < L[li].nlinked)     //get mutation rate info from data line
    {
      /* For ASSIGNMENT */
      if (cc != NULL && cc[0] == 'A')
      {
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          L[li].numgenesunknown = atoi (&cc[1]);
          assert (L[li].numgenesunknown > 0);
          cc = nextnonspace (cc);
          break;
        }
        else
        {
          IM_err (IMERR_INPUTFILEINVALID,
                  "There must be at least one gene with unknown origin: %s",
                  cc);
        }
      }

      i++;
      // sscanf (cc, "%lf", &L[li].u[ui].uperyear.val);
      sscanf (cc, "%lf", &L[li].uperyear_vals[ui]);
      // SP "\t%lg", L[li].u[ui].uperyear.val);
      SP "\t%lg", L[li].uperyear_vals[ui]);
      counturateperyear++;
      cc = nextnonspace (cc);

      /* For ASSIGNMENT */
      if (cc != NULL && cc[0] == 'A')
      {
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          L[li].numgenesunknown = atoi (&cc[1]);
          assert (L[li].numgenesunknown > 0);
          cc = nextnonspace (cc);
          break;
        }
        else
        {
          IM_err (IMERR_INPUTFILEINVALID,
                  "There must be at least one gene with unknown origin: %s",
                  cc);
        }
      }

      if (cc && cc[0] == '(')   /* look for a mutation rate range in parentheses e.g. (0.03,0.05)  */
      {
        cc++;
        sscanf (cc, "%lf",
                //                &(L[li].u[ui].uperyear.pr.min));
                &(L[li].uperyear_prior[ui].min));
        while (cc[0] != ',')
          cc++;
        cc++;
        sscanf (cc, "%lf",
                //                &(L[li].u[ui].uperyear.pr.max));
                &(L[li].uperyear_prior[ui].max));
        while (cc[0] != ')')
          cc++;
        while (!isdigit (cc[0]) && cc[0] != '\0')
          cc++;
        if (cc[0] == '\0')
          cc = NULL;
        if (L[li].uperyear_prior[ui].min > L[li].uperyear_prior[ui].max
            || L[li].uperyear_prior[ui].min > L[li].uperyear_vals[ui]
            || L[li].uperyear_prior[ui].max < L[li].uperyear_vals[ui])
              IM_err(IMERR_MUTSCALARPRIORRANGEFAIL, "locus %d  mutation scalar prior range problem min:%lf max:%lf current val%lf",
                li, L[li].uperyear_prior[ui].min,L[li].uperyear_prior[ui].max,L[li].uperyear_vals[ui]);
        SP "\t(%lg - %lg)", L[li].uperyear_prior[ui].min, L[li].uperyear_prior[ui].max);
        countuprior++;
      }
      ui++;
    }
  }
  else
  {
    L[li].hval = 1;
    SP "\t%lf", L[li].hval);
  }
  SP "\n");

  nurates += L[li].nlinked;

  /* For ASSIGNMENT */
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    if (L[li].numgenesunknown < 1)
    {
      IM_err (IMERR_INPUTFILEINVALID,
              "There must be at least one gene with unknown origin: %d",
              L[li].numgenesunknown);
    }
    L[li].numgenes += L[li].numgenesunknown;
  }

  L[li].numlines = 2 * L[li].numgenes - 1;
  // LEVINE LIKELIHOOD: 2008-07-11
  for (i = 0; i < L[li].numgenes; i++)
  {
    L[li].pairs[i] = -1;
  }

  return;
}                               //parse_locus_info


void 
skip_datafile_toplines (FILE * infile)
{
  char textline[DATAFILEMAXLINELENGTH + 1];
  char ch;
  fgets (textline, DATAFILEMAXLINELENGTH, infile);
  ch = (char) getc (infile);
  while (ch == '#')
  {
    fgets (textline, DATAFILEMAXLINELENGTH, infile);
    ch = (char) getc (infile);
  }
  ungetc (ch, infile);
  fgets (textline, DATAFILEMAXLINELENGTH, infile);
  fgets (textline, DATAFILEMAXLINELENGTH, infile);
  fgets (textline, DATAFILEMAXLINELENGTH, infile);
  fgets (textline, DATAFILEMAXLINELENGTH, infile);
  return;
}                               //skip_datafile_toplines


/*** GLOBAL STUFF *****/
void read_datafile_top_lines (char infilename[], int *fpstri, char fpstr[],
                              char startpoptreestring[])
{
  char textline[DATAFILEMAXLINELENGTH + 1];
  char ch;
  char popnames[MAXPOPS][NAMELENGTH];
  int i;
  int poptreestring_given;
  FILE *infile;

  if ((infile = fopen (infilename, "r")) == NULL)
  {
    printf ("Error opening text file for reading [filename: %s]\n",
            infilename);
    IM_err(IMERR_READFILEOPENFAIL,"data file not found or can't be opened: %s",infilename);
  }

  fgets (textline, DATAFILEMAXLINELENGTH, infile);
  if (strlen (textline) > DATAFILEMAXLINELENGTH - 3)
    {
      IM_err (IMERR_INPUTFILEINVALID, 
              "Length of a line is limited upto %d", 
              DATAFILEMAXLINELENGTH - 2);
    }
  infilelines++;
  SP "\nText from input file: %s", textline);
  ch = (char) getc (infile);
  while (ch == '#')
  {
    infilelines++;
    fgets (textline, DATAFILEMAXLINELENGTH, infile);
    if (strlen (textline) > DATAFILEMAXLINELENGTH - 3)
      {
        IM_err (IMERR_INPUTFILEINVALID, 
                "Length of a line is limited upto %d", 
                DATAFILEMAXLINELENGTH - 2);
      }
    SP "Text from input file: %s", textline);
    ch = (char) getc (infile);
  }
  ungetc (ch, infile);
  SP "\n");
  fscanf (infile, "%d\n", &npops);
  if (npops < 1)
  {
    IM_err (IMERR_INPUTFILEINVALID, 
            "number of population (%d) must be positive", 
            npops);
  }
  else if (npops > MAXPOPS)
  {
    IM_err (IMERR_INPUTFILEINVALID, 
            "number of population (%d) is greater than MAXPOPS [%d]", 
            npops, MAXPOPS);
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    numtreepops = npops + 1;
  }
  else
  {
    numtreepops = 2 * npops - 1;
  }

  for (i = 0; i < npops; i++)
    fscanf (infile, "%s ", popnames[i]);
  fscanf (infile, "\n");
  SP "Number of populations: %d \n", npops);
  SP "- Population Names - \n");
  for (i = 0; i < npops; i++)
    SP "Population %d : %s \n", i, popnames[i]);
  SP "\n");
  if (modeloptions[ADDGHOSTPOP])
    SP " **GHOST Population added to model **\n\n");
  fscanf (infile, "%s\n", startpoptreestring);

  // allow input files for two populations not to have a line with the population tree string 
  if (strlen (startpoptreestring) == 1)
  {
    poptreestring_given = 1;
    assert (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1);
  }
  else if (strlen (startpoptreestring) < 5)     // not a proper tree string.  maybe no tree string is included in the file
  {
    poptreestring_given = 0;
    assert (npops == 2);
    if (npops == 2)
    {
      sscanf (startpoptreestring, "%d", &nloci);
      strcpy (startpoptreestring, "(0,1):2");
    }
    else
    {
      IM_err(IMERR_MISSINGPOPSTRING,"no population string given in file, but # of populations is greater than 2");
    }
  }
  else
  {
    poptreestring_given = 1;
  }
  infilelines += poptreestring_given;
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    SP "Population Tree : Island Model with %d populations\n", npops);
  }
  else
  {
    SP "Population Tree : %s\n", startpoptreestring);
  }

  if (modeloptions[ADDGHOSTPOP])
  {
    add_ghost_to_popstring (startpoptreestring);
    SP "Population Tree with Ghost Population: %s\n", startpoptreestring);
  }
  if (poptreestring_given)
    fscanf (infile, "%d\n", &nloci);
  if (nloci < 1)
  {
    IM_err (IMERR_INFILEFAIL_NLOCI, "The number of loci must be positive: %d is given\n", nloci);
  }
  else if (nloci > MAXLOCI)
  {
    IM_err (IMERR_INFILEFAIL_NLOCI, "The number of loci (%d) is larger than %d\n", nloci, MAXLOCI);
  }
  /* moved this 
  SP "\n\nLocus Information\n");
  SP "-----------------\n");
  SP "\nNumber of loci: %d \n", nloci);
  SP "Locus#\tLocusname");
  for (i = 0; i < npops; i++)
    SP "\tPop%d#", i);
  SP "\tModel\tInheritanceScalar\tMutationRatesPerYear\n");
  */
  infilelines += 3;
  f_close (infile);

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || modeloptions[SINGLEPOPULATION] == 1)
  {
    lastperiodnumber = 1;
    numsplittimes = 0;
  }
  else
  {
    lastperiodnumber = npops - 1;
    numsplittimes = npops - 1;
  }
  SP "\nParameter Priors\n");
  SP "-----------------\n");
  SP "  Population size parameters maximum value : %lf \n", thetaprior);
  if (modeloptions[EXPOMIGRATIONPRIOR])
    SP "  Migration rate parameters exponential distribution mean : %lf \n",
      mprior);
  else
    SP "  Migration rate parameters maximum value: %lf \n", mprior);
  if (modeloptions[SPLITTINGRATEPARAMETER])
  {
    SP "  Splitting rate parameters : %lf \n", splitprior);
  }
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
  {
    SP "  Splitting time : %lf\n", tprior);
  }
  SP "\n");
  return;
}                               //read_datafile_top_lines


void
readdata (char infilename[], char startpoptreestring[], int *fpstri,
          char fpstr[], int **numsitesIS)
{
  /* reads through the data file once at first, to get info on nloci and npops */
  /* then rewinds the file and reads through the data once for each chain */
  /* first line is a comment */
  /* any number of additional comments,  each line begins with '#' */
  /* first data line is number of populations */
  /* second data line is names of each population,  species1 first, followed by species 2 etc */
  /* next line is population tree string,  can be skipped if only 2 populations */
  /* number of loci */
  /* then for each locus: */
  /* name, number of sequences in pop1, number of sequences in pop2, etc etc  length of sequences, mutationrate info 
   * number of sequence unknown
   H for HKY model 
   I for  infinite sites
   S  for stepwise model
   J  for joint infinite sites and stepwise
   inheritance scalar (typically 1 or 0.75 or 0.25) */
  int li, uinext, i;
  char *cc;
  int holdinfilelines;
  char textline[DATAFILEMAXLINELENGTH + 1];
  FILE *infile;

  if ((infile = fopen (infilename, "r")) == NULL)
  {
    IM_err(IMERR_READFILEOPENFAIL,"Error opening text file for reading [filename: %s]\n",infilename);
  }
  skip_datafile_toplines (infile);
  holdinfilelines = infilelines;
  nurates = 0;
  countuprior = 0;
  counturateperyear = 0;

  L = malloc (nloci * sizeof (struct locus));
  uinext = 0;
  SP "\nLocus Information\n");
  SP "-----------------\n");
  SP "\nNumber of loci: %d \n", nloci);
  SP "Locus#\tLocusname");
  for (i = 0; i < npops; i++)
    SP "\tPop%d#", i);
  SP "\tModel\tInheritanceScalar\tMutationRatesPerYear\n");
  for (li = 0; li < nloci; li++)
  {
    fgets (textline, DATAFILEMAXLINELENGTH, infile);
    while ((textline[strlen (textline) - 1] == '\n')
           || (textline[strlen (textline) - 1] == ' '))
      textline[strlen (textline) - 1] = '\0';
    cc = textline;
    /* read in information from infile,  for chains > 0 only read data */
    parse_locus_info (li, &uinext, cc, fpstri, fpstr);
    infilelines++;
    switch (L[li].model)
    {
    case INFINITESITES:
      readseqIS (infile, li, INFINITESITES, numsitesIS);
      break;
    case HKY:
      readseqHKY (infile, li);
      break;
    case STEPWISE:
      readseqSW (infile, li);
      break;
    case STEPWISEP:
      readseqSWP (infile, li);
      break;
    case JOINT_IS_SW:
      readseqIS (infile, li, JOINT_IS_SW, numsitesIS);
      break;
    }

    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      IMA_rearrange_readseq (li);
    }
  }
  if (calcoptions[MUTATIONPRIORRANGE])
  {
    SP " Use prior ranges on mutation rates specified in input file \n");
    if (countuprior <= 1)

    {
      SP " Less than 2 prior ranges given in input file,  mutation rate priors not used \n");
      calcoptions[MUTATIONPRIORRANGE] = 0;
    }
  }
  f_close (infile);

  return;
}                               /* readdata */
