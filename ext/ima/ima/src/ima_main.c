
/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
/*
  IMa 2.0  recent history
  - May 15. 2009 bug fixed re calculating likelihood under HKY model
  - May 4, 2009 bug fixed re  printing burntrend files




*/


#define GLOBVARS

#include "imamp.h"
#include "updateassignment.h"

/* JODY 4/16/09
_UPDATE_GTREE_COVAR removed
*/

extern struct edgemiginfo oldedgemig;
extern struct edgemiginfo oldsismig;
extern struct edgemiginfo newedgemig;
extern struct edgemiginfo newsismig;

#define BURNTRENDSTARTDELAYDEFAULT  10000

/* This is the number of tries of updating assignment, which is equal to number
 * of individuals with unknown origin */
static int snupdatei;              

static FILE *outfile;
char *outfilename;
char *trueassignment;
static char loadfilebase[FNSIZE] = "\0";
static time_t starttime;
static time_t endtime;
static time_t chainstarttime;
static time_t timer;
static time_t lasttime;
static time_t remained_starttime;
static time_t remained_endtime;
static long burnduration, chainduration;
static int burndone;
static long int burnsteps;
static int burntrendstartdelay;
static int burndurationmode, cdurationmode;
static int treestosave;
static int memfortreessaved = 0;
static int savetreeint;
static int recordint;
static int printint;
static double generationtime;
static double scaleumeaninput = 0;
static int swaptries;
static int heatmode;
static char fpstr[50000];       // probably long enough
static int *fpstri;
static char oldoutfilename[FNSIZE];
static double hilocuslike[MAXLOCI];
static int numtreefiles;
static FILE *treeinfosavefile;
static char treeinfosavefilename[FNSIZE];
static long int recordstep = 0;
static double hilike = -1e20, hiprob = -1e20;
static int adaptcheck = 1000;
static FILE *checkdonefile;
static double hval1, hval2;
static int gsampinflength;
static int trenddoublepoint;
static int trendspot = 0;
static int burntrendplotinc;
static int maxedouttreesave = 0;
static char priorfilename[FNSIZE];
static char mcfwritefilename[FNSIZE];
static char mcfreadfilename[FNSIZE];
static char command_line[1001];
static char *infilename;
static long seed_for_ran1;
static int **migcount, *migfrom, *migto;
static char migplotfilename[FNSIZE];
static char migrationnamefilename[FNSIZE];
static int migrationnamefrom,migrationnameto;
static FILE *migrationnamefile;
static FILE *migplotfile;

/*Local function prototypes  */
static void init_after_start_IMA ();
static void init_IMA ();
static void scan_commandline (int argc, char *argv[]);
void print_outputfile_info_string (void);
static void start (int argc, char *argv[]);
static void qupdate (void);
static void savetreeinfo (void);
static void reset_after_burn (void);
static void output_burntrendfile (void);
static int run (void);
static void inctrend (int m, int t, struct value_record *v, double newval);
static void trend_reset (struct value_record *v, int nv);
static void trendrecord (int loadarrayj);
static void recordval (struct value_record *v, double val);
static void record_migrations (void);
static void record (void);
static void loadtreevalues (void);
static void callasciicurves (void);
static void callasciitrend (FILE * outtofile);
static void printoutput (void);
static void intervaloutput (FILE * outto);
static void free_ima_main_stuff ();
static void callprintacceptancerates (FILE * outto);
static void printsteps (FILE * outto, double like);
static void check_to_record (void);
static void record_migration_names();

/* SANGCHUL: Mon Jan 12 20:32:41 EST 2009
 * All global variables can be initialized in function [[init_IMA]].
 * All variables must be finalized in function [[free_IMA]].
 * Sometimes we need to initialze some global variables after we call function
 * [[start]]. Function [[init_after_start_IMA]] does this job.
 * */
void
init_IMA ()
{
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    init_update_assignment ();
  }
  return;
}

void
init_after_start_IMA ()
{
  char Structurama[FNSIZE];

  if (assignmentoptions[PRINTSTRUCTURAMA] == 1)
    {
      strcpy (Structurama, outfilename);
      strcat (Structurama, ".in");
      IMA_convert_IM2Structurama (Structurama);
      IMA_output_structurama_bat (outfilename);
    }

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    recordassignment_header (outfilename);
  }

  /* Function IMA_ninds returns the number of individuals with their label being
   * unknown. */
  /* snupdatei = IMA_nindsunknown (); */
  snupdatei = 1;
  return;
}


void
free_ima_main_stuff ()
{
  int i;
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    free_update_assignment ();
  }
  if (memfortreessaved > 0)
  {
    for (i = 0; i < memfortreessaved; i++)
    {
      XFREE (gsampinf[i]);
    }
    XFREE (gsampinf);
  }
  XFREE (fpstri);
  if (outputoptions[MIGRATEHIST])
  {
    for (i = 0; i < nloci + (nloci > 1); i++)
    {
      XFREE (migcount[i]);
    }
    XFREE (migcount);
    XFREE (migfrom);
    XFREE (migto);
  }

  if (trueassignment != NULL)
  {
    XFREE (trueassignment);
  }

  if (infilename != NULL)
  {
    XFREE (infilename);
  }

  if (outfilename != NULL)
  {
    XFREE (outfilename);
  }

  freeanymemory ();

}                               //free_ima_main_stuff


#define GENERATIONTIMEDEFAULT   1
#define DEFAULTNUMTREES 10000

void
scan_commandline (int argc, char *argv[])
{
  static int i, j, k;
  static char ch, ch1;
  static char pstr[256];
  int Ap, Bp, Cp, Dp, Ep, Fp, Hfp, Hnp, Hkp, Hgp, Hbp, Ip, Jp, Lp, Mp, Op, Pp,
    Qp, Rp, Sp, Tp, Up, Vp, Wp, Yp, Zp;
  int Xp;
  int lenptr;
  double tempf;
  char *opt;
  const char *sep = ",:";
  int ioption;

  time (&starttime);
  time (&remained_starttime);

#ifdef IMPROFILE
  init_clock ();
#endif

  Ap = 0;                       /* Assignment options */
  Bp = 0;                       /* duration of burnin */
  Cp = 0;                       /* calculation options */
  Dp = 0;                       /* number of steps in between tree saves */
  Ep = 0;                       /* prior on splitting rate parameter */
  Fp = 0;                       /* name of mcf file */
  Hfp = 0;                      /* heating model */
  Hnp = 0;                      /* # of chains */
  Hkp = 0;                      /* # of swap attempts */
  Hgp = 0;                      /* heat term1 */
  Hbp = 0;                      /* heat term2 */
  Ip = 0;                       /*input file */
  Jp = 0;                       /* used for programmer options */
  Lp = 0;                       /* duration of chain */
  Mp = 0;                       /* migration rate max */
  Op = 0;                       /* output file */
  Pp = 0;                       /* output options */
  Qp = 0;                       /* Theta max scalar */
  Rp = 0;                       /* run options */
  Sp = 0;                       /* random number seed */
  Tp = 0;                       /* Time maximum */
  Up = 0;                       /* generation time in years */
  Vp = 0;                       /* tree load file name base */
  Wp = 0;                       /* name of file with priors */
  Yp = 0;                       /* mutation rate scalar for loci with mutation rates given in input file - for use with LOADRUN mode  */
  Zp = 0;                       /* screen printout frequency */

  Xp = 0;                       /* True Assignment */
  trueassignment = NULL;

  printf ("executing program ...\n");
  strcpy (command_line, "");
  if ((argc == 2 && (char) toupper (argv[1][1]) == 'H') || argc == 1)
  {
    printf ("IMa 2.0 Program - copyright 2009 by Jody Hey, Rasmus Nielsen and Sang Chul Choi\n");
    printf ("Release date:  April 28, 2009\n");
    printf ("This program is run through a command line interface\n");
    printf ("The program file and the data file should be in the same directory \n");
    printf ("To execute the program, type the program name followed by the necessary command line flags and instructions \n");
    printf ("Command line usage: - upper or lower case letters can be used\n");
    /*  assignment options not included in output for now 4/27/09  jhey */
    printf ("-a  Population assignment options: \n");
    printf ("    0  Turn on check-points\n");
    printf ("    1  Invoke population assignment\n");
    printf ("    2  Relabel update\n");
    printf ("    3  Beerli and Felsenstein update\n");
    printf ("    4  Print info. for DNA Barcoding\n");
    printf ("    5  Island model\n");
    printf ("    6  Local assignment of genes\n");
    printf ("    7  Print input file for Structurama\n");
    printf ("    e.g, -a12 relabel update of assignment with population tree model\n");
    printf ("         -a13 Beerli and Felsenstein update of assignment with population tree model\n");
    printf ("         -a135 Beerli and Felsenstein update of assignment with island model\n");
    printf ("         -a012: turn on checking genealogy integrity additional to option 1,2\n");
    printf ("         -a124: print assignment proportion of a single unknown gene\n");
    printf ("         -a125: relabel update of assignment with island model\n");
    printf ("         -a127: print out the STRUCTURAMA input file additional to option 1,2\n");
    printf ("         -a126: local assignment of genes are allowed\n");
    printf ("-b  Duration of burn  (MCMC mode only)\n");
    printf ("    - if integer, the number of burnin steps \n");
    printf ("    - if floating point, the time in hours between writing of burntrend file\n");
    printf ("         run continues until file " "IMburn" " is no longer present\n");
    printf ("         in the directory, or if present, does not begin with 'y'\n");
    printf ("-c  Calculation options: \n");
    printf ("    0 Likelihood of data functions return a constant - posterior should equal prior \n");
    printf ("    1 Include ranges on mutation rates as priors on mutation rate scalars\n");
    /* do not include this in screen ouput for now,  4/27/09 jh
    printf ("    2 Use diploid information in likelihoods (Levene values)\n"); 
    */
    printf ("-d  number of steps between tree saving (MCMC mode only) (default 100)\n");
    printf ("-e  maximum value of splitting rate parameter (invokes -j0,  -t value is ignored) \n");
    printf ("-f  Name of file with saved Markov chain state generated in previous run - use with -r3\n");
    //printf("-jh  jh personal options \n");
    //printf("      0  alt data format for SW data -  one data line for each allele in each pop,  1st # is allele length, 2nd is # copies \n");
    printf ("-h  Heating terms (MCMC mode only): \n");
    printf ("  -hf heating model: l linear (default); t twostep; a adaptive twostep; g geometric\n");
    printf ("  -hn number of chains \n");
    printf ("  -hk number of chain swap attempts per step \n");
    printf ("  -ha first heating parameter, effect depends on heating model (default 0.05) \n");
    printf ("  -hb second heating parameter, effect depends on heating model  \n");
    printf ("-i  Input file name (no spaces) \n");
    printf ("-j  Model options: \n");
    printf ("    0  Include population splitting rate parameter (specify prior with -e ) \n");
    printf ("    1  Migration only between sister populations (no migration between non-sister populations)\n");
    printf ("    2  One migration for each pair of populations \n");
    printf ("    3  Migration only between sampled populations (ancestral populations have zero migration)\n");
    printf ("    4  add a non-sampled ghost population to the model \n");
    printf ("    5  Separate population size and migration parameters in each period (lots of parameters) \n");
    printf ("    6  no migration in the model\n");
    printf ("    7  migration prior follows exponential distribution with mean given by -m \n");
    printf ("    8  single population\n");
    printf ("-l  Run duration (default: %d genealogies sampled per locus):\n", DEFAULTNUMTREES);
    printf ("     If in MCMC mode (i.e. not loading trees from a previous run) \n");
    printf ("       - if integer, the number of genealogies to save\n");
    printf ("           this value times -d value sets the # of steps in chain after burnin) \n");
    printf ("	      - if floating point, the time in hours between outputs. \n");
    printf ("         run continues until file " "IMrun" " is no longer present\n");
    printf ("         in the directory, or if present, does not begin with 'y'\n");
    printf ("     If in load-tree mode (i.e. using -r0 to load trees from previous run)\n");
    printf ("       - integer indicates number of trees to load from file(s) named with -r\n");
    printf ("-m  migration prior value (maximum for uniform,  mean if exponential distribution is used \n");
    printf ("-o  Output file name (no spaces) default is 'outfile.txt' \n");
    printf ("-p  Output options: \n");
    printf ("    0 turn off trend plots in outfile (default is to print trend plots)\n");
    printf ("    1 turn off plots of marginal curves in outfile (default is to print marginal density plots)\n");
    printf ("    2 print TMRCA histogram for each genealogy (MCMC mode only)\n");
    printf ("    3 print histogram of parameters on demographic scales  (requires mutation rate(s) in data file)\n");
    printf ("    4 print histogram of splitting times divided by prior (do not use with -j0 or when only 2 sampled populations\n");
    printf ("    5 print histograms of population migration rate (2Nm) estimates\n");
    printf ("    6 print pairwise probabilities that one parameter is greater than another \n");
    printf ("    7 print plots of the number of migration events and the mean time of migration events (MCMC mode only)\n");
    printf ("    8 print joint estimate for splitting times (MCMC mode only, for models with 3 or 4 populations)\n");
    printf ("-q  Maximum for population size parameters (4Nu) \n");
    printf ("-r  Run options \n");
    printf ("    0 LOAD-TREE Mode - load trees from previous run(s); also requires -v \n");
    printf ("    1 do not save genealogies to a file (default saves sampled genealogies) \n");
    printf ("    2 save the state of the Markov chain in a file - named with extension mcf (MCMC mode only)\n");
    printf ("    3 start run by loading a previously saved *.mcf file; requires -f (data and priors must be the same) \n");
    printf ("    4 write all mutation related updates rates to stdout during the run (default is to suppress this)\n");
    printf ("    5 print burntrend file at end of burnin period; use with -b followed by integer (MCMC mode only)\n");
    printf ("-s  Random number seed (default is taken from current time)\n");
    printf ("-t  Maximum time of population splitting ( do not use with -e) \n");
    printf ("-u  Generation time in years - for use with -p3 (default is %d) \n", GENERATIONTIMEDEFAULT);
    //    printf ("-w  Name of file with parameter priors  (requires -c3)  default: 'imapriors.txt'\n");
    printf ("-v  Base name (no extension) of *.ti files with tree data  (requires use of -r0) \n");
    printf ("-y  Mutation rate scalar for relevant loci - for use with -p3 \n");
    printf ("-x  beta for rasing the power to likelihood\n");
    printf ("-z  Number of steps between screen output (default 10000) (MCMC mode only)\n");
    exit (0);
  }
  else
  {
/*
command line circumstances:
all flags begin with '-'
-most flags are single letter flags
-some are double letter flags: h
-some flags are followed by a string or a character
others by a single number (int or float)
others by a string of integers

it is ok to have spaces between a flag and its values 

All flags are followed by at least something
no flag is followed by nothing 
*/
    for (i = 1; i < argc; i++)
    {
      strcpy (pstr, argv[i]);
      strcat (command_line, " ");
      strcat (command_line, pstr);

      if (strlen (pstr) < 2)
        IM_err (IMERR_COMMANDLINEFORMAT, " one of the command line strings is too short: %s ",pstr);
      if (pstr[0] != '-')
        IM_err (IMERR_COMMANDLINEFORMAT, "command line flag not preceded by '-' : %s", pstr);
      ch = toupper (pstr[1]);
      if (ch == 'H')
        ch1 = toupper (pstr[2]);
      else
        ch1 = ' ';
      if (strlen (argv[i]) == 2 || (i < argc - 1 && isdigit (argv[i + 1][0])))  // space separates flag from its number
      {
        i++;
        strcpy (pstr, argv[i]);
        strcat (command_line, " ");
        strcat (command_line, pstr);
      }
      else
      {
        if ((ch == 'H'))
          strdelete (pstr, 1, 3);
        else
          strdelete (pstr, 1, 2);
      }
      switch ((char) toupper (ch))
      {
      case 'A':
        j = (int) (strlen (pstr) - 1);
        if (comma_exists (pstr))
        {
          for (opt = strtok (pstr, sep); opt; opt = strtok (NULL, sep))
          {
            ioption = atoi (opt);
            if (ioption < 0 || ioption >= POPULATIONASSIGNMENT_NUMBER)
            {
              IM_err (IMERR_COMMANDLINEFORMAT, "option -a %s", pstr);
            }
            assignmentoptions[ioption] = 1;
          }
          Ap = 1;
        }
        else
        {
          while (j >= 0)
          {
            ioption = atoi (&pstr[j]);
            if (ioption < 0 || ioption >= POPULATIONASSIGNMENT_NUMBER)
            {
              IM_err (IMERR_COMMANDLINEFORMAT, "option -a %s", pstr);
            }
            assignmentoptions[ioption] = 1;
            pstr[j] = '\0';
            j--;
            Ap = 1;
          }
        }
        break;
      case 'B':
        tempf = atof (&pstr[0]);
        /* check to see if the value is floating point, in which case treat it as being in fractions of an hour  and convert to seconds */
        if (strchr (pstr, '.'))
        {
          burnduration = (int) (3600 * tempf);
          burndurationmode = TIMEINF;
          time (&lasttime);
          runoptions[PRINTBURNTREND] = 1;
        }
        else
        {
          burnduration = (int) tempf;
          burndurationmode = TIMESTEPS;
        }
        Bp = 1;
        break;
      case 'C':
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
            IM_err (IMERR_COMMANDLINEFORMAT, "calculation option flag -c should be followed by a digit: %s",pstr);
          calcoptions[atoi (&pstr[j])] = 1;
          pstr[j] = '\0';
          j--;
          Cp = 1;
        }
        break;
      case 'D':
        savetreeint = atoi (&pstr[0]);
        Dp = 1;
        break;
      case 'E':
        splitprior = (double) atof (&pstr[0]);
        Ep = 1;
        break;
      case 'F':
        strcpy (mcfreadfilename, pstr);
        Fp = 1;
        break;
      case 'H':
        switch ((char) toupper (ch1))
        {
        case 'A':
          hval1 = atof (&pstr[0]);
          Hgp = 1;
          break;
        case 'B':
          hval2 = atof (&pstr[0]);
          Hbp = 1;
          break;
        case 'N':
          numchains = atoi (&pstr[0]);
          Hnp = 1;
          break;
        case 'K':
          swaptries = atoi (&pstr[0]);
          Hkp = 1;
          break;
        case 'F':
          Hfp = 1;
          switch ((char) toupper (pstr[0]))
          {
          case 'T':
            heatmode = HTWOSTEP;
            break;
          case 'A':
            heatmode = HADAPT;
            break;
          case 'G':
            heatmode = HGEOMETRIC;
            break;
          default:
            heatmode = HLINEAR;
            break;
          }
          break;
        default:
          IM_err (IMERR_COMMANDLINEFORMAT, "mistake in use of -h flag : %s", pstr);
        }
        break;
      case 'I':
        lenptr = strlen (pstr);
        infilename = (char *) malloc ((lenptr + 1) * sizeof (char));
        strcpy (infilename, pstr);
        Ip = 1;
        break;
      case 'J':
        Jp = 0;
        if (!(toupper(pstr[0]) == 'H'))
        {
          j = (int) (strlen (pstr) - 1);
          while (j >= 0)
          {
            if (!isdigit (pstr[j]))
              IM_err (IMERR_COMMANDLINEFORMAT, "model option flag -j should be followed by a digit: %s", pstr);
            modeloptions[atoi (&pstr[j])] = 1;
            pstr[j] = '\0';
            j--;
          }
        }
        else
        {
          j = (int) (strlen (pstr) - 1);
          if (isdigit(pstr[1])  && pstr[1] == '1')
          {
            migrationnamefrom = atoi (&pstr[2]);
            migrationnameto = atoi (&pstr[3]);
            strcpy(migrationnamefilename, &pstr[4]);
            jheyoptions[WRITEMIGRATIONNAME] = 1;
          }
          else
          {
            while (j >= 0)
            {
              if (isdigit(pstr[j]))
              {
                jheyoptions[atoi (&pstr[j])] = 1;
              }
              pstr[j] = '\0';
              j--;
            }
          }
        };
        break;
      case 'L':
        tempf = atof (&pstr[0]);
        /* check to see if the value is floating point, in which case treat it as being in fractions of an hour  and convert to seconds */
        if (strchr (pstr, '.'))
        {
          chainduration = (int) (3600 * tempf);
          cdurationmode = TIMEINF;
          treestosave = -1;
        }
        else
        {
          treestosave = (int) tempf;
          cdurationmode = TIMESTEPS;
        }
        Lp = 1;
        break;
      case 'M':
        if (modeloptions[SINGLEPOPULATION] == 1)
        {
          IM_err (IMERR_COMMANDLINEFORMAT, "model option [single population] should not go with -m");
        }
        mprior = (double) atof (&pstr[0]);
        if (mprior == 0)        // turn the use of this off  9/15/08
          modeloptions[NOMIGRATION] = 1;
        /*if (mprior == 0)
           mprior = MPRIORMIN;  not needed ?? */
        Mp = 1;
        break;
      case 'O':
        lenptr = strlen (pstr);
        outfilename = (char *) malloc ((lenptr + 1) * sizeof (char));
        strcpy (outfilename, pstr);
        Op = 1;
        break;
      case 'P':
        Pp = 1;
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
          {
            IM_err (IMERR_COMMANDLINEFORMAT, "print option flag -p should be followed by a digit : %s",pstr);
          }
          k = atoi (&pstr[j]);
          outputoptions[k] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
      case 'Q':
        thetaprior = (double) atof (&pstr[0]);
        Qp = 1;
        break;
      case 'R':
        Rp = 1;
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
            IM_err (IMERR_COMMANDLINEFORMAT, "run option flag -r should be followed by a digit : %s ", pstr);
          k = atoi (&pstr[j]);
          runoptions[k] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;

      case 'S':
        seed_for_ran1 = atoi (&pstr[0]);
        if (!seed_for_ran1)
          seed_for_ran1 = 1;
        Sp = 1;
        break;
      case 'T':
        if (modeloptions[SINGLEPOPULATION] == 1)
        {
          IM_err (IMERR_COMMANDLINEFORMAT, "model option [single population] should not go with -t");
        }
        tprior = (double) atof (&pstr[0]);
        Tp = 1;
        break;
      case 'U':
        generationtime = atof (&pstr[0]);
        Up = 1;
        break;
      case 'V':
        strcpy (loadfilebase, pstr);
        Vp = 1;
        break;
      case 'W':
        strcpy (priorfilename, pstr);
        Wp = 1;
        break;
      case 'Y':
        scaleumeaninput = atof (&pstr[0]);
        Yp = 1;
        break;
      case 'Z':
        printint = atoi (&pstr[0]);
        Zp = 1;
        break;
      case 'X':
        if (strchr (pstr, '.') == NULL)
        {
          lenptr = strlen (pstr);
          trueassignment = (char *) malloc ((lenptr + 1) * sizeof (char));
          strcpy (trueassignment, pstr);
          gbeta = 1.0;
        }
        else
        {
          gbeta = atof (&pstr[0]);
          trueassignment = NULL;
        }
        Xp = 1;
        break;
      default:
        IM_err (IMERR_COMMANDLINEFORMAT, &ch);
      }
    }
  }

  if (Xp == 0)
  {
    gbeta = 1.0;
  }

  if (!Ip)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " No data file given on command line");
  }
  if (!Lp && !runoptions[LOADRUN])
  {
    treestosave = (int) DEFAULTNUMTREES;
    cdurationmode = TIMESTEPS;
  }
  if (runoptions[LOADRUN] || modeloptions[NOMIGRATION])
    outputoptions[MIGRATEHIST] = 0;
  if (!Op)
    strcpy (outfilename, "outfile.txt");
#ifdef _DEBUG
  checkoutfileclosed (&outfile, outfilename);   // just make sure that outfilename does not name a file that is already opened
#endif
  if (runoptions[LOADRUN] && !Vp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " -r0 invoked without -v information, i.e. no base name for files containing trees was given on the command line");
  }
  if (!Hnp || runoptions[LOADRUN])
  {
    numchains = DEFCHAINS;      /* default value */
  }
  else
  {
    if ((heatmode == HGEOMETRIC && numchains < 4)
        || (heatmode == HTWOSTEP && numchains < 3) 
        || (heatmode == HADAPT && numchains < 3))
    {
      IM_err (IMERR_MISSINGCOMMANDINFO, "too few chains specified in heating model");
    }
    if (!Hfp)
    {
      heatmode = HLINEAR;
      if (!Hgp)
        hval1 = 0.05;           /* default value */
    }
    else
    {
      if (heatmode > HLINEAR)
      {
        if (heatmode == HTWOSTEP)
        {
          if (!Hgp)
            hval1 = 0.05;
          if (!Hbp)
            hval2 = 2;
        }
        if (heatmode == HADAPT)
        {
          if (!Hgp)
            hval1 = 0.05;
          if (!Hbp)
            hval2 = 1;
        }
        if (heatmode == HGEOMETRIC)
        {
          if (!Hgp)
            hval1 = 0.9;
          if (!Hbp)
            hval2 = 0.8;
        }
      }
      else if (!Hgp)
        hval1 = 0.05;           /* default value */
    }
  }
  recordint = RECORDINTDEFAULT;
  if (!Up)
    generationtime = GENERATIONTIMEDEFAULT;
  if (!Dp)
    savetreeint = SAVETREEINTDEFAULT;
  if (!Zp)
    printint = PRINTINTDEFAULT;
  if (!Bp && !runoptions[LOADRUN])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            "No burn duration information given on command line (use -b)");
  }
  if (modeloptions[NOMIGRATION] || mprior == 0)
  {
    modeloptions[NOMIGRATION] = 1;
    mprior = 0;
  }
  if (Ep)
  {
    modeloptions[SPLITTINGRATEPARAMETER] = 1;
  }
  if (!Tp && !modeloptions[SPLITTINGRATEPARAMETER]
      && !assignmentoptions[POPULATIONASSIGNMENTINFINITE]
      && !modeloptions[SINGLEPOPULATION])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No prior information provided for population splitting times (-t)");
  }
  if (!Mp && modeloptions[NOMIGRATION] != 1)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value for migration parameter (-m)");
  }
  if (!Ep && modeloptions[SPLITTINGRATEPARAMETER])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No prior information provided for population splitting rate (-e)");
  }
  if (!Qp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value of 4Nu parameters");
  }

  if (runoptions[LOADMCSTATE] && !Fp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            "  -r3 invoked without -f information, i.e. no filename given for markov chain state file,  for loading state of markov chain ");
  }
  if (!Hkp)
    swaptries = 1;
  else
  {
    if (numchains > 1 && swaptries > numchains * (numchains - 1) / 2)
      swaptries = numchains * (numchains - 1) / 2;
  }


  if (runoptions[PRINTBURNTREND])
  {
    if (runoptions[LOADMCSTATE])
      burntrendstartdelay = 0;
    else
      burntrendstartdelay = BURNTRENDSTARTDELAYDEFAULT;
  }
  if (!Sp)
    seed_for_ran1 = (long) time (NULL);
  if (strcmp (infilename, outfilename) == 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT, " Input and output file names are identical");
  }
  if (cdurationmode == TIMESTEPS)
  {
    chainduration = treestosave * savetreeint;
  }
  if (jheyoptions[WRITEMIGRATIONNAME])
  {
    migrationnamefile = fopen (migrationnamefilename, "w");
  }

  /* IMa stops its running if there is any conflicted options. */
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    if (modeloptions[NOMIGRATION] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model must be allowed for migration events: do not use -j%d with -a%d", 
              NOMIGRATION, 
              POPULATIONASSIGNMENTINFINITE);
    }
    if (modeloptions[SPLITTINGRATEPARAMETER] == 1
        || modeloptions[NOMIGBETWEENNONSISTERS] == 1
        || modeloptions[PARAMETERSBYPERIOD] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has only a single period: do not use -j%d, -j%d, or -j%d with -a%d", 
              SPLITTINGRATEPARAMETER, NOMIGBETWEENNONSISTERS, PARAMETERSBYPERIOD, 
              POPULATIONASSIGNMENTINFINITE);
    }
    if (outputoptions[THISTDIVIDEBYPRIOR] == 1
        || outputoptions[PRINTJOINTTEST] == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has no split time: do not use -p%d or -p%d with -a%d", 
              THISTDIVIDEBYPRIOR, 
              PRINTJOINTTEST,
              POPULATIONASSIGNMENTINFINITE);
    }
    if (Tp == 1)
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, 
              "Island model has no split time: do not use -t with -a%d", 
              POPULATIONASSIGNMENTINFINITE);
    }
  }

  return;
}                               // scan_commandline 


/*  this prints basic info to a string, fpstr, that later gets printed to the output file */
void
print_outputfile_info_string (void)
{
  fpstri = (int *) malloc (sizeof (int));
  *fpstri = 0;
  SP "IMa version 2.0 - Isolation with Migration Analysis  -  Jody Hey, Rasmus Nielsen, Sang Chul Choi 2009 \n\n");
  SP "\nINPUT AND STARTING INFORMATION \n");
  SP "================================\n");
  SP "\nCommand line string : %s \n", command_line);
  SP "Input filename : %s \n", infilename);
  SP "Output filename: %s \n", outfilename);
  SP "Random number seed : %li \n", seed_for_ran1);
  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
    SP "**NO DATA ** - Data likelihoods set to constant  posterior should equal prior \n");
  if (calcoptions[CALCLEVINELIKELIHOOD])
    {
      SP "**Diploid Configurations used for calculating P(D|G)\n");
    }
  if (!runoptions[LOADRUN])
  {
    SP "- Run Duration - \n");
    switch (burndurationmode)
    {
    case TIMESTEPS:
      SP "     Burn period, # steps: %li \n", burnduration);
      break;
    case TIMEINF:
      SP "     Burn period, # seconds: %li (total burn duration depends on IMburn file status)\n", burnduration);
      break;
    };
    if (runoptions[PRINTBURNTREND])
    {
      SP "          -User option for printing trendline during, or at end of burnin period, invoked\n");
      SP "          -initial burn duration prior to beginning recording burntrend : %d steps\n", burntrendstartdelay);
    }

    switch (cdurationmode)
    {
    case TIMESTEPS:
      SP "     Record period, #saves: %d  #steps each: %li   total #steps: %li \n", treestosave, chainduration / treestosave, chainduration);
      break;
    case TIMEINF:
      SP "      Record period, # seconds per interval: %ld \n",
        chainduration);
      SP "      Run period indefinite in length - determined by user using 'IMrun' file\n");
      break;
    };
    SP "- Metropolis Coupling -\n");
    if (numchains > 1)
    {
      SP "     Metropolis Coupling implemented using %d chains \n", numchains);
      switch (heatmode)
      {
      case HLINEAR:
        SP "     Linear Increment Model   term: %.3f\n", hval1);
        break;
      case HTWOSTEP:
        SP "     Twostep Increment Model   term1: %.3f  term2: %.3f\n",
          hval1, hval2);
        break;
      case HADAPT:
        SP "     Twostep Adaptive Model  starting term1: %.3f  starting term2: %.3f\n", hval1, hval2);
        break;
      case HGEOMETRIC:
        SP "     Geometric Increment Model   term1: %.3f  term2: %.3f\n",
          hval1, hval2);
        break;
      }
    }
    else
      SP "     None \n");
  }
  if (runoptions[LOADMCSTATE])
  {
    SP "\n- Initial Markov chain state space loaded from file: %s\n",
      mcfreadfilename);
  }
  if (runoptions[SAVEMCSTATEFILE])
  {
    SP "State of Markov chain saved to file : %s\n", mcfwritefilename);
  }
  if (!runoptions[DONTSAVETREES])
  {
    SP "All tree information used for surface estimation printed to file: %s\n", treeinfosavefilename);
  }
  if (outputoptions[PRINTTMRCA])
    SP "TMRCA  histograms printed \n");
}                               //  print_outputfile_info_string


// sets up filenames and calls the main initialization functions 
void start (int argc, char *argv[])
{
  int i;

  scan_commandline (argc, argv);
  if (runoptions[SAVEMCSTATEFILE])
  {
    strcpy (mcfwritefilename, outfilename);
    strcat (mcfwritefilename, ".mcf");
  }
  if (!runoptions[DONTSAVETREES])
  {
    strcpy (treeinfosavefilename, outfilename);
    strcat (treeinfosavefilename, ".ti");
  }
  if (cdurationmode == TIMEINF)
  {
    strcpy (oldoutfilename, outfilename);
    strcat (oldoutfilename, ".old");
  }
  for (i = 0; i < MAXLOCI; i++)
    hilocuslike[i] = -1e20;
  print_outputfile_info_string ();
  setseeds (seed_for_ran1);
  setlogfact ();
  setup (infilename, fpstri, fpstr);
  setheat (hval1, hval2, heatmode);
  checkautoc (1, 0, 0);
  gsampinflength = calc_gsampinf_length ();
  if (outputoptions[MIGRATEHIST])
  {
    strcpy (migplotfilename, outfilename);
    strcat (migplotfilename, ".mpt");
    migcount = (int **) malloc ((nloci + (nloci > 1)) * sizeof (int *));
    for (i = 0; i < nloci + (nloci > 1); i++)
    {
      migcount[i] = (int *) malloc (nummigrateparams * sizeof (int));
    }
    migfrom = (int *) malloc (nummigrateparams * sizeof (int));
    migto = (int *) malloc (nummigrateparams * sizeof (int));
    for (i = 0; i < nummigrateparams; i++)
    {
      if (strchr (imig[i].str, ',') == NULL)    // migration parameters note done by period
      {
        if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
        {
          migfrom[i] = atoi (&imig[i].str[1]);
          migto[i] = atoi (&imig[i].str[4]);

        }
        else
        {
          migfrom[i] = atoi (&imig[i].str[1]);
          migto[i] = atoi (&imig[i].str[3]);
        }
      }
      else
      {
        if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
        {
          migfrom[i] = atoi (&imig[i].str[3]);
          migto[i] = atoi (&imig[i].str[6]);

        }
        else
        {
          migfrom[i] = atoi (&imig[i].str[3]);
          migto[i] = atoi (&imig[i].str[5]);
        }
      }
    }
  }
  if (runoptions[LOADMCSTATE])
  {
    readmcf (mcfreadfilename);
  }

}                               /* start */



#define TUPDATEINC  0           //1            // do a single t parameter every TUPDATEINC steps
#define UUPDATEINC  9           // u parameters seem to mix well so do every 10 steps


#define UPDATEGENEALOGYCOVARFRAC  0.2   // proportion of time that updategenealogy_covar() is called rather than regular updategenealogy()



void qupdate (void)
{
  int i;
  int j, k, li, ci, ui;
  int changed;
  int qswapped = 0;
  int topolchange, tmrcachange;
  int periodpick, tupdatemethodpick;
  static int tui = 0;
  //static int nexttinc = 0;
  static int uui = 0;

  /* update genealogies */
  for (ci = 0; ci < numchains; ci++)
  {
    for (li = 0; li < nloci; li++)
    {
      if (uniform () > UPDATEGENEALOGYCOVARFRAC)
      {
        if (ci == 0)
        {
          L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries++;
          L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries++;
          L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries++;
        }

#ifdef IMPROFILE
        startclock (&clock_updategenealogy);
#endif //IMPROFILE
        if (updategenealogy (ci, li, &topolchange, &tmrcachange))
        {
          if (ci == 0)
          {
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp++;
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp += (topolchange > 0);
            L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp += (tmrcachange > 0);
          }
        }

#ifdef IMPROFILE
        addclock (&clock_updategenealogy);
#endif //IMPROFILE
      }
      else
      {
        /* SANGCHUL: Multiple-Branch update does not see eye to eye with
         * assignment update at the moment. */
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1
            || assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
        {
          continue;
        }
        if (ci == 0)
          L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_COVAR].tries++;
  
#ifdef IMPROFILE
        startclock (&clock_updategenealogy_covar);
#endif //IMPROFILE
        if (updategenealogy_covar (ci, li) && ci == 0)
          L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_COVAR].accp++;
#ifdef IMPROFILE
        addclock (&clock_updategenealogy_covar);
#endif //IMPROFILE
      }
    }
  }

  /* update population assignment */
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    if (assignmentoptions[POPULATIONASSIGNMENTLOCAL] == 0)
    {
      for (ci = 0; ci < numchains; ci++)
      {
        for (i = 0; i < snupdatei; i++)
        {
          if (assignmentoptions[POPULATIONASSIGNMENTRELABEL] == 1)
          {
            updateassignmentrelabel (ci);
          }
          if (assignmentoptions[POPULATIONASSIGNMENTBF] == 1)
          {
            updateassignmentbf (ci);
          }
        }
      }
    }
    else
    {
      for (ci = 0; ci < numchains; ci++)
      {
        for (li = 0; li < nloci; li++)
        {
          if (assignmentoptions[POPULATIONASSIGNMENTRELABEL] == 1)
          {
            updateassignmentrelabellocus (ci, li);
          }
        }
      }
    }
    /* for DNA Barcoding */
    if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
    {
      IMA_assigned_increase ();
    }
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0
      && modeloptions[SINGLEPOPULATION] == 0
      && tui == TUPDATEINC)
  {
    for (ci = 0; ci < numchains; ci++)
    {
      if (numsplittimes > 1)
      {
        assert (C[0]->tvals[1] < T[1].pr.max);
      }
      periodpick = randposint (numsplittimes);
      //  tupdatemethodpick = (periodpick < numsplittimes - 1) ? randposint (3) : randposint (2);
      tupdatemethodpick = randposint (2);       // no ry2 update
      if (modeloptions[NOMIGRATION])
        tupdatemethodpick = 0;  // must do ry1 , NW does not work with no migration

      switch (tupdatemethodpick)
      {

      case 0:
#ifdef DO_RY1UPDATE
#ifdef IMPROFILE
        startclock (&clock_changet_RY1);
#endif
        changed = changet_RY1 (ci, periodpick);
#ifdef IMPROFILE
        addclock (&clock_changet_RY1);
#endif
        if (ci == 0)
        {
          T[periodpick].upinf[IM_UPDATE_TIME_RY1].tries++;
          if (changed)
            T[periodpick].upinf[IM_UPDATE_TIME_RY1].accp++;
        }
#endif //DO_RY1UPDATE
        break;

      case 1:
#ifdef DO_NWUPDATE
#ifdef IMPROFILE
        startclock (&clock_changet_NW);
#endif
        /* SANGCHUL: NW update does not see eye to eye with
         * assignment update at the moment. */
        if (assignmentoptions[POPULATIONASSIGNMENTBF] == 1)
        {
          changed = changet_RY1 (ci, periodpick);
        }
        else
        {
          changed = changet_NW (ci, periodpick);
        }
#ifdef IMPROFILE
        addclock (&clock_changet_NW);
#endif
        if (ci == 0)
        {
          T[periodpick].upinf[IM_UPDATE_TIME_NW].tries++;
          if (changed)
            T[periodpick].upinf[IM_UPDATE_TIME_NW].accp++;
        }
#endif /* DO_NWUPDATE */
        break;

      case 2:
#ifdef DO_RY2UPDATE
#ifdef IMPROFILE
        startclock (&clock_changet_RY2);
#endif
        changed = changet_RY2 (ci, periodpick);
#ifdef IMPROFILE
        addclock (&clock_changet_RY2);
#endif

        if (ci == 0)
        {
          T[periodpick].upinf[IM_UPDATE_TIME_RY2].tries++;
          if (changed)
            T[periodpick].upinf[IM_UPDATE_TIME_RY2].accp++;
        }
#endif //DO_RY2UPDATE
        break;

      }
    }
    tui = 0;
  }
  else
  {
    tui++;
  }

  if (uui == UUPDATEINC)
  {
    if (nurates > 1)
      for (ci = 0; ci < numchains; ci++)
      {
        for (j = 0; j < (nurates - (nurates == 2)); j++)
        {
#ifdef IMPROFILE
          startclock (&clock_changeu);
#endif
          ui = changeu (ci, j, &k);
#ifdef IMPROFILE
          addclock (&clock_changeu);
#endif
          if (ci == 0)

          {
            L[ul[j].l].u_rec[ul[j].u].upinf->tries++;
            L[ul[k].l].u_rec[ul[k].u].upinf->tries++;
            if (ui == 1)

            {
              L[ul[j].l].u_rec[ul[j].u].upinf->accp++;
              L[ul[k].l].u_rec[ul[k].u].upinf->accp++;
            }
            /*if (ui == -1)

               {
               C[ci]->G[ul[j].l].u[ul[j].u].uprior_failedupdate++;
               C[ci]->G[ul[k].l].u[ul[k].u].uprior_failedupdate++;
               } */
          }
        }
      }
    else
    {
      if (nloci == 1 && L[0].model == HKY)      /* if there is just one HKY locus kappa needs updating on its own */
        for (ci = 0; ci < numchains; ci++)
          changekappa (ci);
    }
    uui = 0;
  }
  else
  {
    uui++;
  }

  if (numchains > 1)
  {
#ifdef IMPROFILE
    startclock (&clock_swapchains);
#endif
    qswapped =
      swapchains (swaptries, burndone, heatmode, adaptcheck, hval1, hval2);
#ifdef IMPROFILE
    addclock (&clock_swapchains);
#endif
  }

  intervaloutput (stdout);
  if (step >= CHECKAUTOCWAIT)
    checkautoc (0, burndone, burnsteps);

  return;
}                               /* qupdate */



void reset_after_burn (void)
{
  int ci;
  int li, ui, i, j;
  time (&chainstarttime);
  burnsteps = step - 1;
  if (cdurationmode == TIMEINF)
  {
    time (&lasttime);
    time (&timer);
  }
  // reset acceptance rate accumulators

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || modeloptions[SINGLEPOPULATION] == 1)
  {
    /* no split times */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
      for (j = 0; j < T[i].num_uptypes; j++)
        T[i].upinf[j].accp = T[i].upinf[j].tries = 0;
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (ci = 0; ci < numchains; ci++)
    {
      for (j = 0; j < Cupinf[ci].num_uptypes; j++)
      {
        Cupinf[ci].upinf[j].accp = 0;
        Cupinf[ci].upinf[j].tries = 0;
      }
    }
  }

  for (li = 0; li < nloci; li++)
  {
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      for (j = 0; j < L[li].a_rec->num_uptypes; j++)
      {
        L[li].a_rec->upinf[j].accp = 0;
        L[li].a_rec->upinf[j].tries = 0;
      }
      if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
      {
        IMA_assigned_reset ();
      }
    }
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
      for (j = 0; j < L[li].u_rec[ui].num_uptypes; j++)
        L[li].u_rec[ui].upinf[j].accp = L[li].u_rec[ui].upinf[j].tries = 0;
      if (L[li].model == HKY)
        for (j = 0; j < L[li].kappa_rec->num_uptypes; j++)
          L[li].kappa_rec->upinf[j].accp = L[li].kappa_rec->upinf[j].tries =
            0;
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      {
        for (j = 0; j < L[li].A_rec[ui].num_uptypes; j++)
          L[li].A_rec[ui].upinf[j].accp = L[li].A_rec[ui].upinf[j].tries = 0;
      }
      for (j = 0; j < L[li].g_rec->num_uptypes; j++)
        L[li].g_rec->upinf[j].accp = L[li].g_rec->upinf[j].tries = 0;
    }
  }
/*updated the way checkautoc() initializes on 3/25/08 */
  checkautoc (1, burndone, burnsteps);
  adaptcheck = 10000;
  if (!runoptions[DONTSAVETREES])
  {
    if ((treeinfosavefile = fopen (treeinfosavefilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL, "Error creating file for holding tree information");
    }
    fprintf (treeinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (treeinfosavefile, "Header for tree file:  %s\n\n",
             treeinfosavefilename);
    fprintf (treeinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (treeinfosavefile, "%s\n", fpstr);
    fprintf (treeinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (treeinfosavefile, "End of header for tree file:  %s\n\n",
             treeinfosavefilename);
    fprintf (treeinfosavefile,
             "-------------------------------------------\n\n");
    fprintf (treeinfosavefile, "VALUESSTART\n");
    f_close (treeinfosavefile);
  }
  if (runoptions[SAVEMCSTATEFILE])
  {
    writemcf (mcfwritefilename);
  }
  gloglikelihood = C[0]->allpcalc.pdg; 
}                               /* reset_after_burn() */

void output_burntrendfile (void)
{
  FILE *burntrendfile;
  char burntrendfilename[FNSIZE];
  assert (runoptions[PRINTBURNTREND]);
  strcpy (burntrendfilename, outfilename);
  strcat (burntrendfilename, ".burntrend.out");
  if ((burntrendfile = fopen (burntrendfilename, "w")) == NULL)
  {
    IM_err (IMERR_CREATEFILEFAIL,
            "Error opening burntrend output file for writing");
  }
  printf ("\n\n========================\n");
  printf ("Printing Burn Trend File\n");
  printf ("========================\n\n");
  fprintf (burntrendfile,
           "Plots of Runtime Information and Parameter Trends during Burnin \n");
  fprintf (burntrendfile, "========================================\n");
  intervaloutput (burntrendfile);
  fprintf (burntrendfile, "========================================\n\n");
  fprintf (burntrendfile, "Current Step #: %d \n\n", step);
  if (trendspot > 1)
    callasciitrend (burntrendfile);
  else
    fprintf(burntrendfile, "burn period too short to plot trends,  trend recording begins at step %d \n",burntrendstartdelay); 
  fclose (burntrendfile);
}                               // output_burntrendfile

#define CHECKINTERVALSTEPS 1000

/* run() determines the status of a run in terms of whether it is in the burnin phase or not
  and of whether the run should keep going. 
  
  if the burnin period has just ended,  some work is done

  run()  also  checks to see if it is time to print an output file */

int run (void)
{
  static int checkinterval = 0;
  static int printburntrendstep = 0, burnrecordi = 1;
  int tempburndone;
  char ch;
  if (burndone)
  {
    switch (cdurationmode)
    {
    case TIMESTEPS:
      return (step < (chainduration + burnsteps));
      break;
    case TIMEINF:
      if (checkinterval < CHECKINTERVALSTEPS)
      {
        checkinterval++;
        return (1);
      }
      else
      {
        checkinterval = 0;
        if (maxedouttreesave)
          return (0);
        time (&timer);
        if ((timer - lasttime) > chainduration)
        {
          if ((checkdonefile = fopen ("IMrun", "r")) == NULL)
          {
            return (0);
          }
          else
          {
            ch = (char) getc (checkdonefile);
            f_close (checkdonefile);
            if ((char) toupper (ch) != 'Y')
            {
              return (0);
            }
            else
            {
              printoutput ();
              time (&lasttime); // start the clock again 
            }
          }
        }
        return (1);
      }
      break;
    default:
      return (0);
      break;
    }
  }
  else
  {
    tempburndone = 0;
    switch (burndurationmode)
    {
    case TIMESTEPS:
      tempburndone = (step > burnduration);
      break;
    case TIMEINF:
      if (checkinterval < CHECKINTERVALSTEPS)
      {
        checkinterval++;
       // return (1);
      }
      else
      {
        checkinterval = 0;
        time (&timer);
        tempburndone = (timer - lasttime) > burnduration;
      }
      break;
    default:
      return (0);
      break;
    }
    // plot burn trend if called for 
    if (runoptions[PRINTBURNTREND] && (tempburndone || step >= burntrendstartdelay))
    {
      if (burnrecordi == recordint)
      {
        trendrecord (-1);       // record values of parameters that are in mcmc
        burnrecordi = 1;
      }
      else
      {
        burnrecordi++;
      }
      if (tempburndone)
      {
        output_burntrendfile ();

        printburntrendstep = 1;
        /* now check to see if IMburn file is present */
        if (burndurationmode == TIMEINF)
        {
          if ((checkdonefile = fopen ("IMburn", "r")) != NULL)
          {
            ch = (char) getc (checkdonefile);
            f_close (checkdonefile);
            if ((char) toupper (ch) == 'Y')
            {
              tempburndone = 0;
              time (&lasttime); // start the clock again 
            }
          }
        }
      }
      else
      {
        printburntrendstep++;
      }
    }
    if (tempburndone)
    {
      burndone = 1;
      reset_after_burn ();
    }
    return (1);
  }
}                               /* run */

#define  TRENDLASTPT (TRENDDIM - 1)
void inctrend (int m, int t, struct value_record *v, double newval)
{
  int j;
  for (j = m; j < TRENDLASTPT; j++)
  {
    //assert (v->trend[j + 1] != 0);
    v->trend[j] = v->trend[j + 1];
  }
  v->trend[t] = newval;
  //assert(v != 0);
}

void trend_reset (struct value_record *v, int nv)
{
  int i, j;
  for (i = 0; i < nv; i++)
    for (j = 0; j < TRENDDIM; j++)
      v[i].trend[j] = 0;
}                               // trend_reset


/* Using trendrecord()
This function records trend lines
It works on instances of struct value_record 
the value_record must first be initialized (probably in initialize.c) 
Code for a particular value_record or array of value_records 
can be placed in trendrecord at two places:
  in the "if (burndone && reset == 0) " section
  and in the "if (recordinc == recordtrendinc)" section


explanation for how trendline data are accumulated:
 - movespot is the point at which values are deleted from the array, 
 - each new value is added to the end of the array (i.e. at trendspot)
 - all values from one position past movespot up to the end of the array are moved down one position
 - this erases the value at movespot and makes room for a new value at the end. 
 -each time the replacement point (movespot) reaches the end of the array the time 
	period between additions doubles 
values are added more slowly as the run proceeds.  
- this is because the time period doubles when movespot reaches the end 
- the values to the left of movespot have half the density (in time) of those to the right 
*/
void trendrecord (int loadarrayj)
{
  static int /*trendspot = 0,*/ recordtrendinc = 1, recordinc = 0, movespot =
    TRENDLASTPT;
  static int init = 1, reset = 0;
  int j, k, li, ui;

  if (burndone && reset == 0)   // reset all trend-related values after burnin
  {
    init = 1;
    reset = 1;
    if (lpgpd_v->do_trend)
    {
      trend_reset (lpgpd_v, 1);
    }
    for (li = 0; li < nloci; li++)
    {
      if (runoptions[LOADRUN] == 0)
      {
        if (L[li].g_rec->v->do_trend)
        {
          trend_reset (L[li].g_rec->v, 1);
        }
        for (ui = 0; ui < L[li].nlinked; ui++)
        {
          if (L[li].u_rec[ui].v->do_trend)
          {
            trend_reset (L[li].u_rec[ui].v, 1);
          }
        }
        if (L[li].model == HKY)
        {
          if (L[li].kappa_rec->v->do_trend)
          {
            trend_reset (L[li].kappa_rec->v, 1);
          }
        }
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          if (L[li].a_rec->v->do_trend == 1)
          {
            trend_reset (L[li].a_rec->v, 1);
          }
        }
      }
    }
    if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
        || modeloptions[SINGLEPOPULATION] == 1)
    {
      /* no split time */
    }
    else
    {
      for (k = 0; k < lastperiodnumber; k++)
        if (T[k].v->do_trend)
          trend_reset (T[k].v, 1);
    }

/* ADD ADDITONAL trend_reset() calls here */
  }
  if (init == 1)
  {
    trendspot = 0;
    recordtrendinc = 1;
    recordinc = 0;
    movespot = TRENDLASTPT;
    init = 0;
  }
  recordinc++;
  if (recordinc == recordtrendinc)
  {
    if (runoptions[LOADRUN])
    {

      if (lpgpd_v->do_trend)
        inctrend (movespot, trendspot, lpgpd_v,
                  gsampinf[loadarrayj][gsamp_pdgp] +
                  gsampinf[loadarrayj][gsamp_plgp] +
                  gsampinf[loadarrayj][gsamp_probgp]);
      for (j = 0; j < lastperiodnumber; j++)
        if (T[j].v->do_trend)
          inctrend (movespot, trendspot, T[j].v,
                    gsampinf[loadarrayj][gsamp_tp + j]);
    }
    else
    {
      if (lpgpd_v->do_trend)
      {
        inctrend (movespot, trendspot, lpgpd_v,
                  C[0]->allpcalc.probg + C[0]->allpcalc.pdg +
                  C[0]->allpcalc.plg);
      }

      if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
          || modeloptions[SINGLEPOPULATION] == 1)
      {
        /* no split time */
      }
      else
      {
        for (j = 0; j < lastperiodnumber; j++)
          if (T[j].v->do_trend)
            inctrend (movespot, trendspot, T[j].v, C[0]->tvals[j]);
      }
    }
    for (li = 0; li < nloci; li++)
    {
      if (runoptions[LOADRUN] == 0)
      {
        for (ui = 0; ui < L[li].nlinked; ui++)
          if (L[li].u_rec[ui].v->do_trend)
            inctrend (movespot, trendspot, L[li].u_rec[ui].v,
                      C[0]->G[li].uvals[ui]);
        if (L[li].model == HKY)
          if (L[li].model == HKY)
            inctrend (movespot, trendspot, L[li].kappa_rec->v,
                      C[0]->G[li].kappaval);
        if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
        {
          inctrend (movespot, trendspot, L[li].a_rec->v, 
                    IMA_assignment2value (0, li));
        }
      }
    }
/* ADD ADDITONAL inctrend() calls here */
    if (movespot == TRENDLASTPT && trendspot == TRENDLASTPT)
    {
      movespot = 0;
      recordtrendinc *= 2;
    }
    else
    {
      movespot += (movespot < TRENDLASTPT);
    }
    trendspot += (trendspot < TRENDLASTPT);
    recordinc = 0;
  }
  trenddoublepoint = movespot;
}                               /* trendrecord */

/* calculates the bin number of an xy array of a value_record that a value falls in,  increments the count in that bin */
void recordval (struct value_record *v, double val)
{
  int k;
  double logval;
  if (v->do_logplot)
  {
    assert (!(val < 0.0));
    logval = log (val);
    k =
      (int) (GRIDSIZE * (logval + v->plotrange.max) /
             (2.0 * v->plotrange.max));
  }
  else
  {
    k = (int) (GRIDSIZE * val / v->plotrange.max);
  }
  if (k < 0)
  {
    v->beforemin++;
  }
  else if (k >= GRIDSIZE)
  {
    v->aftermax++;
  }
  else
  {
    assert (!(k < 0));          // FIXME: it's been crashing
    v->xy[k].y++;
  }
  return;
}

/* this is used to record the names of loci and gene copies that migrate, from to to
write these names to a file */
void record_migration_names(void)
 {
  int i, j, li;
  int from, to;
  struct edge *gtree;
  from = migrationnamefrom;
  to = migrationnameto;
  for (li = 0; li < nloci; li++)
  {
    gtree = C[0]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > 0)
      {
        if (from == nowedgepop (0, &gtree[i], gtree[i].mig[j].mt) && to == C[0]->G[li].gtree[i].mig[j].mp)
        {
          fprintf(migrationnamefile,"%s ",L[li].name);
          if (i< L[li].numgenes)
            fprintf(migrationnamefile,"%s ",L[li].gNames[i]);
          else
            fprintf(migrationnamefile,"internal ");
        }
        j++;
      }
    }
  }
  fprintf(migrationnamefile,"\n");
 } // record_migration_names
/* to record a numerical value from the markov chain:
----------------------------------------------------
this works on instances of struct value_record

the value_record must be initiatlized (e.g. in initialize.c,  see e.g. init_g_rec)
this includes a call to init_value_record()

insert line(s) code into record()  below,  to make a call to recordval() 

*/

void record_migrations (void)
{
  int i, j, k, li, from, to, foundparam;
  struct edge *gtree;
  for (j = 0; j < nloci + (nloci > 1); j++)
    for (i = 0; i < nummigrateparams; i++)
    {
      migcount[j][i] = 0;
    }
  for (li = 0; li < nloci; li++)
  {
    gtree = C[0]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > 0)
      {
        from = nowedgepop (0, &gtree[i], gtree[i].mig[j].mt);
        to = C[0]->G[li].gtree[i].mig[j].mp;
        k = 0;
        foundparam = 0;
        while (k < nummigrateparams && foundparam == 0)
        {
          if (modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS])
          {
            if ((from == migfrom[k] && to == migto[k])
                || (from == migto[k] && to == migfrom[k]))
              foundparam = 1;
            else
              k++;
          }
          else
          {
            if (from == migfrom[k] && to == migto[k])
              foundparam = 1;
            else
              k++;
          }
        }
        assert (foundparam);
        if (nloci == 1)
        {
          migcount[0][k]++;
          recordval (&migration_counts_times[0][2 * k], gtree[i].mig[j].mt);
        }
        else
        {
          migcount[li + 1][k]++;
          migcount[0][k]++;
          recordval (&migration_counts_times[0][2 * k], gtree[i].mig[j].mt);
          recordval (&migration_counts_times[li + 1][2 * k],
                     gtree[i].mig[j].mt);

        }
        j++;
      }
    }
  }
  for (i = 0; i < nloci + (nloci > 1); i++)
    for (k = 0; k < nummigrateparams; k++)
    {
      recordval (&migration_counts_times[i][2 * k + 1], (double) migcount[i][k]);
    }
}                               //record_migrations

void record (void)
{
  int j, li, ui;
  struct genealogy *G;
#ifdef IMPROFILE
  startclock (&clock_record);
#endif
  if (outputoptions[PRINTTMRCA])
    for (li = 0; li < nloci; li++)
    {
      recordval (L[li].g_rec->v, C[0]->G[li].roottime);
    }
  for (li = 0; li < nloci; li++)
  {
    G = &(C[0]->G[li]);
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
      recordval (L[li].u_rec[ui].v, G->uvals[ui]);
    }
    if (L[li].model == HKY)
    {
      recordval (L[li].kappa_rec->v, G->kappaval);
    }
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || modeloptions[SINGLEPOPULATION] == 1)
  {
    /* we could have lastperiodnumber be 0 */
  }
  else
  {
    for (j = 0; j < lastperiodnumber; j++)
    {
      assert (C[0]->tvals[j] > T[j].pr.min && C[0]->tvals[j] < T[j].pr.max);
      recordval (T[j].v, C[0]->tvals[j]);
    }
  }
  if (outputoptions[MIGRATEHIST])
  {
    record_migrations ();
  }

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
  {
    /* no split time */
  }
  else
  {
    if (npops > 2 && npops < 5 && outputoptions[PRINTJOINTTEST])
      setup_multi_t_arrays ();
  }
  if (!outputoptions[DONTPRINTASCIITREND])
  {
    trendrecord (-1);
  }
#ifdef IMPROFILE
  addclock (&clock_record);
#endif
  return;
}                               /* record */

void savetreeinfo (void)        // use floats to save space
{
  int i;
  if (treessaved == 0)
  {
    if (cdurationmode == TIMESTEPS)
    {
      gsampinf = (float **) malloc (treestosave * sizeof (float *));
      for (i = 0; i < treestosave; i++)
        gsampinf[i] = (float *) malloc (gsampinflength * sizeof (float));
      memfortreessaved = treestosave;
    }
    else                        /* cdurationmode == TIMEINF */
    {
      gsampinf = (float **) malloc (MAXTREESTOSAVE * sizeof (float *));
      for (i = 0; i < MAXTREESTOSAVE; i++)
        gsampinf[i] = (float *) malloc (gsampinflength * sizeof (float));
      memfortreessaved = MAXTREESTOSAVE;
    }
  }
  if (treessaved >= (MAXTREESTOSAVE - 1) && cdurationmode == TIMEINF)
  {
    printf (" maximum possible trees saved \n");
    maxedouttreesave = 1;
  }
  else
  {
    savegsampinf (gsampinf[treessaved]);
  }
}                               /* savetreeinfo */

void loadtreevalues (void)
{
  char filenamewildcard[FNSIZE];
  char *ctp, *textline, *dataline, *c, tempc;
  int charspervalue = 12;
  FILE *sfile;
  int i, j, numtrees, totalnumtrees;
  int filefound, nofile;
  char defaultdir[2] = ".";
  struct dirent *dir_entry;
  DIR *dp;
  int numfiles = 0;
  int numtoload, nottoload, loadall, loaded, notloaded;
  float load_notload_ratio;
  if (strlen (loadfilebase))
  {
    strcpy (filenamewildcard, loadfilebase);
  }
  else
  {
    strcpy (filenamewildcard, outfilename);
    strtrunc (filenamewildcard, '.');
    strtrunc (filenamewildcard, '-');
  }
  strcat (filenamewildcard, "*.ti");
  SP "Base filename for loading files with sampled genealogies: %s\n", filenamewildcard);
  SP "Files loaded with sampled genealogies:\n");
  numtreefiles = 0;
textline = (char *) malloc (300 * sizeof (char));
dataline = (char *) malloc (gsampinflength * charspervalue * sizeof (char));
  ctp = &textline[0];
  numtrees = totalnumtrees = 0;
  if ((dp = opendir (defaultdir)) == NULL)
  {
    IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
  }
// first run through files to see how many trees there are in them 
  do
  {
    do
    {
      nofile = ((dir_entry = readdir (dp)) == NULL);
      if (!nofile)
        filefound = IsWildcardMatch (filenamewildcard, dir_entry->d_name, 0);
      else
        filefound = 0;
      numfiles += filefound != 0;
    } while (!nofile && !filefound);
    if (numfiles == 0 && nofile == 1)
      IM_err (IMERR_TIFILE, " not .ti files found");
    if (!nofile)
    {
      numtreefiles++;
      if ((sfile = fopen (dir_entry->d_name, "r")) == NULL)
      {
        IM_err (IMERR_TIFILE, " cannot open .ti file");
      }
      while (fgets (textline, 300, sfile)
             && strstr (textline, "VALUESSTART") == NULL && !feof (sfile));
      numtrees = 0;
      while ((tempc = fgetc (sfile)) != EOF)    // count lines
      {
        numtrees += (tempc == '\n');
      }
      if (numtrees < 1)
      {
        printf ("*no trees loaded from file %s\n", dir_entry->d_name);
        SP "*no trees loaded from file %s\n", dir_entry->d_name);
      }
      else
      {
        printf ("loaded %d genealogies from tree file  %s\n", numtrees,
                dir_entry->d_name);
        SP "loaded %d genealogies from tree file  %s\n", numtrees,
          dir_entry->d_name);
      }
      //fclose(sfile);
      totalnumtrees += numtrees;
    }
  } while (!nofile);
  if (treestosave > 0)
    numtoload = IMIN (totalnumtrees, treestosave);
  else
    numtoload = totalnumtrees;
  numtoload = IMIN (numtoload, MAXTREESTOSAVE);
  memfortreessaved = numtoload;
  gsampinf = (float **) malloc (numtoload * sizeof (float *));
  loadall = numtoload >= totalnumtrees;
  nottoload = totalnumtrees - numtoload;
  load_notload_ratio = (float) numtoload / (float) nottoload;
  loaded = 0;
  notloaded = 0;
  SP "\nHEADER INFORMATION FROM FIRST TREE FILE\n");
  SP "========================================\n");
// now go through again and save the trees
  closedir (dp);
  if ((dp = opendir (defaultdir)) == NULL)
  {
    IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
  }
  numfiles = 0;
  do
  {
    do
    {
      nofile = ((dir_entry = readdir (dp)) == NULL);
      if (!nofile)
        filefound = IsWildcardMatch (filenamewildcard, dir_entry->d_name, 0);
      else
        filefound = 0;
      numfiles += filefound != 0;
    } while (!nofile && !filefound);
    if (numfiles == 0 && nofile == 1)
      IM_err (IMERR_TIFILE, " not .ti files found");
    if (!nofile)
    {
      if ((sfile = fopen (dir_entry->d_name, "r")) == NULL)
      {
        IM_err (IMERR_TIFILE, " cannot open .ti file");
      }
      while (fgets (textline, 300, sfile) && strstr (textline, "VALUESSTART") == NULL && !feof (sfile))
        if (numfiles == 1)
        {
          SP "||%s", textline);
        }
      while (fgets (dataline, gsampinflength * charspervalue, sfile)
             && !feof (sfile))
      {
        if (loadall || (loaded == 0)
            || ((float) loaded / (float) notloaded) <= load_notload_ratio)
        {
          gsampinf[loaded] = (float *) malloc (gsampinflength * sizeof (float));
          c = dataline;
          for (i = 0; i < gsampinflength; i++)
          {
            sscanf (c, "%f", &gsampinf[loaded][i]);
            c = nextwhite (c);
          }
          loaded++;
        }
        else
          notloaded++;
      }
      fclose (sfile);
    }
  } while (!nofile);
  SP "END OF HEADER INFORMATION FROM FIRST TREE FILE\n");
  SP "==============================================\n\n");
  
  SP "   # of files: %d    total number of loaded genealogies: %d out of a total of: %d \n", numtreefiles, loaded, totalnumtrees);
  if (numtrees < 1)
    IM_err (IMERR_TIFILE, " no trees loaded from .ti file(s)");
  closedir (dp);
  treessaved = loaded;
  for (j = 0; j < treessaved; j++)
  {
    // use full range of t, ignore t.pr.min > 0
    for (i = 0; i < lastperiodnumber; i++)
    {
      recordval (T[i].v, gsampinf[j][gsamp_tp + i]);

    }
    if (!outputoptions[DONTPRINTASCIITREND])
      trendrecord (j);
  }
  XFREE (textline);
  XFREE (dataline);
}                               /* loadtreevalues */

void printsteps (FILE * outto, double like)
{
  int ci;

  ci = 0;
  if (!burndone)
  {
    fprintf (outto, "=BURNIN-PERIOD===============================\n");
    fprintf (outto, "STEP # %d  p(D|G): %.3lf p(G): %.3lf\n", step,
             like, C[ci]->allpcalc.probg);
  }
  else
  {
    fprintf (outto, "=============================================\n");
    if (treessaved > 0)
      fprintf (outto,
               "STEP # %ld # Genealogies Saved: %d p(D|G): %.1lf p(G): %.1f\n",
               step - burnsteps, treessaved, like, C[ci]->allpcalc.probg);
    else
      fprintf (outto, "STEP # %ld  p(D|G): %.3lf p(G): %.3lf\n",
               step - burnsteps, like, C[ci]->allpcalc.probg);
  }
  return;
}

/* To print acceptance rates:
	reclist[] is an array of pointers to struct chainstate_record_updates_and_values
	set values of reclist[] to those structures for which you want to print acceptance rates
	call printacceptancerates ()
	*/

void callprintacceptancerates (FILE * outto)
{
  int i, j, li;
  // length of this array must be fairly long, although it is technically possible to have MAXLOCI * MAXLINKED records,  but very unlikely
  struct chainstate_record_updates_and_values *reclist[MAXLOCI + MAXLINKED];
// t values 
  for (i = 0; i < numsplittimes; i++)
    reclist[i] = (T + i);

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || modeloptions[SINGLEPOPULATION] == 1)
  {
    assert (numsplittimes == 0);
    /* no split time */
  }
  else
  {
    printacceptancerates (outto, numsplittimes, reclist,
                          "Update Rates -- Population Splitting Times");
  }
// genealogy updates
  for (li = 0; li < nloci; li++)
    reclist[li] = L[li].g_rec;
  printacceptancerates (outto, nloci, reclist, "Update Rates -- Genealogies");
// mutation rate scalars
  if (nurates > 1
      && (runoptions[PRINTMUTATIONUPDATESTOSCREEN] || outto != stdout))
  {
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        reclist[i] = &L[li].u_rec[j];
        i++;
      }
// kappa values for HKY model
    printacceptancerates (outto, i, reclist,
                          "Update Rates -- Mutation Rate Scalars");
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        if (L[li].umodel[j] == HKY)
        {
          reclist[i] = L[li].kappa_rec;
          i++;
        }
      }
    if (i > 0)
      printacceptancerates (outto, i, reclist,
                            "Update Rates -- HKY Model Kappa parameter");
// STR ancestral allele states 
    for (i = 0, li = 0; li < nloci; li++)
      for (j = 0; j < L[li].nlinked; j++)
      {
        if (L[li].umodel[j] == STEPWISE)
        {
          reclist[i] = &L[li].A_rec[j];
          i++;
        }
      }
    if (i > 0)
    {
      printacceptancerates (outto, i, reclist,
                            "Update Rates -- STR Genealogy Allele States");
    }
  }
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (li = 0; li < nloci; li++)
      reclist[li] = L[li].a_rec;

    if (assignmentoptions[POPULATIONASSIGNMENTLOCAL] == 1)
    {
      printacceptancerates (outto, nloci, reclist,
                            "Update Rates -- Assignment");
    }
    // assignment updating information 
    // set reclest
    // call printacceptancerates 
    //
    printacceptancerates_multichain (outto);
  }

  return;
}                               //callprintacceptancerates


/* set up arrays pointing to information to put in curve plots,  then call asciicurve
   some things that are plotted are based on struct value_record and others on struct i_param
   this is why we cannot simply call asciicurve() with a pointer to a single type of structure  */
void callasciicurves (void)
{
  struct plotpoint **curvexy;
  char **curvestr;
  int *curve_do_logplot;
  int numcurve = 0;
  int i, j, li, ui;
  int *nrecstep;
  //_CrtCheckMemory( );
// find out how many curves
  numcurve += modeloptions[SPLITTINGRATEPARAMETER];
  numcurve += numpopsizeparams;
  for (i = 0; i < nummigrateparams; i++)
    if (imig[i].pr.max > MPRIORMIN)
      numcurve++;
  numcurve += numsplittimes;
  if (runoptions[LOADRUN] == 0 && nurates > 1)
    numcurve += nurates;
  if (runoptions[LOADRUN] == 0)
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY)
        numcurve++;
// allocate
  curvexy = (struct plotpoint **)malloc (numcurve * sizeof (struct plotpoint *));
  curvestr = (char **) malloc (numcurve * sizeof (char *));
  nrecstep = (int *) malloc (numcurve * sizeof (int));
  curve_do_logplot = (int *) malloc (numcurve * sizeof (int));
// assign
  j = 0;
  if (modeloptions[SPLITTINGRATEPARAMETER] == 1)
  {
    curvexy[j] = isplit.xy;
    curvestr[j] = &isplit.str[0];
    curve_do_logplot[j] = 0;
    nrecstep[j] = 1;
    j++;
  }
  for (i = 0; i < numpopsizeparams; i++)
  {
    curvexy[j] = itheta[i].xy;
    curvestr[j] = &itheta[i].str[0];
    curve_do_logplot[j] = 0;
    nrecstep[j] = 1;
    j++;
  }
  for (i = 0; i < nummigrateparams; i++)
    if (imig[i].pr.max > MPRIORMIN)
    {
      curvexy[j] = imig[i].xy;
      curvestr[j] = &imig[i].str[0];
      curve_do_logplot[j] = 0;
      nrecstep[j] = 1;
      j++;
    }
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || modeloptions[SINGLEPOPULATION] == 1)
  {
    /* no split time */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
    {
      curvexy[j] = T[i].v->xy;
      curvestr[j] = &T[i].v->str[0];
      curve_do_logplot[j] = T[i].v->do_logplot;
      nrecstep[j] = recordstep;
      j++;
    }
  }
  if (runoptions[LOADRUN] == 0 && nurates > 1)
  {
    for (li = 0; li < nloci; li++)
      for (ui = 0; ui < L[li].nlinked; ui++)
      {
        curvexy[j] = L[li].u_rec[ui].v->xy;
        curvestr[j] = &L[li].u_rec[ui].v->str[0];
        curve_do_logplot[j] = L[li].u_rec[ui].v->do_logplot;
        nrecstep[j] = recordstep;
        j++;
      }
  }
  if (runoptions[LOADRUN] == 0)
    for (li = 0; li < nloci; li++)
      if (L[li].model == HKY)
      {
        curvexy[j] = L[li].kappa_rec->v->xy;
        curvestr[j] = &L[li].kappa_rec->v->str[0];
        curve_do_logplot[j] = L[li].kappa_rec->v->do_logplot;
        nrecstep[j] = recordstep;
        j++;
      }
  assert (numcurve == j);
  for (j = 0; j < numcurve; j++)
    asciicurve (outfile, curvexy[j], curvestr[j], curve_do_logplot[j],
                nrecstep[j]);
//free
  XFREE (curvexy);
  XFREE (curvestr);
  XFREE (curve_do_logplot);
  XFREE (nrecstep);
}                               //callasciicurve 

// makes calls to asciitrend
void callasciitrend (FILE * outtofile)
{
  int i, li, ui;
  asciitrend (outtofile, lpgpd_v, trenddoublepoint, trendspot);

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1
      || modeloptions[SINGLEPOPULATION] == 1)
  {
    /* no split time */
  }
  else
  {
    for (i = 0; i < lastperiodnumber; i++)
      asciitrend (outtofile, T[i].v, trenddoublepoint, trendspot);
  }
  if (nurates > 1 && runoptions[LOADRUN] == 0)
  {
    for (li = 0; li < nloci; li++)
      for (ui = 0; ui < L[li].nlinked; ui++)
        asciitrend (outtofile, L[li].u_rec[ui].v, trenddoublepoint, trendspot);
  }
  for (li = 0; li < nloci; li++)
    if (L[li].model == HKY && runoptions[LOADRUN] == 0)
      asciitrend (outtofile, L[li].kappa_rec->v, trenddoublepoint, trendspot);

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    for (li = 0; li < nloci; li++)
    {
      asciitrend (outtofile, L[li].a_rec->v, trenddoublepoint, trendspot);
    }
  }
  return;
}                               // callasciitrend 

void 
printoutput (void)         // mostly calls functions in output.c
{
  int i, ci;
  double seconds;
  static int lasttreesaved = -1;
  int p;
  float *holdpeakloc;
  double multitpeak[MAXPOPS - 1];

  if (runoptions[LOADRUN] == 0)
  {
    /* save tree info in *.ti file */
    if (!runoptions[DONTSAVETREES] && treessaved > 0)
      savetreefile (treeinfosavefilename, treeinfosavefile, &lasttreesaved, gsampinflength);
    if (cdurationmode == TIMEINF)
    {
      remove (oldoutfilename);
      rename (outfilename, oldoutfilename);
    }
  }
  if ((outfile = fopen (outfilename, "w")) == NULL)
  {
    IM_err (IMERR_CREATEFILEFAIL, "Error opening text file for writing");
  }
  printrunbasics (outfile, runoptions[LOADRUN], fpstr, burnsteps, recordint,
                  recordstep, savetreeint, endtime, starttime, hilike,
                  hiprob, step);
  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    print_num_genes_pops (outfile);
  }
  fprintf (outfile, "Average loglikelihood: %lf\n", gloglikelihood);
  if (runoptions[LOADRUN] == 0)
  {
    callprintacceptancerates (outfile);
    if (numchains > 1)
      printchaininfo (outfile, heatmode, adaptcheck, hval1, hval2);
    callprintautoctable (outfile, step);
  }
  if (outputoptions[PARAMGREATERTHAN])
  {
    print_greater_than_tests (outfile);
  }
  if (!modeloptions[EXPOMIGRATIONPRIOR])        // have to do the math for case of migration with exponential prior
    print_means_variances_correlations (outfile);
  init_surface_calc ();
/*  get marginal peaks */
  if (!calcoptions[DONTCALCLIKELIHOODMUTATION])
  {
    p =
      numpopsizeparams + nummigrateparams +
      modeloptions[SPLITTINGRATEPARAMETER];
    holdpeakloc = (float *) malloc (p * sizeof (float));
    printf ("surface calculations . . .\n");
    findmarginpeaks (outfile, holdpeakloc);
/* not implemented
    likecalc = (double) - jointp(holdpeakloc - 1); // have to use -1 because jointp expects the first value in position 1,  this in turn in case jointp() called by NR functions
    FP"\n  Margin peak :");
    for (i=0;i< p;i++)
        FP"\t%.4f",holdpeakloc[i]);
    FP"\n");
    FP"Posterior probability at margin peak : %.4lf \n",likecalc); 
    closeopenout (&outfile, outfilename);  */

  }

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
  {
    print_num_genes_pops (outfile);
  }
  else
  {
/* get joint splittime peak */
    if (npops > 2 && npops < 5 && !runoptions[LOADRUN] && outputoptions[PRINTJOINTTEST])
    {
      return_joint_t (multitpeak);
      FP "\nEstimated joint splitting time from multi-dimensional histogram\n");
      FP "  number of bins per dimension %d\n", NUMTARRAYBINS);

      if (modeloptions[SPLITTINGRATEPARAMETER])
        FP "  highest t represented in bin array: %.3lf\n",
          TIMERECORDPRIORFRAC * T[lastperiodnumber - 1].pr.max);
      else
        FP "  highest t represented in bin array: %.3lf\n",
          T[lastperiodnumber - 1].pr.max);
      FP "  Posterior probability of estimated joint value of splitting time: %7.4lf\n", joint_t_prob (&multitpeak[0]));
      FP "---------------------------------------------------------------\n");
      for (i = 0; i < numsplittimes; i++)
        FP "   %s\t%.3lf\n", T[i].str, multitpeak[i]);
      FP "\n\n");
    }
  }
  printhistograms (outfile, recordstep, generationtime, scaleumeaninput);
  ci = 0;
  if (!outputoptions[DONTPRINTASCIICURVE])
  {
    FP "\n\nASCII Curves - Approximate Posterior Densities \n");
    FP "===================================================\n");
    callasciicurves ();
  }
  if (!outputoptions[DONTPRINTASCIITREND])
  {
    ci = 0;
    if (trendspot <= 1)
    {
      FP "burn period too short to plot trends,  trend recording begins at step %d \n",burntrendstartdelay); 
    }
    else
    {
      FP "\n\nASCII Plots of Parameter Trends \n");
      FP "===================================================\n");
      FP " - note points to the left of '!' on X axis have twice the density in time relave to points to the right\n\n");
      callasciitrend (outfile);
    }

  }
  time (&endtime);
  seconds = difftime (endtime, starttime);
  FP "Time Elapsed : %d hours, %d minutes, %d seconds \n\n",
    (int) seconds / (int) 3600,
    ((int) seconds / (int) 60) - ((int) 60 * ((int) seconds / (int) 3600)),
    (int) seconds - (int) 60 *((int) seconds / (int) 60));
  FP "\nEND OF OUTPUT\n");
#ifdef IMPROFILE
  print_timings (outfile);
#endif
  f_close (outfile);
  free_surface_calc ();

  if (outputoptions[MIGRATEHIST] && nummigrateparams > 0)
  {

    if ((migplotfile = fopen (migplotfilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL,
              "Error opening file for plotting migration amounts and times");
    }
    printmigrationhistograms (migplotfile, recordstep);
    f_close (migplotfile);
  }


  if (!calcoptions[DONTCALCLIKELIHOODMUTATION])
    XFREE (holdpeakloc);
  if (runoptions[SAVEMCSTATEFILE])
  {
    writemcf (mcfwritefilename);
  }
  if (jheyoptions[WRITEMIGRATIONNAME])
  {
    f_close(migrationnamefile);
    migrationnamefile = fopen (migrationnamefilename, "a");
  }

  return;
}                               /* printoutput */


void intervaloutput (FILE * outto)
{
  int li;
  int j;
  double seconds;
  double like;
  int ci;
  ci = 0;
  checkhighs (0, printint, &hilike, &hiprob, &like, step);
  if (((step / (int) printint) * (int) printint == step && step > 0)
      || outto != stdout)
  {
    printsteps (outto, like);

    callprintacceptancerates (outto);

    printcurrentvals (outto);
    callprintautoctable (outto, step);
    if (numchains > 1)
      printchaininfo (outto, heatmode, adaptcheck, hval1, hval2);
    if (treessaved > 0)
    {
      time (&remained_endtime);
      seconds = difftime (remained_endtime, remained_starttime);
    }
    
    /* For ASSIGNMENT */
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      printsteps (outto, like);
      print_num_genes_pops (outto);
      for (ci = 0; ci < numchains; ci++)
      {
        for (j = 0; j < Cupinf[ci].num_uptypes; j++)
        {
          Cupinf[ci].upinf[j].accp = 0;
          Cupinf[ci].upinf[j].tries = 0;
        }
      }
      for (li = 0; li < nloci; li++)
      {
        for (j = 0; j < L[li].a_rec->num_uptypes; j++)
          L[li].a_rec->upinf[j].accp = L[li].a_rec->upinf[j].tries = 0;
        if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
        {
          IMA_assigned_reset ();
        }
      }
    }

  }
  return;
}                               /* intervaloutput */

// check if it is time to call record() and savetreeinfo(), and call if it is
void check_to_record (void)
{
  static int i;
  static int j;
  static int init = 0;

  if (init == 0)
  {
    i = recordint;
    j = savetreeint;
    init = 1;
  }
  if (i == recordint)
  {
    record ();                  // record values of parameters that are in mcmc
    recordstep++;
    i = 1;
  }
  else
  {
    i++;
  }

  if (j == savetreeint)
  {
    savetreeinfo ();            // record values associated with genealogies
    if (jheyoptions[WRITEMIGRATIONNAME])
      record_migration_names();
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      recordassignmentloci (outfilename, 0);
    }
    treessaved++;
    j = 1;
  }
  else
  {
    j++;
  }
}                               // check_to_record

int run_main (int argc, char *argv[])
{
#ifdef HPDBG
  int tmpFlag = _CrtSetDbgFlag (_CRTDBG_REPORT_FLAG);
//  tmpFlag |= _CRTDBG_DELAY_FREE_MEM_DF;  // slow
  //tmpFlag |= _CRTDBG_CHECK_ALWAYS_DF;
  _CrtSetDbgFlag (tmpFlag);
  // Send all reports to STDOUT
  _CrtSetReportMode (_CRT_WARN, _CRTDBG_MODE_FILE);
  _CrtSetReportFile (_CRT_WARN, _CRTDBG_FILE_STDOUT);
  _CrtSetReportMode (_CRT_ERROR, _CRTDBG_MODE_FILE);
  _CrtSetReportFile (_CRT_ERROR, _CRTDBG_FILE_STDOUT);
  _CrtSetReportMode (_CRT_ASSERT, _CRTDBG_MODE_FILE);
  _CrtSetReportFile (_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
  _CrtCheckMemory ();
  //heapdebugfile= fopen("debugheap.out","w");
  //fprintf(heapdebugfile,"allocate \n");

#endif //HPDBG
  init_IMA ();
  start (argc, argv);
  if (runoptions[LOADRUN])
  {
    loadtreevalues ();
    recordstep = treessaved;
    if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      IMA_align_genealogy (outfilename, treessaved);
    }
    printoutput ();
    free_ima_main_stuff ();
    return 0;
  }

  printf ("Starting Markov chain.\n");
  init_after_start_IMA ();
  step = 0;
  recordstep = 0;
  while (run ())
  {
    qupdate ();
    if (burndone)
    {
      check_to_record ();
#ifdef COUNT_PRIOR_ASSIGNMENT
      IMA_rgf_tickSaasn (0);
#endif /* COUNT_PRIOR_ASSIGNMENT */
      gloglikelihood += C[0]->allpcalc.pdg; 
      gloglikelihood /= 2.0;
    }
    step++;

  }
#ifdef COUNT_PRIOR_ASSIGNMENT
  IMA_rgf_saveSaasn ();
#endif /* COUNT_PRIOR_ASSIGNMENT */

  if (assignmentoptions[POPULATIONASSIGNMENT] == 1
      && assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 0)
  {
    IMA_align_genealogy (outfilename, treessaved);
  }
  printoutput ();

  if (jheyoptions[WRITEMIGRATIONNAME])
    f_close(migrationnamefile);
  free_ima_main_stuff ();
  return 0;
//fclose(heapdebugfile);
}                               /* main */
