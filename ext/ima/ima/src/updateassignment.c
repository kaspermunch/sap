/* IMa2  2007-2009 Jody Hey, Rasmus Nielsen and Sang Chul Choi */

#include <search.h>
#undef GLOBVARS
#include "imamp.h"
#include "imagsl.h"
#include "update_gtree_common.h"
#include "utilities.h"
#include "updateassignment.h"

extern char outfilename[];
extern char *trueassignment;

static int snind;
static int snind_unknown;
static int *sngenes_ind;     
static int ssizePi;
static double **sassigned;
static int **sind2gi;        
static int sngenes;
static int sci;
static double *snnminus1;
static int *snsample_pop;
static int *slineages;    /* array of maximum lineages */
static int *smigrates;    /* array of maximum populations for a lineage */

static struct edgemiginfo saoem;
static struct edgemiginfo saosm;
static struct edgemiginfo sanem;
static struct edgemiginfo sansm;
static struct edgemiginfo saoem2;
static struct edgemiginfo sanem2;
static struct edgemiginfo *saoems;
static struct edgemiginfo *sanems;
static struct edgemiginfo *saosms;
static struct edgemiginfo *sansms;
static im_popntree saC;            /* _s_tatic _a_ssignment _C_hain */
static im_savedlocus *saT;         /* _s_tatic _a_ssignment _T_ree */

static im_ginfo *saGiold;
static im_ginfo *saGioldz1;
static im_ginfo *saGioldz2;
static im_ginfo *saGinew;
static im_ginfo *saGinewz1;
static im_ginfo *saGinewz2;
static im_bfupdate sbfold;
static im_bfupdate sbfnew;

static int snasn;
static int *saasn;

static void IMA_init_edgemiginfo ();
static void IMA_free_edgemiginfo ();
static int genename_cmp (const void *c1, const void *c2);
static int addpair (struct genename **gns, int lns, char *n, int li, int i);
static void freepair (struct genename *gns, int lns);
static void IMA_memory_initnnminus1 ();
static void IMA_memory_freennminus1 ();
static int IMA_number_largest_branches ();
static void   removeSpaces (char *s1, const char *s2);
static int IMA_genealogy_sister (int ci, int li, int ei);
static double IMA_genealogy_time (int ci, int li, int ei);
static int findperiodoftime (int ci, double ktime);
static int IMA_edge_labeloftime (int ci, int li, int ei, double ktime);
static int IMA_genealogy_nmigration (int ci, int li, int ei);
static int IMA_genealogy_nmigrationupto (int ci, int li, int ei, double t);
static void IMA_genenode_print_node (FILE *fp, struct edge *gtree, int e);
static void IMA_genenode_print (FILE *fp, struct edge *gtree, int e);
static void IMA_genenode_print_nodeseq (FILE *fp, struct edge *gtree, int e, int si);
void IMA_genenode_printseq (FILE *fp, struct edge *gtree, int e, int si);
static int    IMA_node_period (int ci, int li, int ei);
static int    IMAedge_downperiod (int ci, int li, int ei);
static double IMAedge_length (int ci, int li, int ei);
static int    IMAmig_period (int ci, int li, int ei, int mi);
static void   BitSPrintF (char *s, UByteP A, int size);
void print_popnlist ();
static double IMA_assignment2value_mean (int ci, int li);
static double relabel (int ci, int *mpart);
static void IMA_edge_fillmiginfo (int ci, int li, int ei, 
                                  struct edgemiginfo *em);
void IMA_edge_copynewmig (int ci, int li, struct edgemiginfo *em);
static void IMA_restore_edgemig (int ci, int li, int ei, im_migstruct *savedmig,
                                int ncmmsavedmig);
static double relabellocus (int ci, int li, int *mpart);
static int IMA_edge_movable (int ci, int li, int ei);
static double IMA_popn_move (int *topop, int ci, int pop, int apop);
static double addmigration_relabel (int ci, int li, int newpop, 
                                    int oldmigcount,
                                    int *newmigcount,
                                    struct edgemiginfo *olde, 
                                    struct edgemiginfo *newe);
static void copy_probcalc_all (struct probcalc *dest, struct probcalc *srce);
static double bfmove (int ci);
static void bfmove_chooseA1 (int ci, int li, im_bfupdate *pbf);
static double bfmove_selectA1 (int ci, int li, im_bfupdate *pbf, double lq);
static int bfmove_chooseAnEdge (int ci, int li, im_bfupdate *pbf, double ct, int pc);
static int bfmove_selectAnEdge (int ci, int li, im_bfupdate *pbf, im_bfupdate *qbf);
static void IMA_sbf_correctnl (int ci, int li, int sis, int *n, int **l,
                               int *seqzu, int *seqzl);
static char IMA_sbf_correctn3 (int ci, int li, int sis,
                               int *seqzu, int *seqzl);
static double bfmove_nextI (int ci, int li,
                            double lq, im_event *ek, 
                            int np, double rates[], 
                            im_ginfo *gk, 
                            im_bfupdate *pbf,
                            im_bfupdate *qbf, /* FIXME: qbf is not used. */
                            int *ki,
                            int z1,
                            char *is_coalesced1,
                            int *pc1);
static int bfmove_ratesI (double (*r)[],
                          im_ginfo *gk,
                          im_ginfo *gkz1,
                          im_bfupdate *pbf,
                          int ki,
                          int pc1);
static void bfmove_ratesII (double (*rates)[], im_ginfo *gk, im_bfupdate *pbf,
                            int ki,
                            int pc1, int pc2, int pc3);
static void bfmove_ratesIII (double (*rates)[], im_ginfo *gk, im_bfupdate *pbf,
                             int ki,
                             int pc1, int pc3);
static void IMA_genealogy_detach (int ci, int li, im_bfupdate *pbf,
                                  int z1);
static int IMA_sbf_c (int ci, int li, int edge, int pi, 
                      int *seqz,
                      int *n, int **l, int *a);
static int IMA_sbf_m (int ci, int li, int edge, int pi, int pj,
                      int *n, int **l, int *a);
static int IMA_sbf_nsum (im_bfupdate *pbf, int ci, int li,
                         int pc1, int pc3);
static void IMA_sbf_lineages (im_bfupdate *pbf, int ci, int li, int z1);
static void IMA_intervals_collectsbfnew (int ci, int li, im_edgemiginfo *newm3);
static double bfq (im_bfupdate *pbf, im_bfupdate *qbf, 
                   int ci, int li, 
                   im_ginfo *gk,
                   im_ginfo *gkz1);
static double IMA_choose_migtopop (int *topop, int ci, int li, int ki, int pi, im_ginfo *gk);
static void IMA_edgemiginfo_copymig (im_edgemiginfo *em, int ci, int li, int ei);
static void IMA_edgemiginfo_copypartmig (im_edgemiginfo *em, 
                                         int ci, int li, int ei, double upt);
static void IMA_edgemiginfo_appendmig (im_edgemiginfo *em, int ci, int li, int ei);
static void IMA_genealogy_copymig (int ci, int li, int ei, im_edgemiginfo *em);
static void IMA_genealogy_appendmig (int ci, int li, int sis, im_edgemiginfo *em, 
                                     int edge, int down);
static int IMA_intervals_collect (int ci, int li, double upt, double dnt);
                                  
static int IMA_genealogy_derivetheta (char w, int ci, int li, int z1);
int IMA_ginfo_computetheta (im_ginfo *gk);
static void IMA_memory_initginfo (im_ginfo **g);
static void IMA_memory_freeginfo (im_ginfo **g);
static void IMA_memory_resetginfo (im_ginfo *g);
static void IMA_memory_initsaGi ();
static void IMA_memory_freesaGi ();
static int ran_discrete (int n, ...);
static int ran_discretea (int n, double pr[]); 
static int ran_discreteloga (int n, double pr[]); 
static int normalizeLoga (int n, double pr[]); 
int IMA_create_intervals ();
int IMA_delete_intervals ();
static int im_event_cmp (const void *a, const void *b);
int IMA_convertStr2Asn (int ncol, int *a, char *s);
int IMA_io_readp (char *fname, int ***m, int *mr, int *mc);
int IMA_gusfield_distanceRemoved (int *A, int **m, int nrow, int ncol, double *r);
int IMA_gusfield_distance (int *A, int **m, int nrow, int ncol); 
int IMA_gusfield_squaredDistance (int *A, int **m, int nrow, int ncol); 
double IMA_gusfield_varianceDistance (int *A, int **m, int nrow, int ncol, double mean); 
/*
static int skip_a_line (FILE *fp);
*/

int IMA_gdistanceRemoved (int n, int *Ai, int *Bi, char *Ci);
int IMA_gdistanceOrder (int n, int *Ai, int *Bi, int m, int *od);
int IMA_gdistance (int n, int *Ai, int *Bi);
int IMA_gdistanceSimple (int n, int *Ai, int *Bi);
int IMA_max_intarray (int n, int *A);

int IMA_permLexSuccessor (int n, int *p);
int IMA_permLexRank (int n, int *p);
int IMA_permLexUnrank (int *p, int n, int r);
int IMA_permFprintf (FILE *fp, int n, int *p);
int IMA_factorial (int n);
static int IMA_ptree_nsispop (int ci);
static int IMA_ptree_sispop (int ci, int **s);
static void IMA_savelocusinfo (int ci, int li);
static void IMA_restorelocusinfo (int ci, int li);
static void   IMA_genealogy_save_reset (int li);
static void   IMA_genealogy_save_init ();
static void   IMA_genealogy_save_free ();
static int    IMA_genealogy_save_edge (int ci, int li, int ei);
static void   IMA_genealogy_restore_edges (int ci, int li);
static void   IMA_genealogy_joinsisdown (int ci, int li, int ei);
static void   IMA_genealogy_splitsisdown (int ci, int li, int ei, int gparent, int gsister);
static char*  IMA_getIndividualName (int i);
static int    IMA_getGenotype (int ii, int li, int *ig1, int *ig2);
static void   IMA_set_popntree_nosplit ();
static void   IMA_set_popntree ();
static void   IMA_unset_popntree ();
/* assertgenealogy */
#ifndef NDEBUG
static void assertgenealogylocdiploid (int ci, int li);
static void assertgenealogyindividual (int ci);
static void assertgenealogyedgelabel (int ci, int li, int ei);
static void assertgenealogylabel (int ci);
static void assertgenealogyloclabel (int ci, int li);
static void assertgenealogyultrametric (int ci);
static void assertgenealogylocultrametric (int ci, int li);
static void assertgenealogylocbranch (int ci, int li);
static void assertgenealogybranch (int ci);
static void IMA_check_popn_up (int kprevperiod, int prevpop, int kperiod, int pop);
static void IMA_check_popn_move (int notpop, int kperiod, int pop);
#endif /* NDEBUG */
static void   IMA_ind_init ();
static void   IMA_ind_fin ();
static void   IMA_ind_set ();
static void   IMA_ind_edge_init_i ();
static void   IMA_ind_edge_copy_i ();
static int    IMA_ind_n (int k);
int IMA_rgf_setSnasn (int nind, int npops);
void
IMA_init_edgemiginfo ()
{
  int i;
  assert (!(npops < 0));
  saoem.mtimeavail = NULL;
  saoem.mp = NULL;
  saosm.mtimeavail = NULL;
  saosm.mp = NULL;
  sanem.mtimeavail = NULL;
  sanem.mp = NULL;
  sansm.mtimeavail = NULL;
  sansm.mp = NULL;
  saoem2.mtimeavail = NULL;
  saoem2.mp = NULL;
  sanem2.mtimeavail = NULL;
  sanem2.mp = NULL;
  
  saoem.mtimeavail = malloc (npops * sizeof (double));
  saoem.mp = malloc (npops * sizeof (int));
  saosm.mtimeavail = malloc (npops * sizeof (double));
  saosm.mp = malloc (npops * sizeof (int));
  sanem.mtimeavail = malloc (npops * sizeof (double));
  sanem.mp = malloc (npops * sizeof (int));
  sansm.mtimeavail = malloc (npops * sizeof (double));
  sansm.mp = malloc (npops * sizeof (int));
  saoem2.mtimeavail = malloc (npops * sizeof (double));
  saoem2.mp = malloc (npops * sizeof (int));
  sanem2.mtimeavail = malloc (npops * sizeof (double));
  sanem2.mp = malloc (npops * sizeof (int));

  if (saoem.mtimeavail == NULL 
      || saoem.mp == NULL
      || saosm.mtimeavail == NULL 
      || saosm.mp == NULL
      || sanem.mtimeavail == NULL 
      || sanem.mp == NULL
      || sansm.mtimeavail == NULL 
      || sansm.mp == NULL
      || saoem2.mtimeavail == NULL 
      || saoem2.mp == NULL
      || sanem2.mtimeavail == NULL 
      || sanem2.mp == NULL)
    {
      IM_err (IMERR_MEM, "saoem/saosm/sanem/sansm/saoem2/sanem2");
    }

  saoems = malloc (2 * nloci * sizeof (struct edgemiginfo)); 
  sanems = malloc (2 * nloci * sizeof (struct edgemiginfo)); 
  saosms = malloc (2 * nloci * sizeof (struct edgemiginfo)); 
  sansms = malloc (2 * nloci * sizeof (struct edgemiginfo)); 
  for (i = 0; i < 2 * nloci; i++)
    {
      saoems[i].mp = malloc (npops * sizeof (int));
      saoems[i].mtimeavail = malloc (npops * sizeof (double));
      sanems[i].mp = malloc (npops * sizeof (int));
      sanems[i].mtimeavail = malloc (npops * sizeof (double));
      saosms[i].mp = malloc (npops * sizeof (int));
      saosms[i].mtimeavail = malloc (npops * sizeof (double));
      sansms[i].mp = malloc (npops * sizeof (int));
      sansms[i].mtimeavail = malloc (npops * sizeof (double));
    }
  return;
}

void
IMA_free_edgemiginfo ()
{
  int i;

  free (saoem.mtimeavail);
  saoem.mtimeavail = NULL;
  free (saoem.mp);
  saoem.mp = NULL;
  free (saosm.mtimeavail);
  saosm.mtimeavail = NULL;
  free (saosm.mp);
  saosm.mp = NULL;
  free (sanem.mtimeavail);
  sanem.mtimeavail = NULL;
  free (sanem.mp);
  sanem.mp = NULL;
  free (sansm.mtimeavail);
  sansm.mtimeavail = NULL;
  free (sansm.mp);
  sansm.mp = NULL;

  free (saoem2.mtimeavail);
  saoem2.mtimeavail = NULL;
  free (saoem2.mp);
  saoem2.mp = NULL;
  free (sanem2.mtimeavail);
  sanem2.mtimeavail = NULL;
  free (sanem2.mp);
  sanem2.mp = NULL;

  for (i = 0; i < 2 * nloci; i++)
    {
      free (saoems[i].mp);
      saoems[i].mp = NULL;
      free (saoems[i].mtimeavail);
      saoems[i].mtimeavail = NULL;
      free (sanems[i].mp);
      sanems[i].mp = NULL;
      free (sanems[i].mtimeavail);
      sanems[i].mtimeavail = NULL;
      free (saosms[i].mp);
      saosms[i].mp = NULL;
      free (saosms[i].mtimeavail);
      saosms[i].mtimeavail = NULL;
      free (sansms[i].mp);
      sansms[i].mp = NULL;
      free (sansms[i].mtimeavail);
      sansms[i].mtimeavail = NULL;
    }
  free (saoems);
  saoems = NULL;
  free (sanems);
  sanems = NULL;
  free (saosms);
  saosms = NULL;
  free (sansms);
  sansms = NULL;

  return;
}

int
genename_cmp (const void *c1, const void *c2)
{
  const struct genename *C1 = (const struct genename *) c1;
  const struct genename *C2 = (const struct genename *) c2;
  return strcmp (C1->name, C2->name);
}

int
addpair (struct genename **gns, int lns, char *n, int li, int ei)
{
  int v;
  struct genename target, *result;
  assert (!(ei < 0));
  target.name = n;
  if (lns > 0 && *gns != NULL)
    {
      qsort (*gns, (size_t) lns, sizeof (struct genename), genename_cmp);
      result =
        bsearch (&target, *gns, (size_t) lns, sizeof (struct genename), genename_cmp);
    }
  else
    {
      result = NULL;
    }
  if (result)
    {
      assert (!(result->index < 0));
      assert (L[li].pairs[ei] < 0);
      assert (L[li].pairs[result->index] < 0);
      L[li].pairs[ei] = result->index;
      L[li].pairs[result->index] = ei;
      v = lns;
    }
  else
    {
      /* create one if there is none, */
      if (*gns == NULL)
        {
          *gns = (struct genename *) malloc (sizeof (struct genename));
          (*gns)->name = strdup (n);
          (*gns)->index = ei;
          assert (ei == 0);
          v = 1;
        }
      /* increase one if it is not empty, */
      else
        {
          *gns = realloc (*gns, (lns + 1) * sizeof (struct genename));
          v = lns + 1;
          (*gns)[lns].name = strdup (n);
          (*gns)[lns].index = ei;
        }
    }
  return v;
}

void
freepair (struct genename *gns, int lns)
{
  int i;

  if (gns != NULL)
    {
      for (i = 0; i < lns; i++)
        {
          XFREE (gns[i].name);
          gns[i].name = NULL;
        }
      XFREE (gns);
      gns = NULL;
    }
  return;
}
void
IMA_genealogy_pairing ()
{
  /* can be a function */
  /* set i, pairs */
  int lengns;
  struct genename *gns;
  int ngenes;
  int li;
  int ei;
  char name[MAXLENGENENAME];

  for (li = 0; li < nloci; li++)
    {
      gns = NULL;
      lengns = 0;
      ngenes = L[li].numgenes;
      for (ei = 0; ei < ngenes; ei++)
        {
          strcpy (name, L[li].gNames[ei]);
          lengns = addpair (&gns, lengns, name, li, ei);
        }
      freepair (gns, lengns);
      gns = NULL;
    }

  for (li = 0; li < nloci; li++)
    {
      memcpy (L[li].pairsrelabel, L[li].pairs, sizeof (int) * MAXGENES);
    }
  /* Copy pairsrelabel from pairs: pairs had been reset to -1 when we ignore
   * diploidy of data 
   *   for (li = 0; li < nloci; li++)
   *     {
   *       memset (L[li].pairs, -1, sizeof (int) * MAXGENES);
   *     }
   */

  return;
}

void
IMA_genealogy_assignpopulation (int ci)
{
  int ii;
  int ij;
  int nei;
  struct genealogy *G; 
  struct edge *gtree; 
  int ei;
  int j;
  int k;
  int li;
  int pi;
  int *m;

  G = NULL;
  gtree = NULL;
  if (assignmentoptions[POPULATIONASSIGNMENT] == 0)
    {
      for (li = 0; li < nloci; li++)
        {
          G = &(C[ci]->G[li]);
          gtree = G->gtree;
          j = 0;
          k = 0;
          for (pi = 0; pi < npops; pi++)
            {
              while (j < L[li].samppop[pi] + k)
                {
                  gtree[j].pop = pi;
                  j++;
                }
              k += L[li].samppop[pi];
            }
        }
    }
  else if (assignmentoptions[POPULATIONASSIGNMENT] == 1)
    {
      m = malloc (snind * sizeof (int));
      for (ii = 0; ii < snind; ii++)
        {
          m[ii] = -2;
        }

      /* We find which individuals are known and which are not. */
      for (li = 0; li < nloci; li++)
        {
          G = &(C[ci]->G[li]);
          gtree = G->gtree;
          j = 0;
          k = 0;
          for (pi = 0; pi < npops; pi++)
            {
              while (j < L[li].samppop[pi] + k)
                {
                  ii = gtree[j].i;
                  if (m[ii] == -2)
                    {
                      m[ii] = pi;
                    }
                  else if (m[ii] == -1)
                    {
                      IM_err (IMERR_INPUTFILEINVALID, 
                              "Some genes of individual (ID: %d) have unknown label",
                              ii); 
                    }
                  else if (m[ii] != pi)
                    {
                      IM_err (IMERR_INPUTFILEINVALID, 
                              "Some genes of individual (ID: %d) are assigned to other label",
                              ii); 
                    }
                  j++;
                }
              k += L[li].samppop[pi];
            }
          while (j < L[li].numgenesunknown + k)
            {
              ii = gtree[j].i;
              if (m[ii] == -2)
                {
                  m[ii] = -1;
                }
              else if (m[ii] != -1)
                {
                  IM_err (IMERR_INPUTFILEINVALID, 
                          "Some genes of individual (ID: %d) have a known label",
                          ii);
                }
              j++;
            }
          if (j != L[li].numgenes)
            {
              IM_err (IMERR_INPUTFILEINVALID, "Invalid initial assignment"); 
            }
        }

      /* Assign individuals of unknown label to a radom population */
      snind_unknown = 0;
      for (ii = 0; ii < snind; ii++)
        {
          assert (m[ii] > -2);
          if (m[ii] == -1)
            {
              snind_unknown++; 
              m[ii] = randposint (npops);
            }
        }
 
      for (li = 0; li < nloci; li++)
        {
          G = &(C[ci]->G[li]);
          gtree = G->gtree;
          for (ei = 0; ei < L[li].numgenes; ei++)
            {
              ii = gtree[ei].i;
              assert (!(m[ii] < 0) && m[ii] < npops);
              gtree[ei].pop = m[ii];
            }
        }

      /* Check if all genes of an individual are assigned to the same label. */
      assert (snind > 0);
      for (ii = 0; ii < snind; ii++)
        {
          nei = sngenes_ind[ii];
          for (ij = 0; ij < nei; ij++)
            {
              li = saC.indlist[ii][ij].li;
              ei = saC.indlist[ii][ij].ei;
              G = &C[ci]->G[li];
              gtree = G->gtree;
              if (gtree[ei].i != ii)
                {
                  IM_err (IMERR_INPUTFILEINVALID, 
                          "Individual ID (%d) and its gene's individual ID (%d) do not match at locus %s %d-th gene",
                          ii, gtree[ei].i, L[li].name, ei + 1); 
                }
              if (gtree[ei].pop != m[ii])
                {
                  IM_err (IMERR_INPUTFILEINVALID, 
                          "The %d-th Gene of Individual (%d) at locus (%s) is assigned to a different population from other genes of the individual",
                          ei + 1, ii, L[li].name); 
                }
            }
        }
      free (m);
      m = NULL;
    }

  return;
}

void 
recordassignment_header (char *fname)
{
  int ci;
  int li;
  int ii;
  int ei;
  int li2;
  char strName[FNSIZE];
  char *assignfilename;
  int lenname;
  FILE *assignfile = NULL;
  FILE *namefile = NULL;
  FILE *fpassigngenealogy = NULL;

  /* we find individuals for each locus. */
  lenname = 10 + strlen (fname);
  assignfilename = malloc (lenname * sizeof (char));

  strcpy (assignfilename, fname);
  strcat (assignfilename, ".in.p");
  assignfile = fopen (assignfilename, "w");

  strcpy (assignfilename, fname);
  strcat (assignfilename, ".tex");
  namefile = fopen (assignfilename, "w");

  strcpy (assignfilename, fname);
  strcat (assignfilename, ".asn");
  fpassigngenealogy = fopen (assignfilename, "w");

  if (assignfile == NULL || namefile == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL, "Assignment output file cannot be created: %s", assignfilename);
    }

  ci = 0;
  fprintf (assignfile, "Gen");
  fprintf (assignfile, "\tlnL");
  fprintf (fpassigngenealogy, "Gen");
  fprintf (fpassigngenealogy, "\tlnL");

  fprintf (namefile, "\\begin{tabular}{l}\n");
  for (ii = 0; ii < snind; ii++)
    {
      li = 0;
      ei = sind2gi[ii][li];
      if (!(ei < 0))
        {
          removeSpaces (strName, L[li].gNames[ei]);
          fprintf (assignfile, "\tP(%s)", strName);
          fprintf (namefile, "%s\\tabularnewline\n", strName);
          fprintf (fpassigngenealogy, "\tP(%s)", strName);
        }
      else
        {
          /* If the individual does not have genes at first locus, we look for
           * it from the rest of loci. */
          for (li2 = 1; li2 < nloci; li2++)
            {
              ei = sind2gi[ii][li2];
              if (!(ei < 0))
                {
                  removeSpaces (strName, L[li2].gNames[ei]);
                  fprintf (assignfile, "\tP(%s)", strName);
                  fprintf (namefile, "%s\\tabularnewline\n", strName);
                  fprintf (fpassigngenealogy, "\tP(%s)", strName);
                  break;
                }
            }
        }
    }
  fprintf (namefile, "\\end{tabular}\n");

  fprintf (assignfile, "\n");
  fprintf (fpassigngenealogy, "\n");

  fclose (assignfile);
  assignfile = NULL;
  fclose (namefile);
  namefile = NULL;
  fclose (fpassigngenealogy);
  fpassigngenealogy = NULL;
  free (assignfilename);
  assignfilename = NULL;
  return;
}

int 
IMA_ninds ()
{
  assert (snind > 0);
  return snind;
}

int 
IMA_nindsunknown ()
{
  assert (!(snind_unknown < 0));
  return snind_unknown;
}

void 
recordassignment (char *fname, int ci, int gen)
{
  int li;
  int ii;
  int ei;
  FILE *assignfile = NULL;
  FILE *fpassigngenealogy = NULL;
  char assignfilename[FNSIZE];
  double pdg;
  int *a = NULL;
  int *pa;

  a = malloc (snind * sizeof (int));

  strcpy (assignfilename, fname);
  strcat (assignfilename, ".in.p");
  assignfile = fopen (assignfilename, "a");
  if (assignfile == NULL)
    { 
      IM_err (IMERR_APPENDFILEFAIL, "Assignment output file cannot be opened: %s", assignfilename);
    }
  fprintf (assignfile, "%d", gen);

  strcpy (assignfilename, fname);
  strcat (assignfilename, ".asn");
  fpassigngenealogy = fopen (assignfilename, "a");
  if (fpassigngenealogy == NULL)
    { 
      IM_err (IMERR_APPENDFILEFAIL, "Raw assignment output file cannot be opened: %s", assignfilename);
    }
  fprintf (fpassigngenealogy, "%d", gen);

  pdg = C[ci]->allpcalc.pdg;
  fprintf (assignfile, "\t%lf", pdg);
  fprintf (fpassigngenealogy, "\t%lf", pdg);

  /* how do I know individuals of a locus? */
  pa = a;
  for (ii = 0; ii < snind; ii++)
    {
      li = saC.indlist[ii][0].li;
      ei = saC.indlist[ii][0].ei;
      *pa = C[ci]->G[li].gtree[ei].pop;
      pa++;
    }

  /* A raw assignment is printed out without RGF conversion. */
  pa = a;
  for (ii = 0; ii < snind; ii++)
    {
      fprintf (fpassigngenealogy, "\t%d", *pa);
      pa++;
    }
  fprintf (fpassigngenealogy, "\n");

  /* An assignment is printed out with RGF conversion. */
  IMA_rgs_convert (a, snind);
  pa = a;
  for (ii = 0; ii < snind; ii++)
    {
      fprintf (assignfile, "\t%d", *pa);
      pa++;
    }
  fprintf (assignfile, "\n");

  fclose (assignfile);
  assignfile = NULL;
  fclose (fpassigngenealogy);
  fpassigngenealogy = NULL;
  free (a);
  a = NULL;
 
  return;
}

void 
recordassignmentloci (char *fname, int ci)
{
  static int gen = 0;

  recordassignment (fname, ci, gen);
  gen++;
  return;
}

void
start_setup_A ()
{
  /* For possible setup for Assignment */
}

/* initialized before start */
void
init_update_assignment (void)
{
  sbfold.e1 = NULL;
  sbfold.e2 = NULL;
  sbfold.e3 = NULL;
  sbfold.n1 = 0;
  sbfold.n2 = 0;
  sbfold.n3 = 0;
  sbfold.nmax1 = 0;
  sbfold.nmax2 = 0;
  sbfold.nmax3 = 0;
  sbfold.n = NULL;
  sbfold.l = NULL;
  sbfold.seqz1 = NULL;

  /* sbfold.seqz2 = NULL; */

  sbfnew.e1 = NULL;
  sbfnew.e2 = NULL;
  sbfnew.e3 = NULL;
  sbfnew.n1 = 0;
  sbfnew.n2 = 0;
  sbfnew.n3 = 0;
  sbfnew.nmax1 = 0;
  sbfnew.nmax2 = 0;
  sbfnew.nmax3 = 0;
  sbfnew.n = NULL;
  sbfnew.l = NULL;
  sbfnew.seqz1 = NULL;
  /* sbfnew.seqz2 = NULL; */

  snsample_pop = NULL;
  saGiold = NULL;
  saGioldz1 = NULL;
  saGioldz2 = NULL;
  saGinew = NULL;
  saGinewz1 = NULL;
  saGinewz2 = NULL;

  snasn = 0;
  return;
}

void
free_update_assignment (void)
{
  int li;
  int i;

  free_genealogy_weights (&saC.savedallgweight);
  free_probcalc (&saC.savedallpcalc);
  for (li = 0; li < nloci; li++)
    {
      free_genealogy_weights (&saT[li].savedgweight);
    }

  IMA_free_edgemiginfo ();
  IMA_genealogy_save_free ();
  IMA_ind_fin ();

  IMA_unset_popntree();

  if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
    {
      assert (nloci == 1);

      for (i = 0; i < L[0].numgenesunknown; i++)
        {
          free (sassigned[i]);
          sassigned[i] = NULL;
        }
      free (sassigned);
      sassigned = NULL;
    }

  IMA_delete_intervals ();
  IMA_memory_freesaGi ();
  IMA_memory_freennminus1 ();

#ifdef COUNT_PRIOR_ASSIGNMENT
  free (saasn);
  saasn = NULL;
#endif /* COUNT_PRIOR_ASSIGNMENT */

  return;
}

void
IMA_initmemory_aftersetpopntree ()
{
  int i;
  int n;

  if (modeloptions[NOMIGRATION] == 0)
    {
      n = IMA_number_largest_branches ();
    }
  IMA_init_edgemiginfo ();
  IMA_genealogy_save_init ();

  IMA_ind_init ();
  IMA_ind_set ();
  /* called after setup_popntree */
  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
    {
      IMA_set_popntree_nosplit ();
    }
  else
    {
      IMA_set_popntree ();
    }

  if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
    {
      assert (nloci == 1);
      sassigned = malloc (L[0].numgenesunknown * sizeof (double *));
      for (i = 0; i < L[0].numgenesunknown; i++)
        {
          sassigned[i] = malloc (npops * sizeof (double));
          memset (sassigned[i], 0, npops * sizeof (double));
        }
    }
  
  IMA_create_intervals ();
  IMA_memory_initsaGi ();
  IMA_memory_initnnminus1 ();

#ifdef COUNT_PRIOR_ASSIGNMENT
  snasn = IMA_rgf_setSnasn (snind, npops);
  assert (snasn > 0);
  saasn = malloc (snasn * sizeof (int));
  memset (saasn, 0, snasn * sizeof (int));
#endif /* COUNT_PRIOR_ASSIGNMENT */

  return;
}

void
IMA_memory_initnnminus1 ()
{
  int i;
  int li;
  int largestsamp;
  for (i = 0, li = 0; li < nloci; li++)
    {
      if (i < L[li].numgenes)
        {
          i = L[li].numgenes;
        }
    }
  largestsamp = i;
  snnminus1 = malloc ((largestsamp + 1) * sizeof (double));
  snnminus1[0] = 0;
  for (i = 1; i <= largestsamp; i++)
    {
      snnminus1[i] = (double) (i) * ((double) i - 1);
    }
  return;
}

void 
IMA_memory_freennminus1 ()
{
  free (snnminus1);
  snnminus1 = NULL;
}

int
IMA_number_largest_branches ()
{
  int v;
  int li;

  v = 0;
  for (li = 0; li < nloci; li++)
    {
      if (v < L[li].numlines)
        {
          v = L[li].numlines;
        }
    }
  return v;
}
void
IMA_assigned_increase ()
{
  int i;
  int ei;
  int li;
  struct genealogy *G;
  struct edge *gtree;

  assert (nloci == 1);
  for (li = 0; li < nloci; li++)
    {
      G = &(C[0]->G[li]);
      gtree = G->gtree;

      for (i = 0; i < L[li].numgenesunknown; i++)
        {
          ei =  L[li].numgenesknown + i;
          sassigned[i][gtree[ei].pop] += 1.0;
        }
    }
  return;
}

void
IMA_assigned_reset ()
{
  int i;
  int j;
  int li;

  assert (nloci == 1);
  for (li = 0; li < nloci; li++)
    {
      for (i = 0; i < L[li].numgenesunknown; i++)
        {
          for (j = 0; j < npops; j++)
            {
              sassigned[i][j] = 0.0;
            }
        }
    }
  return;
}

void
print_num_genes_pops_stdout ()
{
  print_num_genes_pops (stdout);
}

void
print_num_genes_pops (FILE *fp)
{
  int ii;
  int ei;
  int i;
  int j;
  int li;
  int ngenes[MAXPOPS];
  double sumpops;
  int ci = 0;
  struct edge *gtree = NULL;
  struct genealogy *G = NULL;

  for (i = 0; i < npops; i++)
    {
      ngenes[i] = 0;
    }

  for (ii = 0; ii < snind; ii++)
    {
      li = saC.indlist[ii][0].li;
      ei = saC.indlist[ii][0].ei;
      ngenes[C[0]->G[li].gtree[ei].pop]++;
    } 
     
  fprintf (fp, "Number of genes assigned to populations: ");
  for (i = 0; i < npops; i++)
    {
      fprintf (fp, "[%d] %d\t", i, ngenes[i]);
    }
  fprintf (fp, "\n");
  if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
    {
      assert (nloci == 1);
      li = 0;
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      fprintf (fp, "Locus [%d]: ", li);
      for (i = L[li].numgenesknown; i < L[li].numgenes; i++)
        {
          fprintf (fp, "gene[%d] -> [%d], ", i, gtree[i].pop);
        }
      fprintf (fp, "\n\n");
      fprintf (fp, "Percentage\n");
      fprintf (fp, "----------\n");

      li = 0;
      for (i = 0; i < L[li].numgenesunknown; i++)
        {
          ei = L[li].numgenesknown + i;
          fprintf (fp, "gene[%d]: ", ei);

          sumpops = 0.0;
          for (j = 0; j < npops; j++)
            {
              sumpops += sassigned[i][j];
            }
          for (j = 0; j < npops; j++)
            {
              fprintf (fp, "[%d] %lf, ", j, sassigned[i][j] / sumpops);
            }
        }
      fprintf (fp, "\n");
    }

  return;
}
void
removeSpaces (char *s1, const char *s2)
{
  int j;
  int i;
  int len;

  len = strlen (s2);
  j = 0;
  for (i = 0; i < len; i++)
    {
      if (isalnum(s2[i]) 
          || s2[i] == '_' 
          || s2[i] == '-'
          || s2[i] == '@')
        {
          s1[j] = s2[i];
          j++;
        }
    }
  s1[j] = '\0';
  return;
}

int
IMA_genealogy_sister (int ci, int li, int ei)
{
  int v;
  int gparent;

  v = -1;
  gparent = C[ci]->G[li].gtree[ei].down;
  if (!(gparent < 0))
    {
      if (C[ci]->G[li].gtree[gparent].up[0] == ei)
        {
          v = C[ci]->G[li].gtree[gparent].up[1];
        }
      else
        {
          v = C[ci]->G[li].gtree[gparent].up[0];
        }
    }
  return v;
}

double
IMA_genealogy_time (int ci, int li, int ei)
{
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;
  int gchild;
  double t;

  G = &(C[ci]->G[li]);
  gtree = G->gtree;

  gchild = gtree[ei].up[0]; /* this may be also up[1]. */
  if (gchild < 0)
    {
      t = 0.0;
    }
  else
    {
      assert (gchild < L[li].numlines);
      t = gtree[gchild].time;
    }
  return t;
}
int
findperiodoftime (int ci, double ktime)
{
  int ti; 

  ti = 0;
  while (ktime > C[ci]->tvals[ti])
    ti++;

  return ti;

}

int
IMA_edge_labeloftime (int ci, int li, int ei, double ktime)
{
  int i;
  int nowpop;
  int kperiod;
  struct genealogy *G = NULL;    /* locus li */
  struct edge *gtree = NULL; /* genealogy of locus li */

  G = &(C[ci]->G[li]);
  gtree = C[ci]->G[li].gtree;
  assert (IMA_genealogy_time (ci, li, ei) < ktime);
  assert (gtree[ei].time > ktime);
  i = 0;
  while (gtree[ei].mig[i].mt > -0.5 && gtree[ei].mig[i].mt < ktime)
    {
      i++;
    }
  if (i > 0)
    {
      nowpop = gtree[ei].mig[i - 1].mp;
    }
  else
    {
      nowpop = gtree[ei].pop;
    }
  assert (!(nowpop < 0));
  kperiod = findperiodoftime (ci, ktime);
  nowpop = saC.popndown[nowpop][kperiod];
  assert (!(nowpop < 0));
  assert (nowpop < 2 * npops - 1);
  return nowpop;
}

int 
IMA_genealogy_nmigration (int ci, int li, int ei)
{
  int n;
  int gparent;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;
  
  assert (!(ci < 0));
  assert (!(li < 0));
  
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  gparent = gtree[ei].down;

  if (gparent < 0)
    {
      if (gtree[ei].cmm > 0)
        {
          assert (gtree[ei].mig[0].mt < -0.5);
        }
    }

  n = 0;

  if (gtree[ei].cmm > 0)
    {
      while (gtree[ei].mig[n].mt > -0.5)
        {
          n++; 
        }
      assert (!(n < 0));
    }
  
  gtree->nmig = n;
  return n;
}

int 
IMA_genealogy_nmigrationupto (int ci, int li, int ei, double t)
{
  int n;
  int gparent;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;
  
  assert (!(ci < 0));
  assert (!(li < 0));
  assert (t > 0.0);
  
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  gparent = gtree[ei].down;

  if (gparent < 0)
    {
      if (gtree[ei].cmm > 0)
        {
          assert (gtree[ei].mig[0].mt < -0.5);
        }
    }

  n = 0;

  if (gtree[ei].cmm > 0)
    {
      while (gtree[ei].mig[n].mt > -0.5  && gtree[ei].mig[n].mt < t)
        {
          n++; 
        }
      assert (!(n < 0));
    }
  
  return n;
}
 
/*
 * Print genealogical trees in Newick formated strings
 */
void 
IMA_genenode_print_node (FILE *fp, struct edge *gtree, int e)
{
  int mi;

  fprintf (fp, "%d[%d]", e, gtree[e].pop);

  if (modeloptions[NOMIGRATION] == 0)
    {
                for (mi = 0; mi < gtree[e].cmm; mi++)
                  {
                         if (gtree[e].mig[mi].mt < -0.5)
                                {
                                  break;
                                }

          fprintf (fp, "-[%d@%.1lf]", gtree[e].mig[mi].mp, gtree[e].mig[mi].mt);
                  }
    }

  return;
}

void 
IMA_genenode_print (FILE *fp, struct edge *gtree, int e)
{
  int left;
  int right;
  int parent;

  left = gtree[e].up[0];
  right = gtree[e].up[1];
  parent = gtree[e].down;

  assert ((left < 0 && right < 0) || !(left < 0 && right < 0));

  if (left < 0 && right < 0)
    {
      IMA_genenode_print_node (fp, gtree, e);
    }
  else
    {
      fprintf (fp, "(");
      IMA_genenode_print (fp, gtree, left);
      fprintf (fp, ",");
      IMA_genenode_print (fp, gtree, right);
      fprintf (fp, ")");
      IMA_genenode_print_node (fp, gtree, e);
    }
  if (parent < 0)
    {
      fprintf (fp, ";\n");
    }
  else
    {
      fprintf (fp, ":%lf", gtree[e].time - gtree[left].time);
    }
  return;
}

void 
IMA_genenode_print_nodeseq (FILE *fp, struct edge *gtree, int e, int si)
{
  fprintf (fp, "%d[%d]", e, gtree[e].seq[si]);

  return;
}


void 
IMA_genenode_printseq (FILE *fp, struct edge *gtree, int e, int si)
{
  int left;
  int right;
  int parent;

  left = gtree[e].up[0];
  right = gtree[e].up[1];
  parent = gtree[e].down;

  assert ((left < 0 && right < 0) || !(left < 0 && right < 0));

  if (left < 0 && right < 0)
    {
      IMA_genenode_print_nodeseq (fp, gtree, e, si);
    }
  else
    {
      fprintf (fp, "(");
      IMA_genenode_printseq (fp, gtree, left, si);
      fprintf (fp, ",");
      IMA_genenode_printseq (fp, gtree, right, si);
      fprintf (fp, ")");
      IMA_genenode_print_nodeseq (fp, gtree, e, si);
    }
  if (parent < 0)
    {
      fprintf (fp, ";\n");
    }
  else
    {
      fprintf (fp, ":%lf", gtree[e].time - gtree[left].time);
    }
  return;
}

void 
x()
{
  IMA_genealogy_printall (0);
  return;
}

char* 
IMA_genealogy_print (int ci, int li)
{
  char *ns = NULL;
  FILE *fp;
  struct genealogy *G = NULL;    /* locus li */
  struct edge *gtree = NULL; /* genealogy of locus li */
  
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  fp = fopen ("1.tre", "a");
  IMA_genenode_print (fp, gtree, G->root);
  fclose (fp);
  fp = NULL;
  
  return ns;
}

char* 
IMA_genealogy_printall (int ci)
{
  char *ns = NULL;
  FILE *fp;
  struct genealogy *G = NULL;    /* locus li */
  struct edge *gtree = NULL; /* genealogy of locus li */
  int li;
  
  fp = fopen ("1.tre", "w");
  for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      IMA_genenode_print (fp, gtree, G->root);
    }
  fclose (fp);
  fp = NULL;

  return ns;
}

char* 
IMA_genealogy_printsite (int ci, int li, int si)
{
  char *ns = NULL;
  FILE *fp;
  struct genealogy *G = NULL;    /* locus li */
  struct edge *gtree = NULL; /* genealogy of locus li */
  
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  fp = fopen ("1.tre", "w");
  IMA_genenode_printseq (fp, gtree, G->root, si);
  fclose (fp);
  fp = NULL;
  
  return ns;
}

char* 
IMA_genealogy_printseq (int ci, int li)
{
  int si;
  char *ns = NULL;
  FILE *fp;
  struct genealogy *G = NULL;    /* locus li */
  struct edge *gtree = NULL; /* genealogy of locus li */
  
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  fp = fopen ("1.tre", "w");
  for (si = 0; si < L[li].numsites; si++)
    {
      IMA_genenode_printseq (fp, gtree, G->root, si);
    }
  fclose (fp);
  fp = NULL;
  
  return ns;
}

void
x1 (im_bfupdate *pbf)
/* IMA_sbf_print (im_bfupdate *pbf) */
{
  int i;
  int j;

  printf ("---------\n");
  printf ("Event e1 [%d]:\n", pbf->n1);
  for (i = 0; i < pbf->n1; i++)
    {
      printf ("[%2d] %c - [p:%2d, pi:%2d, pj:%2d] edge [%2d] at %lf\n", 
              i, 
              pbf->e1[i].type,  
              pbf->e1[i].p,
              pbf->e1[i].pi,
              pbf->e1[i].pj,
              pbf->e1[i].ei,
              pbf->e1[i].t);
    }
  printf ("\n"); 

  printf ("---------\n");
  printf ("Event e2 [%d]:\n", pbf->n2);
  for (i = 0; i < pbf->n2; i++)
    {
      printf ("[%2d] %c - [p:%2d, pi:%2d, pj:%2d] edge [%2d] at %lf\n", 
              i, 
              pbf->e2[i].type,  
              pbf->e2[i].p,
              pbf->e2[i].pi,
              pbf->e2[i].pj,
              pbf->e2[i].ei,
              pbf->e2[i].t);
    }
  printf ("\n"); 
 
  printf ("---------\n");
  printf ("Event e3 [%d]:\n", pbf->n3);
  for (i = 0; i < pbf->n3; i++)
    {
      printf ("[%2d] %c - [p:%2d, pi:%2d, pj:%2d] edge [%2d] at %lf\n", 
              i, 
              pbf->e3[i].type,  
              pbf->e3[i].p,
              pbf->e3[i].pi,
              pbf->e3[i].pj,
              pbf->e3[i].ei,
              pbf->e3[i].t);
    }
  printf ("\n"); 

  printf ("---------\n");
  if (pbf->m1 != NULL)
    {
      printf ("z1: edge [%2d], sis [%2d]\n",
              pbf->m1->edgeid, pbf->sis1);
/* 
      printf ("z1: edge [%2d], sis [%2d], down [%2d], downdown [%2d]\n",
              pbf->m1->edgeid, pbf->sis1, pbf->down1, pbf->downdown1);
      printf ("    upt [%lf], sisupt [%lf], sisdnt [%lf]\n",
              pbf->upt1, pbf->sisupt1, pbf->sisdnt1);
*/
      IMA_edgemiginfo_print (pbf->m1);
      printf ("\n");
    }
  else
    {
      printf ("z1:\n");
    }

/* DIPLOID
  printf ("---------\n");
  if (pbf->m2 != NULL)
    {
      printf ("z2: edge [%2d], sis [%2d], down [%2d], downdown [%2d]\n",
              pbf->m2->edgeid, pbf->sis2, pbf->down2, pbf->downdown2);
      printf ("    upt [%lf], sisupt [%lf], sisdnt [%lf]\n",
              pbf->upt2, pbf->sisupt2, pbf->sisdnt2);
      IMA_edgemiginfo_print (pbf->m2);
      printf ("\n");
    }
  else
    {
      printf ("z2:\n");
    }
*/

  printf ("---------\n");
  if (pbf->m3 != NULL)
    {
      printf ("z3:\n");
      IMA_edgemiginfo_print (pbf->m3);
      printf ("\n");
    }
  else
    {
      printf ("z3:\n");
    }

  printf ("---------\n");
  printf ("lz1:\n");
  for (i = 0; i < numpopsizeparams; i++)
    {
      printf ("[%d] %2d: ", i, pbf->nz1[i]);
      for (j = 0; j < pbf->nz1[i]; j++)
        {
          printf (" %2d", pbf->lz1[i][j]);
        }
      printf ("\n");
    }

/* DIPLOID
  printf ("---------\n");
  printf ("lz2:\n");
  for (i = 0; i < numpopsizeparams; i++)
    {
      printf ("[%d] %2d:", i, pbf->nz2[i]);
      for (j = 0; j < pbf->nz2[i]; j++)
        {
          printf (" %2d", pbf->lz2[i][j]);
        }
      printf ("\n");
    }
*/

  printf ("---------\n");
  printf ("l:\n");
  for (i = 0; i < numpopsizeparams; i++)
    {
      printf ("[%d]: %2d\n", i, pbf->n[i]);
    }
   
  return;
}

void
IMA_edgemiginfo_print (im_edgemiginfo *m)
{
  int mi;

  if (m == NULL)
    {
      return;
    }

  printf ("    miginfo: edge [%2d] from %lf to %lf\n", m->edgeid, m->upt, m->dnt);
  printf ("             popn [%2d] to [%2d]\n", m->pop, m->fpop);
  mi = 0;
  while (m->mig[mi].mt > -0.5)
    {
      printf ("    mp [%2d] at %lf\n", m->mig[mi].mp, m->mig[mi].mt);
      mi++;
    }
  printf ("    mi vs. mpall: %2d == %2d\n", mi, m->mpall);
  return;
}

int
x2 (int old, int li)
{  
  im_ginfo *gk;

  if (old == 0)
    {
      printf ("saGiold\n");
      printf ("-------\n");
      gk = &saGiold[li];
      x2print (gk);
      printf ("saGioldz1\n");
      printf ("-------\n");
      gk = &saGioldz1[li];
      x2print (gk);
    }
  else
    {
      printf ("saGinew\n");
      printf ("-------\n");
      gk = &saGinew[li];
      x2print (gk);
      printf ("saGinewz1\n");
      printf ("-------\n");
      gk = &saGinewz1[li];
      x2print (gk);
    }
  return 0;
}

int
x2print (im_ginfo *gk)
{  
  
  int pi;
  int pj;
  int ti;

  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      printf ("[%d]: cc (%3d) fc (%lf) - theta (%lf)\n", 
              pi, gk->cc[pi], gk->fc[pi], gk->thetas[pi]);
    }

  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      for (pj = 0; pj < numpopsizeparams; pj++)
        {
          if (pi == numpopsizeparams - 1
              || pj == numpopsizeparams - 1)
            {
              continue;
            }
          if (pi == pj)
            {
              continue;
            }
          printf ("[%d -> %d]: cm (%3d) fm (%lf) - mms (%lf)\n", 
                  pi, pj, gk->cm[pi][pj], gk->fm[pi][pj], gk->mms[pi][pj]);
        }
    }

  /* We need number of periods where migrations are allowed. */
  /* npops - 1 for tree model and 1 for tree model.          */
  /* lastperiodnumber is for this                            */
  for (ti = 0; ti < lastperiodnumber; ti++)
    {
      for (pi = 0; pi < numpopsizeparams; pi++)
        {
          if (saC.popnmig[ti][pi] != NULL)
            {
              printf ("[%d, %d], ms (%lf)\n", ti, pi, gk->ms[ti][pi]);
            }
        }
    }
  return 0;
}

int
IMA_node_period (int ci, int li, int ei)
{
  int gchild;
  double t;
  int ti; 
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;

  G = &(C[ci]->G[li]);
  gtree = G->gtree;

  gchild = gtree[ei].up[0]; /* this may be also up[1]. */
  if (gchild < 0)
    {
      t = 0.0;
    }
  else
    {
      assert (gchild < L[li].numlines);
      t = gtree[gchild].time;
    }

  ti = 0;
  while (t > C[ci]->tvals[ti])
    {
      ti++;
    }

  return ti;
}

int
IMAedge_downperiod (int ci, int li, int ei)
{
  double t;
  int ti; 
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;

  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  t = gtree[ei].time;

  ti = 0;
  while (t > C[ci]->tvals[ti])
    {
      ti++;
    }

  return ti;
}

double
IMAedge_length (int ci, int li, int ei)
{
  double len;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
    
  assert (ei != G->root);
  len = gtree[ei].time - IMA_genealogy_time (ci, li, ei);
  return len;
}

int
IMAmig_period (int ci, int li, int ei, int mi)
{
  double t;
  int ti; 
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;

  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  t = gtree[ei].mig[mi].mt;
  assert (!(t < 0.0));

  ti = 0;
  while (t > C[ci]->tvals[ti])
    {
      ti++;
    }

  return ti;
}

void 
BitPrintF (UByteP A, int size)
{
  int i;
  fprintf (stdout, "Set: ");
  for (i = 0; i < size*8; i++) 
    { 
      if (BitIsTrue (A, i)) 
        { 
          fprintf (stdout, "%d ", i);
        } 
    }
  fprintf (stdout, "\n");
  return;
}

void 
BitSPrintF (char *s, UByteP A, int size)
{
  int i;
  strcpy (s, "#");
  for (i = 0; i < size*8; i++) 
    { 
      if (BitIsTrue (A, i)) 
        { 
                    sprintf (s, "%s&%d", s, i); 
        } 
    }
  sprintf (s, "%s", s);
  return;
}
void
print_popnlist ()
{
  int i;
  int j;
  unsigned int k;
  int npnodes = 2 * npops - 1;
  assert (!(npops < 0));
  assert (npops < MAXPOPS);
  for (i = 0; i < npnodes; i++)
    {
      printf ("%d:(", i);
      for (j = 0; j < npops; j++)
        {
          if (saC.popnlist[i][j] != NULL)
            {
              printf ("P%d:", j);
              for (k = 0; k < saC.popnlist[i][j]->size; k++)
                {
                  printf ("%d,", saC.popnlist[i][j]->data[k]);
                }
            }
        }
      printf (")\n");
    }
  return;
}

double 
IMA_assignment2value (int ci, int li)
{
  double v;  
  v = IMA_assignment2value_mean (ci, li);
  return v;
}
double
IMA_assignment2value_mean (int ci, int li)
{
  double v;
  int ei;
  int ngenes;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;

  v = 0.0;
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  ngenes = L[li].numgenes;
  for (ei = 0; ei < ngenes; ei++)
    {
      v += (double) gtree[ei].pop;
    }
  v /= (double) ngenes; 
  G->asn = v;
  return v; 
}


int
skip_a_line (FILE *fp)
{
  char c;
  c = 'x';
  while (c != '\n')
    {
      c = fgetc(fp);
    }
  return 0;
}

int
read_int (FILE *fp)
{
  int v;
  char c;
  char *buf;
  int len;
  int i;
  
  len = 5;
  buf = malloc (len * sizeof (char));
  
  c = fgetc (fp);
  i = 0;
  while (c != '\t' && c != '\n')
    {
      buf[i++] = c;
      if (i == len)
        {
          len += 5;
          buf = realloc (buf, len * sizeof (char));
        }
      c = fgetc (fp);
    }
  buf[i] = '\0';
  v = atoi (buf);
  
  free (buf);
  buf = NULL;
  return v;
}

double 
read_double (FILE *fp)
{
  double v;
  char c;
  char *buf;
  int len;
  int i;
  
  len = 5;
  buf = malloc (len * sizeof (char));
  
  c = fgetc (fp);
  i = 0;
  while (c != '\t' && c != '\n')
    {
      buf[i++] = c;
      if (i == len)
        {
          len += 5;
          buf = realloc (buf, len * sizeof (char));
        }
      c = fgetc (fp);
    }
  buf[i] = '\0';
  v = atof (buf);
  
  free (buf);
  buf = NULL;
  
  return v;
}

int
updateassignmentrelabel (int ci)
{
  int accp;                   /* accept (1) or reject (0)            */
  double logweight;           /* MH ratio                            */
  double logU;                /* log uniform                         */
  double tpw;                 /* log [p(G*)/p(G)]                    */
  double migweight;           /* hastings ratio for migration update */
  struct genealogy *G;        /* locus li                            */
  int mpart;
  int li;
  int nmigaftertreeweight;

  /* declarations should precede any statement for compatibility with Visual C++
   * compiler */
  G = NULL;
  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }

  if (ci == 0)
    {
      for (li = 0; li < nloci; li++)
        {
          L[li].a_rec->upinf[IM_UPDATE_ASSIGNMENT_RELABEL].tries++;
        }
    }
  Cupinf[ci].upinf[IM_UPDATE_ASSIGNMENT_RELABEL].tries++;

  /* save the current state of a chain */
  for (li = 0; li < nloci; li++)
    {
      IMA_genealogy_save_reset (li);
      IMA_savelocusinfo (ci, li);
      copy_treeinfo (&saT[li].savedgweight, &C[ci]->G[li].gweight);
    }
  copy_treeinfo (&saC.savedallgweight, &C[ci]->allgweight);
  copy_probcalc (&saC.savedallpcalc, &C[ci]->allpcalc);

  /* the key proposal function of relabel method */
  migweight = relabel (ci, &mpart);
  assert (modeloptions[NOMIGRATION] == 0);

  /* compute Metropolis-ratio */
  setzero_genealogy_weights (&C[ci]->allgweight);
  for (li = 0; li < nloci; li++)
    {
      setzero_genealogy_weights (&C[ci]->G[li].gweight);
      treeweight (ci, li);
      sum_treeinfo (&C[ci]->allgweight, &C[ci]->G[li].gweight);
    }
  /* check if the number of migration events after the call of treeweight is the
   * same as the number of migration events from update */
  nmigaftertreeweight = 0;
  for (li = 0; li < nloci; li++)
    {
       nmigaftertreeweight += C[ci]->G[li].mignum; 
    }
  assert (nmigaftertreeweight == mpart);
  /* We use allgweight to compute allpcalc */
  initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);

  /* Log[p(G*)] - Log[p(G)] is [[tpw]] */
  tpw = C[ci]->allpcalc.probg - saC.savedallpcalc.probg;

  logU = log(uniform ());
  logweight = beta[ci] * tpw + migweight;
  assert (beta[ci] * tpw + migweight > -1e200
          && beta[ci] * tpw + migweight < 1e200);

  accp = 0;
  if (logweight > 0.0 || logweight > logU)
    {
      /* accept the update */
      accp = 1;
      if (ci == 0)
        {
          for (li = 0; li < nloci; li++)
            {
              L[li].a_rec->upinf[IM_UPDATE_ASSIGNMENT_RELABEL].accp++;
            }
        }
      Cupinf[ci].upinf[IM_UPDATE_ASSIGNMENT_RELABEL].accp++;
    }
  else
    {
      for (li = 0; li < nloci; li++)
        {
          /* restore edges changed during update */
          IMA_genealogy_restore_edges (ci, li); 
          IMA_restorelocusinfo (ci, li);
          copy_treeinfo (&C[ci]->G[li].gweight, &saT[li].savedgweight);
        }
      copy_probcalc (&C[ci]->allpcalc, &saC.savedallpcalc);
      copy_treeinfo (&C[ci]->allgweight, &saC.savedallgweight);
    }

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
  return accp;
}

double
relabel (int ci, int *mpart)
{
  int li;
  int ei;
  int pop;
  double logprob;
  int newpop;
  int oldmigcount[MAXLOCI];
  int newmigcount[MAXLOCI];
  struct genealogy *G;
  struct edge *gtree; 
  int ii;
  int ij;
  int nei;

  assert (modeloptions[NOMIGRATION] == 0);
  if (assignmentoptions[POPULATIONASSIGNMENTASSIGNED] == 1)
    {
      assert (nloci == 1);
    }
  for (li = 0; li < nloci; li++)
    {
      IMA_genealogy_save_reset (li);
    }

  /* pick an individual to relabel */
  ii = -1;
  li = randposint (nloci);
  while (ii < 0)
    {
      if (L[li].numgenesunknown > 0)
        {
          ei = L[li].numgenesknown + randposint (L[li].numgenesunknown);
          ii = C[ci]->G[li].gtree[ei].i;
        }
      li = randposint (nloci);
    }
  assert (ii < snind);

  assert (!(ii < 0));
  assert (ii < snind);
  nei = sngenes_ind[ii];
  for (ij = 0; ij < nei; ij++)
    {
      li = saC.indlist[ii][ij].li;
      ei = saC.indlist[ii][ij].ei;
      IMA_genealogy_save_edge (ci, li, ei);
    }

  /* relabel the selected genes: newpop is chosen with logprob */
  /* we need a list of edge mig info structures */
  li = saC.indlist[ii][0].li;
  ei = saC.indlist[ii][0].ei;
  G = &C[ci]->G[li];
  gtree = G->gtree;
  pop = gtree[ei].pop;
  logprob = IMA_popn_move (&newpop, ci, pop, C[ci]->rootpop);
  assert (logprob == 0.0);

  for (ij = 0; ij < nei; ij++)
    {
      li = saC.indlist[ii][ij].li;
      ei = saC.indlist[ii][ij].ei;
      IMA_edge_fillmiginfo (ci, li, ei, &saoems[ij]);
    }

  /* not sure if I need this storeoldedges */
  /* storeoldedges (ci, li, edge, oldsis, freededge); */
  for (li = 0; li < nloci; li++)
    {
      oldmigcount[li] = C[ci]->G[li].mignum;
    }
  for (ij = 0; ij < nei; ij++)
    {
      li = saC.indlist[ii][ij].li;
      ei = saC.indlist[ii][ij].ei;
      C[ci]->G[li].gtree[ei].mig[0].mt = -1;
      logprob += addmigration_relabel (ci, li, newpop, 
                                       oldmigcount[li],
                                       &newmigcount[li],
                                       &saoems[ij], &sanems[ij]);
      oldmigcount[li] = newmigcount[li];
    }
 
  *mpart = 0;
  for (li = 0; li < nloci; li++)
    {
      *mpart += oldmigcount[li];
    }

  for (ij = 0; ij < nei; ij++)
    {
      li = saC.indlist[ii][ij].li;
      ei = saC.indlist[ii][ij].ei;
      assert (sanems[ij].edgeid == ei);
      assert (sanems[ij].li == li);
      IMA_edge_copynewmig (ci, li, &sanems[ij]);
    }

  for (ij = 0; ij < nei; ij++)
    {
      li = saC.indlist[ii][ij].li;
      ei = saC.indlist[ii][ij].ei;
      gtree = C[ci]->G[li].gtree;
      gtree[ei].pop = newpop;
    }

  return logprob;
}
void
IMA_edge_fillmiginfo (int ci, int li, int ei, 
                      struct edgemiginfo *em)
{
  int i;
  int j;
  int down;                        /* down edge             */
  struct genealogy *G;             /* locus li              */
  struct edge *gtree;              /* genealogy of locus li */

  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  down = gtree[ei].down;

  /*******************************************************/
  /* oldedgemig                                          */
  /*******************************************************/
  IMA_reset_edgemiginfo (em);
  em->li = li;
  em->sisid = -1; /* do we need sisid? */
  em->edgeid = ei;
  em->upt = IMA_genealogy_time (ci, li, ei);
  em->dnt = gtree[ei].time;
  em->pop = gtree[ei].pop;
  em->temppop = gtree[ei].pop;
  em->fpop = gtree[down].pop;
  /* b, e, mtimeavail, and mtall are set. */
  fillmiginfoperiods (ci, em);
  /* migration events are copied */
  em->mig[0].mt = -1.0;
  for (j = 0; j <= em->b; j++)
    {
      em->mp[j] = 0;
    }
  em->mpall = 0;

  i = 0;
  j = em->b;
  while (gtree[ei].mig[i].mt > -0.5)
    {
      while (gtree[ei].mig[i].mt > C[ci]->tvals[j])
        {
          j++;
          em->mp[j] = 0;
        }
      em->mp[j]++;
      em->mpall++;
      em->mig[i] = gtree[ei].mig[i];
      i++;
    }
  em->mig[i].mt = -1;

  return;
}
void
IMA_edge_copynewmig (int ci, int li, struct edgemiginfo *em)
{
  int i;
  int ei;
  struct edge *gtree;              /* genealogy of locus li */
  
  ei = em->edgeid;
  gtree = C[ci]->G[li].gtree;

  i = 0;
  while (em->mig[i].mt > 0)
    {
      checkmigt (i, &gtree[ei]);
      copymig (&(gtree[ei].mig[i]), &(em->mig[i]));
      i++;
    }
  gtree[ei].mig[i].mt = -1;
  return;
}

void
IMA_restore_edgemig (int ci, int li, int ei, im_migstruct *savedmig, int ncmmsavedmig)
{
  int mi;
  struct edge *gtree = NULL;

  gtree = C[ci]->G[li].gtree;
  for (mi = 0; mi < ncmmsavedmig; mi++)
    {
      gtree[ei].mig[mi].mt = savedmig[mi].mt;
      gtree[ei].mig[mi].mp = savedmig[mi].mp;
      if (savedmig[mi].mt < -0.5)
        {
          break;
        }
    }

  return;
}

int
updateassignmentrelabellocus (int ci, int li)
{
  double logweight;           /* MH ratio                            */
  double logU;                /* log uniform                         */
  double tpw;                 /* log [p(G*)/p(G)]                    */
  int accp;                   /* accept (1) or reject (0)            */
  double migweight;           /* hastings ratio for migration update */
  struct genealogy *G = NULL; /* locus li                            */
  int mpart;

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogyloc (ci, li);
    }

  if (ci == 0)
    {
      L[li].a_rec->upinf[IM_UPDATE_ASSIGNMENT_RELABEL].tries++;
    }
  Cupinf[ci].upinf[IM_UPDATE_ASSIGNMENT_RELABEL].tries++;

  G = &(C[ci]->G[li]);
  IMA_savelocusinfo (ci, li);
  copy_treeinfo (&saT[li].savedgweight, &G->gweight);
  copy_treeinfo (&saC.savedallgweight, &C[ci]->allgweight);
  copy_probcalc (&saC.savedallpcalc, &C[ci]->allpcalc);

  /* the key proposal function of relabel method */
  migweight = relabellocus (ci, li, &mpart);
  if (migweight == -DBL_MAX)
    {
      assert (modeloptions[NOMIGRATION] == 1);
      if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
        {
          assertgenealogyloc (ci, li);
        }
      return 0;
    }

  /* compute Metropolis ratio */
  setzero_genealogy_weights (&G->gweight);
  treeweight (ci, li);
  assert (G->mignum == mpart);
  sum_subtract_treeinfo (&C[ci]->allgweight, &G->gweight, &saT[li].savedgweight);
  initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);
  tpw = C[ci]->allpcalc.probg - saC.savedallpcalc.probg;
  logU = log(uniform ());
  logweight = beta[ci] * tpw + migweight;
  assert (beta[ci] * tpw + migweight > -1e299
          && beta[ci] * tpw + migweight < 1e299);

  accp = 0;
  if (logweight > 0.0 || logweight > logU)
    {
      /* accept the update */
      accp = 1;
      if (ci == 0)
        {
          L[li].a_rec->upinf[IM_UPDATE_ASSIGNMENT_RELABEL].accp++;
        }
      Cupinf[ci].upinf[IM_UPDATE_ASSIGNMENT_RELABEL].accp++;
    }
  else
    {
      /*****************************************************************/
      /* DELETE THIS?                                                  */
      /* NOTE: restoreedges (ci, li, edge, oldsis, freededge, newsis); */
      /*****************************************************************/
      IMA_restorelocusinfo (ci, li);
      IMA_genealogy_restore_edges (ci, li);
      copy_treeinfo (&G->gweight, &saT[li].savedgweight);
      copy_probcalc (&C[ci]->allpcalc, &saC.savedallpcalc);
      copy_treeinfo (&C[ci]->allgweight, &saC.savedallgweight);
    }

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogyloc (ci, li);
    }
  return accp;
}

double
relabellocus (int ci, int li, int *mpart)
{
  int ei1;
  int ei2;
  int tries;
  int ismovable1;
  int ismovable2;
  int pop1;
  int fpop1;
  int pop2;
  int fpop2;
  int period1;
  int period2;
  double ktime1;
  double ktime2;
  double logprob;
  int newpop;
  int edge;
  int oldmigcount;
  int newmigcount;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;

  G = &(C[ci]->G[li]);
  gtree = G->gtree;

  IMA_genealogy_save_reset (li);
  /* pick an individual to relabel */
  tries = 0;
  if (modeloptions[NOMIGRATION] == 1)
    {
      /* No migration model is not used here. */
      assert (0);
      ismovable1 = 0;
      ismovable2 = 0;
      while (ismovable1 == 0 || ismovable2 == 0)
        {
          edge = L[li].numgenesknown + randposint (L[li].numgenesunknown);
          ei1 = edge;
          ei2 = L[li].pairs[edge];
          ismovable1 = IMA_edge_movable (ci, li, ei1);
          if (ei2 < 0)
            {
              ismovable2 = 1;
            }
          else
            {
              ismovable2 = IMA_edge_movable (ci, li, ei2);
            }
          if (tries > IM_MAXTRIES_RELABEL)
            {
              break;
            }
          tries++;
        }
    }
  else
    {
      edge = L[li].numgenesknown + randposint (L[li].numgenesunknown);
      ei1 = edge;
      ei2 = L[li].pairs[edge];
      tries++;
    }
  if (tries > IM_MAXTRIES_RELABEL)
    {
      assert (modeloptions[NOMIGRATION] == 1);
      logprob = -DBL_MAX;
      *mpart = 0;
      return logprob;
    }
 
  /* these could be faster because we need to save only edges and its population
   * labels */
  IMA_genealogy_save_edge (ci, li, ei1);
  if (!(ei2 < 0))
    {
      IMA_genealogy_save_edge (ci, li, ei2);
    }

  /* relabel the selected genes: newpop is chosen with logprob */
  if (modeloptions[NOMIGRATION] == 1)
    {
      pop1 = gtree[ei1].pop;
      fpop1 = gtree[gtree[ei1].down].pop;
      if (ei2 < 0)
        {
          logprob = IMA_popn_move (&newpop, ci, pop1, fpop1);
        }
      else
        {
          pop2 = gtree[ei2].pop;
          fpop2 = gtree[gtree[ei2].down].pop;
          ktime1 = gtree[ei1].time;
          period1 = findperiodoftime (ci, ktime1);
          ktime2 = gtree[ei2].time;
          period2 = findperiodoftime (ci, ktime2);
          if (period1 < period2)
            {
              logprob = IMA_popn_move (&newpop, ci, pop1, fpop1);

            }
          else
            {
              logprob = IMA_popn_move (&newpop, ci, pop2, fpop2);
            }
        }
      *mpart = 0;
    }
  else
    {
      pop1 = gtree[ei1].pop;
      fpop1 = gtree[gtree[ei1].down].pop;
      logprob = IMA_popn_move (&newpop, ci, pop1, C[ci]->rootpop);
      assert (logprob == 0.0);

      IMA_edge_fillmiginfo (ci, li, ei1, &saoem);
      if (!(ei2 < 0))
        {
          IMA_edge_fillmiginfo (ci, li, ei2, &saoem2);
        }

      /* not sure if I need this storeoldedges */
      //storeoldedges (ci, li, edge, oldsis, freededge);

      oldmigcount = G->mignum;
      gtree[ei1].mig[0].mt = -1;
      logprob += addmigration_relabel (ci, li, newpop, 
                                       oldmigcount,
                                       &newmigcount,
                                       &saoem, &sanem);
      if (!(ei2 < 0))
        {
          oldmigcount = newmigcount;
          gtree[ei2].mig[0].mt = -1;
          logprob += addmigration_relabel (ci, li, newpop, 
                                           oldmigcount,
                                           &newmigcount,
                                           &saoem2, &sanem2);
        }
      *mpart = newmigcount;
      IMA_edge_copynewmig (ci, li, &sanem);
      if (!(ei2 < 0))
        {
          IMA_edge_copynewmig (ci, li, &sanem2);
        }
    }
 
  gtree[ei1].pop = newpop;
  if (!(ei2 < 0))
    {
      gtree[ei2].pop = newpop;
    }    
 
  return logprob;
}

double
IMA_popn_move (int *topop, int ci, int pop, int apop)
{
  double v;
  int n;
  int ip;

  ci = 0; /* For the use of changing population trees later. */
  *topop = pop;
  n = saC.popnlist[apop][0]->size;
  assert (!(apop < npops));
  while (*topop == pop)
    {
      ip = randposint (n);
      *topop = saC.popnlist[apop][0]->data[ip];
    }
  assert (n > 1);
  /* v = -log ((double) (n - 1)); */
  v = 0.0;
  return v;
}

int
IMA_edge_movable (int ci, int li, int ei)
{
  int v;
  int pop;
  int fpop;
  struct edge *gtree;
  
  v = 0;
  assert (!(ei < 0));
  assert (ei < L[li].numgenes);
  gtree = C[ci]->G[li].gtree;
  pop = gtree[ei].pop;
  fpop = gtree[gtree[ei].down].pop;
  assert (!(pop < 0));
  assert (pop < npops);
  assert (!(fpop < 0));
  if (pop != fpop)
    {
      assert (!(fpop < npops));
      v = 1;
    }
  return v;
}

double
addmigration_relabel (int ci, int li, int newpop, 
                      int oldmigcount,
                      int *newmigcount,
                      struct edgemiginfo *olde, 
                      struct edgemiginfo *newe)
{
  int edge;
  double weight;
  double mproposenum, mproposedenom;
  double mparamf, mparamb;
  struct edge *gtree; 
  int mcount;

  gtree = C[ci]->G[li].gtree;
  assert (C[ci]->G[li].mignum >= 0 );
  /* determine current migration rate to use for update */
  
  /* To reflect Jody's comment: we use overall migration rates of a genealogy.
   * This may help assignment update a little bit. 
   * mparamf = calcmrate (olde->mpall, olde->mtall); */
  mparamf = calcmrate (C[ci]->G[li].mignum, C[ci]->G[li].tlength);

  IMA_reset_edgemiginfo (newe);
  assert(olde->edgeid < L[li].numgenes);
  newe->edgeid = edge = olde->edgeid;
  /* We need li to call checkmig in function simmpath. */
  newe->li = olde->li = li;
  newe->pop = newpop;
  newe->temppop = newpop;
  assert (olde->fpop == gtree[gtree[edge].down].pop);
  newe->fpop = olde->fpop;
  newe->upt = 0.0; /* because we only deal with external branches */
  newe->dnt = gtree[edge].time;
  newe->mig[0].mt = -1;
  fillmiginfoperiods (ci, newe);
  newe->mpall = mwork (ci, newe, newe->e, mparamf);
  assert (newe->mtall > 0);    // some of the edge must occur before the last split time

  mproposedenom = getmprob (ci, mparamf, newe, NULL);
  assert (mproposedenom > -1e200 && mproposedenom < 1e200);

  /* calculate probability of reverse update    */
  mcount = oldmigcount - olde->mpall + newe->mpall;
  
  /* find the migation rate for the backward update */
  /* To reflect Jody's comment: we use overall migration rates of a genealogy.
   * This may help assignment update a little bit. 
   * mparamb = calcmrate (newe->mpall, newe->mtall); */
  mparamb = calcmrate (mcount, C[ci]->G[li].tlength);
  assert (mparamb > 0);

  mproposenum = getmprob (ci, mparamb, olde, NULL);
  assert (mproposenum > -1e200 && mproposenum < 1e200);
  weight = mproposenum - mproposedenom;
  *newmigcount = mcount;
  return weight;
}
void
copy_probcalc_all (struct probcalc *dest, struct probcalc *srce)
{
  memcpy (dest->qintegrate, srce->qintegrate,
          numpopsizeparams * sizeof (double));
  if (nummigrateparams > 0)
    {
      memcpy (dest->mintegrate, srce->mintegrate,
              nummigrateparams * sizeof (double));
    }
  dest->sintegrate = srce->sintegrate;
  dest->pdg = srce->pdg;
  dest->plg = srce->plg;
  dest->probg = srce->probg;
  return;
}
int
updateassignmentbf (int ci)
{
  int accp;
  double lweight;
  double logU;
  double mweight;
  double hweight;
  struct genealogy *G;
  int li;

  G = NULL;
  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }

  /* Count up the number of tries of MCMC cycles. */
  if (ci == 0)
    {
      for (li = 0; li < nloci; li++)
        {
          L[li].a_rec->upinf[IM_UPDATE_ASSIGNMENT_BF].tries++;
        }
    }
  Cupinf[ci].upinf[IM_UPDATE_ASSIGNMENT_BF].tries++;

  /* Save the current state of a chain */
  for (li = 0; li < nloci; li++)
    {
      IMA_genealogy_save_reset (li);
      IMA_savelocusinfo (ci, li);
      copy_treeinfo (&saT[li].savedgweight, &C[ci]->G[li].gweight);
    }
  copy_treeinfo (&saC.savedallgweight, &C[ci]->allgweight);
  copy_probcalc_all (&saC.savedallpcalc, &C[ci]->allpcalc);

  /* Propose a new state of genealogy and assignment. */
  hweight = bfmove (ci);

  /* Find all internal node sequences and mutations of a full genealogy. */
  /* This is for debugging purposes.                                     */
  for (li = 0; li < nloci; li++)
    {
      if (L[li].model == INFINITESITES)
        {
          accp = IMA_genealogy_findIntSeq (ci, li);
          assert (accp == 1);
        }
    }

  /* Compute Metropolis-ratio. */
  setzero_genealogy_weights (&C[ci]->allgweight);
  for (li = 0; li < nloci; li++)
    {
      setzero_genealogy_weights (&C[ci]->G[li].gweight);
      treeweight (ci, li);
      sum_treeinfo (&C[ci]->allgweight, &C[ci]->G[li].gweight);
    }
  initialize_integrate_tree_prob (ci, &C[ci]->allgweight, &C[ci]->allpcalc);

  /* Compute likelihood. */
  /* This has to be called after treeweight. */
  if (calcoptions[DONTCALCLIKELIHOODMUTATION] == 0)
    {
      bflikelihood (ci);
    }
  else
    {
      C[ci]->allpcalc.pdg = 0.0;
    }

  logU = log (uniform ());
  mweight = (C[ci]->allpcalc.probg - saC.savedallpcalc.probg) 
            + gbeta * ((C[ci]->allpcalc.pdg - saC.savedallpcalc.pdg)
            + (C[ci]->allpcalc.plg - saC.savedallpcalc.plg));
  lweight = beta[ci] * mweight + hweight;
  assert (beta[ci] * mweight + hweight > -1e200
          && beta[ci] * mweight + hweight < 1e200);

/* check values of acceptance probability: 
  printf ("pdg: %.3lf - %.3lf (%.3lf), pg: %.3lf - %.3lf (%.3lf), hweight: %.3lf, lweight: %.3lf vs. logU: %.3lf, [%d]\n", 
          C[ci]->allpcalc.pdg, saC.savedallpcalc.pdg, C[ci]->allpcalc.pdg - saC.savedallpcalc.pdg,
          C[ci]->allpcalc.probg, saC.savedallpcalc.probg, C[ci]->allpcalc.probg - saC.savedallpcalc.probg,
          hweight, lweight, logU, lweight > logU);
*/
 
  accp = 0;

  /* Find all internal node sequences and mutations of a full genealogy. */
  /* This is for debugging purposes.                                     */
  for (li = 0; li < nloci; li++)
    {
      if (L[li].model == INFINITESITES)
        {
          accp = IMA_genealogy_findIntSeq (ci, li);
          assert (accp == 1);
        }
      else if (L[li].model == HKY || L[li].model == STEPWISE)
        {
          /* No code. */
        }
    }

  if (lweight > 0.0 || lweight > logU)
    {
      /* accept the update */
      accp = 1;
      if (ci == 0)
        {
          for (li = 0; li < nloci; li++)
            {
              L[li].a_rec->upinf[IM_UPDATE_ASSIGNMENT_BF].accp++;
            }
        }
      Cupinf[ci].upinf[IM_UPDATE_ASSIGNMENT_BF].accp++;

      /* Always revert to an original. */
      for (li = 0; li < nloci; li++)
        {
          if (L[li].model == HKY)
            {
              copyfraclike (ci, li);
              storescalefactors (ci, li);
            }
        }
    }
  else
    {
      for (li = 0; li < nloci; li++)
        {
          /* restore edges changed during update */
          IMA_genealogy_restore_edges (ci, li); 
          IMA_restorelocusinfo (ci, li);
          copy_treeinfo (&C[ci]->G[li].gweight, &saT[li].savedgweight);
        }
      copy_treeinfo (&C[ci]->allgweight, &saC.savedallgweight);
      copy_probcalc_all (&C[ci]->allpcalc, &saC.savedallpcalc);
    }

  for (li = 0; li < nloci; li++)
    {
      if (L[li].model == HKY)
        {
          restorescalefactors (ci, li);
        }
    }

  if (assignmentoptions[POPULATIONASSIGNMENTCHECKPOINT] == 1)
    {
      assertgenealogy (ci);
    }
  return accp;
}

int 
bflikelihood (int ci)
{
  int li;
  int ai;
  struct genealogy *G = NULL;

  assert (calcoptions[DONTCALCLIKELIHOODMUTATION] == 0);

  C[ci]->allpcalc.pdg = 0.0;
  for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      switch (L[li].model)
        {
        case HKY:
          if (assignmentoptions[JCMODEL] == 1)
            {
              G->pdg = likelihoodJC (ci, li, G->uvals[0]);
            }
          else
            {
              G->pdg = likelihoodHKY (ci, li, G->uvals[0], G->kappaval, -1, -1, -1, -1);
            }
          G->pdg_a[0] = G->pdg;
          break;
        case INFINITESITES:
          /* We need both pdg and pdb_a to be set. */
          G->pdg = likelihoodIS (ci, li, G->uvals[0]);
          G->pdg_a[0] = G->pdg;
          break;
        case STEPWISEP:
          G->pdg = 0.0;
          for (ai = 0; ai < L[li].nlinked; ai++)
          {
            /* ui = L[li].uii[ai];       */
            /* assert (pdgnew[ui] == 0); */
            G->pdg_a[ai] = likelihoodSWP (ci, li, ai, G->uvals[ai]);
            G->pdg += G->pdg_a[ai];
          }
          break;
        case STEPWISE:
          G->pdg = 0.0;
          for (ai = 0; ai < L[li].nlinked; ai++)
          {
            G->pdg_a[ai] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
            G->pdg += G->pdg_a[ai];
          }
          break;
        case JOINT_IS_SW:
          assert (0);
          ai = 0;
          G->pdg_a[ai] = likelihoodIS (ci, li, G->uvals[0]);
          for (ai = 1; ai < L[li].nlinked; ai++)
          {
            /* ui = L[li].uii[ai];       */
            /* assert (pdgnew[ui] == 0); */
            G->pdg_a[ai] = likelihoodSW (ci, li, ai, G->uvals[ai], 1.0);
            G->pdg += G->pdg_a[ai];
          }
          break;
        }
      C[ci]->allpcalc.pdg += G->pdg;
    }

  C[ci]->allpcalc.plg = 0.0;
  if (calcoptions[CALCLEVINELIKELIHOOD] == 1)
    {
      C[ci]->allpcalc.plg += likelihoodDG (ci, li);
    }

  return 0; 
}

double
bfmove (int ci)
{
  /* iterator variables */
  int li;
  int ei;
  int pj;
  int ii;
  int ij;
  int ej;
  int i;
  int mi;
  int si;
  /* numbers of counts */
  int nei;
  int ngenes;
  int nsum;
  int np;
  int ki;
  int prevki;
  int oldpop;
  int newpop;
  struct genealogy *G;
  struct edge *gtree; 
  int z1;
  int z3;
  int pc1;
  int pc3;
  char is_coalesced1;
  char is_coalesced3;
  im_edgemiginfo *m1;
  im_edgemiginfo *m3;
  double totalrate;
  double rates[2 * MAXPOPS - 1];
  double newct;
  int down;
  int down1;
  int re;
  double dnt_genealogy;
  int ri;
  int nsim;
  int ti;
  double lhweight;
  im_event *ek;
  im_ginfo *gk;
  im_ginfo *gkz1;
  double prevt;
  double dt;
  double u;
  double rt;
  int nr;
  double lqstar;
  double lq;

  int debug_nsum;

  lhweight = 0.0;

  /* Find all internal node sequences and mutations of a full genealogy. */
  for (li = 0; li < nloci; li++)
    {
      if (L[li].model == INFINITESITES)
        {
          i = IMA_genealogy_findIntSeq (ci, li);
          assert (i == 1);
        }
    }

  /* Pick an individual to relabel. */
  ii = randposint (snind);
  nei = sngenes_ind[ii];
  for (ij = 0; ij < nei; ij++)
    {
      li = saC.indlist[ii][ij].li;
      ei = saC.indlist[ii][ij].ei;
      IMA_genealogy_save_edge (ci, li, ei);
    }

  /* Prepare all migration events storage. */
  for (li = 0; li < nloci; li++)
    {
      IMA_reset_edgemiginfo (&saoems[2 * li]);
      IMA_reset_edgemiginfo (&saoems[2 * li + 1]);
      IMA_reset_edgemiginfo (&saosms[2 * li]);
      IMA_reset_edgemiginfo (&sanems[2 * li]);
      IMA_reset_edgemiginfo (&sanems[2 * li + 1]);
      IMA_reset_edgemiginfo (&sansms[2 * li]);
    }
  for (ij = 0; ij < nei; ij++)
    {
      li = saC.indlist[ii][ij].li;
      ei = saC.indlist[ii][ij].ei;
      ej = saoems[2 * li].edgeid;
      if (ej < 0)
        {
          IMA_edgemiginfo_copymig (&saoems[2 * li], ci, li, ei);
        }
      else if (saoems[2 * li + 1].edgeid < 0)
        {
          /* Diploid is not handled yet. */
          assert (0);
          IMA_edgemiginfo_copymig (&saoems[2 * li + 1], ci, li, ei);
          gtree = C[ci]->G[li].gtree;
          down = gtree[ej].down;
          if (down == gtree[ei].down)
            {
              i = randposint (2);
              IMA_edgemiginfo_appendmig (&saoems[2 * li + i], ci, li, down);
            }
        }
      else
        {
          assert (0);
        }
    }

  /* Label all the genes. */
  li = saC.indlist[ii][0].li;
  ei = saC.indlist[ii][0].ei;
  oldpop = C[ci]->G[li].gtree[ei].pop;
  IMA_popn_move (&newpop, ci, oldpop, C[ci]->rootpop);
  for (ij = 0; ij < nei; ij++)
    {
      li = saC.indlist[ii][ij].li;
      ei = saC.indlist[ii][ij].ei;
      C[ci]->G[li].gtree[ei].pop = newpop;
    }

  /* Simulate branches for each genealogy. */
  for (li = 0; li < nloci; li++)
    {
      ngenes = L[li].numgenes;
      G = &C[ci]->G[li];
      gtree = G->gtree;

      /* If no genes in locus [[li]] belong to individual [[ii]]. */
      /* Decide old branches. -> wrong words!     */
      /* We do not have a longer branch any more. */
      if (saoems[2 * li].edgeid < 0)
        {
          /* No long branches any more. */
          /* assert (0); We can delete this line but I do not allow partial
          genes of an individual. */
          continue;
        }

      assert (saoems[2 * li].upt == 0.0);
      sbfold.m1 = &saoems[2 * li];
      z1 = sbfold.m1->edgeid;

      m1 = &sanems[2 * li]; /* sbfnew.m1 could be m1. */
      m1->li = li;
      m1->edgeid = z1;
      m1->upt = 0.0;
      m1->pop = newpop;
      dnt_genealogy = IMA_genealogy_time (ci, li, G->root);

      /* Detach active lineages. */
      sbfold.is_finiterootbranch = 'F';
      IMA_genealogy_detach (ci, li, &sbfold, z1);

      /* Get intervals of events of a partial genealogy. */
      assert (!gsl_fcmp(dnt_genealogy, G->roottime, 1e-7));
      IMA_intervals_collect (ci, li, 0.0, G->roottime);

      /* Derive parameters from a partial genealogy. */
      IMA_genealogy_derivetheta ('o', ci, li, z1);
      gk = &saGiold[li];
      gkz1 = &saGioldz1[li];

      nsim = 0;
simulation:
      nsim++;
      if (L[li].model != INFINITESITES)
        {
          assert (L[li].model == HKY || L[li].model == STEPWISE);
          assert (nsim == 1);
        }
      /* Prepare the number of lineages with which active lineages can coalesce. */
      IMA_sbf_lineages (&sbfold, ci, li, z1);


      /* Initialize active lineages. */
      m1->dnt = -1.0;
      m1->fpop = -1;
      m1->mpall = 0;
      m1->mig[0].mt = -1.0;
      m1->mig[0].mp = -1;
      sbfnew.sis1 = -1;
      sbfnew.m1 = NULL;
      sbfnew.m3 = NULL;
      sbfnew.is_finiterootbranch = 'F';
      sbfnew.canz31 = 'F';
      sbfold.canz31 = 'F';
      sbfold.m3 = NULL;
      is_coalesced1 = 'F';
      ti = 0;
      prevt = 0.0;
      ki = 0;
      pc1 = gtree[z1].pop;
      assert (pc1 == newpop);

      /* Start simulation down towards the genealogy. */
      ek = sbfold.e1; 
      while (is_coalesced1 == 'F')
        {
          dt = ek->t - prevt;
          assert (dt > 0.0);
          while (dt > 0.0 && is_coalesced1 == 'F')
            {
              np = bfmove_ratesI (&rates, gk, gkz1,
                                  &sbfold, ki,
                                  pc1);
              totalrate = 0.0;
              for (pj = 0; pj < np; pj++)
                {
                  assert (!(rates[pj] < 0.0));
                  totalrate += rates[pj];
                }
              if (totalrate > 0.0)
                {
                  u = uniform ();
                  rt = -log (u) / totalrate;
                  dt -= rt; 
                }
              else if (totalrate == 0.0)
                {
                  /* assert (0); */
                  /* assert (0); I do not know when this would happen: Answer -
                   * This can happen when there is no lineage to attach an
                   * active lineage in the last period. I am not sure whether
                   * this could continue to happen once we have zero total rate.
                   * If that is the case, we may let this happen and resimulate
                   * once we get to the beyond of the bottom of a partial
                   * genealogy. */
                  dt = -1.0;
                }
              else
                {
                  assert (0);
                }
              if (dt < 0)
                {
                  /* move to next interval: we hit and exit the current while,
                   * then we increase ti and ek */
                  continue;
                }
              else
                {
                  newct = ek->t - dt;
                }
              re = ran_discretea (np, rates);

              /* Case (A). Only one lineage z1 coalesces. */
              if ((np == 1 && is_coalesced1 == 'F')
                  || (np == 2 && is_coalesced1 == 'F' && re == 0))
                {
                  is_coalesced1 = 'T';
                  /* We need a list of edges: n and l.                     */
                  /* We computes discrete probabilities for all the edges. */
                  i = bfmove_chooseAnEdge (ci, li, &sbfold, newct, pc1);
                  if (L[li].model == INFINITESITES)
                    {
                      sbfnew.sis1 = sbfold.lz1[pc1][i];
                    }
                  else if (L[li].model == HKY || L[li].model == STEPWISE)
                    {
                      sbfnew.sis1 = sbfold.l[pc1][i];
                    }
                  else
                    {
                      assert (0);
                    }
                  /* sbfnew.down1 = sbfold.down1; */
                  sbfnew.m1 = m1;
                  m1->mig[m1->mpall].mt = -1.0;
                  m1->dnt = newct;
                  m1->fpop = pc1;
                }
              /* Case (F). Lineage z1 migrates */
              else if (np == 2 && is_coalesced1 == 'F' && re == 1)
                {
                  IMA_choose_migtopop (&ri, ci, li, ki, pc1, gk);
                  assert (pc1 != ri);
                  pc1 = ri;
                  mi = m1->mpall;
                  m1->mig[mi].mt = newct;
                  m1->mig[mi].mp = pc1;
                  m1->mpall++;
                }
              else
                {
                  assert (0);
                }
            }
          if (is_coalesced1 == 'T')
            {
              break;
            }

          /* New time is further back in time than the bottom of an event. */
          bfmove_nextI (ci, li, 0.0, ek, np, rates, 
                        NULL, &sbfold, &sbfnew,
                        &ki,
                        z1,
                        &is_coalesced1,
                        &pc1);
          prevt = ek->t;
          ti++;
          if (ti < sbfold.n1)
            {
              ek++;
            }
          else
            {
              break;
            }
        }

      /* Simulation beyond the bottom of a genealogy. */
      sbfnew.m3 = NULL;
      if (is_coalesced1 == 'F') 
        {
          assert (sbfold.n1 == ti);
          /* Should we do two different actions for the last event type? */
          assert (sbfold.e1[ti - 1].type == 'c' 
                  || sbfold.e1[ti - 1].type == 'e');
          assert (sbfold.e1[ti - 1].t == prevt);
          assert (findperiodoftime (ci, prevt) == ki);

          is_coalesced3 = 'F';
          newct = prevt;
          pc3 = sbfold.e1[ti - 1].pj;
          assert (!(pc1 < 0));
          assert (!(pc3 < 0));
          assert (sbfold.e1[ti - 1].ei == G->root);
          z3 = G->root;
          m3 = &sansms[2 * li];
          m3->li = li;
          m3->edgeid = z3;
          m3->upt = prevt;
          m3->pop = pc3;

          sbfold.canz31 = 'T';
          if (L[li].model == INFINITESITES)
            {
              for (si = 0; si < L[li].numsites; si++)
                {
                  if (G->mut[si] == z3)
                    {
                      continue;
                    }
                  if (gtree[z3].seq[si] != sbfold.seqz1[si])
                    {
                      sbfold.canz31 = 'F';
                      break;
                    }
                }
            }
          else if (L[li].model == HKY || L[li].model == STEPWISE)
            {
              /* No code: we should attach an active lineage. */
            }

          nsum = IMA_sbf_nsum (&sbfold, ci, li, pc1, pc3);
          assert (nsum == 2); /* Do not consider diploid data type. */
          if (is_coalesced1 == 'F' && sbfold.canz31 == 'F')
            {
              m1->dnt = -1.0;
              m1->fpop = -1;
              m1->mpall = 0;
              m1->mig[0].mt = -1.0;
              goto simulation;
            }
          assert (sbfold.canz31 == 'T');

          /* ki is set before. */
          while (nsum == 2)
            {
              assert (is_coalesced1 == 'F');
              assert (is_coalesced3 == 'F');
              pc1 = saC.popndown[pc1][ki];
              pc3 = saC.popndown[pc3][ki];

              bfmove_ratesIII (&rates, gk, &sbfold, ki, pc1, pc3);
              nr = 3;
              totalrate = 0.0;
              for (ri = 0; ri < nr; ri++)
                {
                  totalrate += rates[ri];
                }
              assert (totalrate > 0.0);
              u = uniform ();
              rt = -log (u) / totalrate;
              prevki = ki;
              newct = newct + rt;
              ki = findperiodoftime (ci, newct);
              if (prevki == ki)
                {
                  /* No Code */
                }
              else if (prevki < ki)
                {
                  /* ki can be larger than prevki by more than one. */
                  newct = C[ci]->tvals[prevki];
                  ki = prevki + 1;
                  continue;
                }
              else
                {
                  assert (0);
                }
              re = ran_discretea (nr, rates);
              switch (re)
                {
                case 0:
                  IMA_choose_migtopop (&ri, ci, li, ki, pc1, gk);
                  pc1 = ri;
                  mi = m1->mpall;
                  m1->mig[mi].mt = newct;
                  m1->mig[mi].mp = pc1;
                  m1->mpall++;
                  break;
                case 1:
                  IMA_choose_migtopop (&ri, ci, li, ki, pc3, gk);
                  pc3 = ri;
                  mi = m3->mpall;
                  m3->mig[mi].mt = newct;
                  m3->mig[mi].mp = pc3;
                  m3->mpall++;
                  break;
                case 2: /* Lineages z1 coalesces with z3. */
                  is_coalesced3 = 'T';
                  is_coalesced1 = 'T';
                  /* What are all these about? */
                  sbfnew.sis1 = z3;
                  /* sbfnew.down1 = sbfold.down1; */
                  sbfnew.m1 = m1;
                  m1->mig[m1->mpall].mt = -1.0;
                  m1->dnt = newct;
                  m1->fpop = pc1;
                  sbfnew.m3 = m3;
                  m3->mig[m3->mpall].mt = -1.0;
                  m3->dnt = newct;
                  m3->fpop = pc3;
                  assert (pc1 == pc3);
                  bfmove_chooseA1 (ci, li, &sbfnew);
                  nsum = 1;
                  break;
                }
            }
        }

      /* Get [[sbfnew]] and [[saGinew]]. */
      IMA_intervals_collectsbfnew (ci, li, sbfnew.m3);
      IMA_genealogy_derivetheta ('n', ci, li, z1);
      if (L[li].model == INFINITESITES)
        {
          /* No code. */
        }
      else if (L[li].model == HKY || L[li].model == STEPWISE)
        {
          /* FIXME: What if there is no or single inactive lineage? */
          bfmove_selectAnEdge (ci, li, &sbfnew, &sbfold);
          /* bfmove_selectA1 (ci, li, &sbfold); */
        }

      /* Compute Hasting ratio. */
      /* G     -> Gp     -> Gstar */
      lqstar = bfq (&sbfold, &sbfnew, ci, li, gk, gkz1);
      /* Gstar -> Gstarp -> G     */
      lq = bfq (&sbfnew, &sbfold, ci, li, &saGinew[li], &saGinewz1[li]);
      lhweight += (lq - lqstar);

      /* Reattach the branch or branches. */
      down1 = gtree[sbfnew.m1->edgeid].down;
      if (sbfnew.m3 != NULL)
        {
          assert (sbfnew.sis1 == G->root);
          IMA_genealogy_appendmig (ci, li, 
                                   sbfnew.sis1, sbfnew.m3, 
                                   sbfnew.m1->edgeid, down1);
        }
      else
        {
          IMA_genealogy_bisectsis (ci, li, 
                                   sbfnew.sis1, sbfnew.m1->edgeid, 
                                   down1, 
                                   sbfnew.m1, 
                                   sbfnew.seqz1);
        }
      IMA_genealogy_copymig (ci, li, sbfnew.m1->edgeid, sbfnew.m1);
      gtree[G->root].time = TIMEMAX;
      gtree[G->root].down = -1;
      gtree[G->root].mig[0].mt = -1.0;

      /* Update locus information. */ 
      gtree[z1].exist = 'T';
      down1 = gtree[z1].down;
      gtree[down1].exist = 'T';
      G->roottime = IMA_genealogy_time (ci, li, G->root); 
      sbfold.n1 = 0;
      sbfold.n2 = 0;
      sbfold.n3 = 0;
      sbfnew.n1 = 0;
      sbfnew.n2 = 0;
      sbfnew.n3 = 0;

      if (L[li].model == INFINITESITES)
        {
          i = IMA_genealogy_findIntSeq (ci, li);
          assert (i == 1);
        }
      else if (L[li].model == HKY) 
        {
          /* No code. */
        }
      else if (L[li].model == STEPWISE)
        {
          /* No code. */
          assert (sbfnew.A1 >= L[li].minA[0] && sbfnew.A1 <= L[li].maxA[0]);
          gtree[down1].A[0] = sbfnew.A1;
        }
      /* printf ("[%d] root time - %lf\n", step, G->roottime);  */
    }

  return lhweight;
}

void
bfmove_chooseA1 (int ci, int li, im_bfupdate *pbf)
{
  double u;
  struct genealogy *G;
  struct edge *gtree; 
  int nA;
  int bi;
  int newA;
  double dlikeA;
  double sisupt;
  int sis;
  int edge;
  double ct;
  int vj;
  double sisnewtlen;
  int d;
  double bessiv;

  if (L[li].model == STEPWISE)
    {
      G = &C[ci]->G[li];
      gtree = G->gtree;

      /**********************************************************/
      /*                  SOMEWHAT COMPLICATED                  */
      /**********************************************************/
      /* All unlinked loci! */
      edge = pbf->m1->edgeid;
      sis = pbf->m3->edgeid;
      sisupt = IMA_genealogy_time (ci, li, sis);
      sisnewtlen = pbf->m3->dnt - sisupt;
      ct = pbf->m1->dnt;
      
      u = G->uvals[0];
      nA = L[li].maxA[0] - L[li].minA[0] + 1;
      for (bi = 0; bi < nA; bi++)
        {
          newA = L[li].minA[0] + bi;
          dlikeA = 0.0; /* FIXME: */

          d = gtree[sis].A[0] - newA;
          bessiv = bessi (d, sisnewtlen * u);
          assert (bessiv > 0.0);
          dlikeA += (-(sisnewtlen * u) + log (bessiv));

          d = gtree[edge].A[0] - newA;
          bessiv = bessi (d, ct * u);
          assert (bessiv > 0.0);
          dlikeA += (-(ct * u) + log (bessiv));

          /* A rooting edge does not have a likelihood under SW model. */

          pbf->likes11[0][bi] = dlikeA;
        }
      normalizeLoga (nA, pbf->likes11[0]);
      vj = ran_discreteloga (nA, pbf->likes11[0]);
      pbf->likei11 = vj;
      pbf->A1 = L[li].minA[0] + vj;
      pbf->likei1 = 0;
      pbf->likes1[0] = 0.0;
      /**********************************************************/
      /*                  SOMEWHAT COMPLICATED                  */
      /**********************************************************/
    }
  return;
}

double
bfmove_selectA1 (int ci, int li, im_bfupdate *pbf, double lq)
{
  double u;
  struct genealogy *G;
  struct edge *gtree; 
  int nA;
  int bi;
  int newA;
  double dlikeA;
  double sisupt;
  int sis;
  int edge;
  double ct;
  int vj;
  double sisnewtlen;
  int d;
  double bessiv;

  if (L[li].model == STEPWISE)
    {
      assert (pbf->m3 != NULL);

      G = &C[ci]->G[li];
      gtree = G->gtree;

      /**********************************************************/
      /*                  SOMEWHAT COMPLICATED                  */
      /**********************************************************/
      /* All unlinked loci! */
      edge = pbf->m1->edgeid;
      sis = pbf->m3->edgeid;
      sisupt = IMA_genealogy_time (ci, li, sis);
      sisnewtlen = pbf->m3->dnt - sisupt;
      ct = pbf->m1->dnt;
      
      u = G->uvals[0];
      nA = L[li].maxA[0] - L[li].minA[0] + 1;
      for (bi = 0; bi < nA; bi++)
        {
          newA = L[li].minA[0] + bi;
          dlikeA = 0.0; /* FIXME: */

          d = gtree[sis].A[0] - newA;
          bessiv = bessi (d, sisnewtlen * u);
          assert (bessiv > 0.0);
          dlikeA += (-(sisnewtlen * u) + log (bessiv));

          d = gtree[edge].A[0] - newA;
          bessiv = bessi (d, ct * u);
          assert (bessiv > 0.0);
          dlikeA += (-(ct * u) + log (bessiv));

          /* A rooting edge does not have a likelihood under SW model. */
          pbf->likes11[0][bi] = dlikeA;

          if (pbf->A1 == L[li].minA[0] + bi)
            {
              pbf->likei11 = bi;
            }
        }
      normalizeLoga (nA, pbf->likes11[0]);
      lq += pbf->likes11[0][pbf->likei11];
      assert (pbf->A1 == L[li].minA[0] + pbf->likei11);
      /**********************************************************/
      /*                  SOMEWHAT COMPLICATED                  */
      /**********************************************************/
    }
  return lq;
}

int 
bfmove_chooseAnEdge (int ci, int li, im_bfupdate *pbf, double ct, int pc)
{
  int i;
  int vi;
  int sis;
  struct genealogy *G;
  struct edge *gtree; 
  int edge;
  int down;
  int downdown;
  int up;
  double sisupt;
  double sisdnt;
  double sisoldtl;
  double sisnewtl;
  double downtl;
  double u;
  int nA;
  int bi;
  int newA;
  double dlikeA;
  int d;
  double bessiv;
  int iszero;
  int vj;

  int ai;

  G = &C[ci]->G[li];
  gtree = G->gtree;
  edge = pbf->m1->edgeid;
  down = gtree[edge].down;

  if (L[li].model == INFINITESITES)
    {
      assert (pbf->nz1[pc] > 0);
      vi = randposint (pbf->nz1[pc]);
    }
  else if (L[li].model == HKY || L[li].model == STEPWISE)
    {
      assert (pbf->n[pc] > 0);
#ifdef NOTCONSIDERINGDATA
      vi = randposint (pbf->n[pc]);
      break;
#endif /* NOTCONSIDERINGDATA */
      if (pbf->n[pc] > 0)
        {
          /* Compute all likelihoods for possible full genealogies. */
          
          /* We attach an active lineage to one of lineages. Each of the lineages is
           * a candidate for lineage [[sis]]. We freely attach it, 
           * compute likelihood, and detach it. We keep doing this until we fill all
           * likelihoods for all full possible genealogies. */
          for (i = 0; i < pbf->n[pc]; i++)
            {
              sis = pbf->l[pc][i]; /* FIXME: l is not set yet. */
              /* Attach it. */ 
              /* Attaching may be done using a function. */
              downdown = gtree[sis].down;
              sisupt = IMA_genealogy_time (ci, li, sis);
              assert (sisupt < ct); 

              sisdnt = gtree[sis].time;
              if (G->root == sis)
                {
                  assert (pbf->n[pc] == 1);
                  /* sisdnt is not ... */
                  /* sisdnt = gtree[sis].time; */
                }
              else
                {
                  assert (ct < sisdnt);
                }

              if (G->root == sis)
                {
                  assert (downdown == -1);
                  G->root = down;
                  G->roottime = ct;
                }
              else
                {
                  if (gtree[downdown].up[IM_EDGE_LEFT] == sis)
                    {
                      gtree[downdown].up[IM_EDGE_LEFT] = down;
                    }
                  else
                    {
                      gtree[downdown].up[IM_EDGE_RIGHT] = down;
                    }
                }
              gtree[down].down = downdown;
              gtree[sis].down = down;
              gtree[edge].down = down;
              gtree[down].up[IM_EDGE_LEFT] = sis;
              gtree[down].up[IM_EDGE_RIGHT] = edge;
              gtree[edge].time = ct;
              gtree[sis].time = ct;
              gtree[down].time = sisdnt;
              sisoldtl = sisdnt - sisupt;
              sisnewtl = ct - sisupt;
              downtl = sisdnt - ct;

              /* Compute likelihood. */
              if (L[li].model == HKY)
                {
                  pbf->likes1[i] = likelihoodHKY (ci, li, 
                                                  G->uvals[0], G->kappaval, 
                                                  -1, -1, -1, -1);
                }
              else if (L[li].model == STEPWISE)
                {
                  /**********************************************************/
                  /*                  SOMEWHAT COMPLICATED                  */
                  /**********************************************************/
                  /* All unlinked loci! */
                  u = G->uvals[0];
                  pbf->likes1[i] = 0.0;
                  nA = L[li].maxA[0] - L[li].minA[0] + 1;
                  for (bi = 0; bi < nA; bi++)
                    {
                      newA = L[li].minA[0] + bi;
                      dlikeA = 0.0; /* FIXME: */
                      if (G->root != down)
                        {
                          dlikeA -= gtree[sis].dlikeA[0];
                        }
                      else
                        {
                          assert (downdown == -1);
                        }

                      d = gtree[sis].A[0] - newA;
                      bessiv = bessi (d, sisnewtl * u);
                      if (!(bessiv > 0.0))
                        {
                          iszero = 1;
                          assert (0);
                          /* Impossible situation: IM_BESSI_MIN; */
                        }
                      else
                        {
                          dlikeA += (-(sisnewtl * u) + log (bessiv));
                        }
                      d = gtree[edge].A[0] - newA;
                      bessiv = bessi (d, ct * u);
                      if (!(bessiv > 0.0))
                        {
                          iszero = 1;
                          assert (0);
                          /* Impossible situation: IM_BESSI_MIN; */
                        }
                      else
                        {
                          dlikeA += (-(ct * u) + log (bessiv));
                        }
 
                      if (!(downdown < 0))
                        {
                          d = newA - gtree[downdown].A[0];
                          bessiv = bessi (d, downtl * u);
                          if (!(bessiv > 0.0))
                            {
                              iszero = 1;
                              assert (0);
                              /* Impossible situation: IM_BESSI_MIN; */
                            }
                          else
                            {
                              dlikeA += (-(downtl * u) + log (bessiv));
                            }
                        }
                      pbf->likes11[i][bi] = dlikeA;
                      if (bi == 0)
                        {
                          pbf->likes1[i] = dlikeA;
                        }
                      else
                        {
                          LogSum2 (pbf->likes1[i], pbf->likes1[i], dlikeA);
                        }
                    }
                  /**********************************************************/
                  /*                  SOMEWHAT COMPLICATED                  */
                  /**********************************************************/
                }
              else if (L[li].model == STEPWISEP)
                {
                  pbf->likes1[i] = 0.0;
                  for (ai = 0; ai < L[li].nlinked; ai++)
                    {
                      pbf->likes1[i] += likelihoodSWP (ci, li, 
                                                       ai, G->uvals[ai]);
                    }
                }


              /* Detach it. */
              up = sis;
              if (G->root == down)
                {
                  assert (downdown == -1);
                  gtree[up].down = -1;
                  G->root = up;
                  gtree[up].time = sisdnt;
                }
              else
                {
                  assert (!(downdown < 0));
                  gtree[up].down = downdown;
                  if (gtree[downdown].up[IM_EDGE_LEFT] == down)
                    {
                      gtree[downdown].up[IM_EDGE_LEFT] = up;
                    }
                  else
                    {
                      gtree[downdown].up[IM_EDGE_RIGHT] = up;
                    }
                  gtree[up].time = gtree[down].time;
                }
            }

          if (pbf->n[pc] > 1)
            {
              normalizeLoga (pbf->n[pc], pbf->likes1);
              vi = ran_discreteloga (pbf->n[pc], pbf->likes1);
              pbf->likei1 = vi;
            }
          else if (pbf->n[pc] == 1)
            {
              vi = 0;
              pbf->likes1[0] = 0.0;
              pbf->likei1 = vi;
            }
          else
            {
              assert (0);
            }

          if (L[li].model == STEPWISE)
            {
              nA = L[li].maxA[0] - L[li].minA[0] + 1;
              normalizeLoga (nA, pbf->likes11[vi]);
              vj = ran_discreteloga (nA, pbf->likes11[vi]);
              pbf->likei11 = vj;
            }
        }
      else if (pbf->n[pc] == 1)
        {
          assert (0); /* Never be called. */
          vi = 0;
          pbf->likes1[0] = 0.0;
          pbf->likei1 = vi;
  
          /* FIXME: likes11[0] is the array that we need to fill. */
          if (L[li].model == STEPWISE)
            {
              nA = L[li].maxA[0] - L[li].minA[0] + 1;
              normalizeLoga (nA, pbf->likes11[vi]);
              vj = ran_discreteloga (nA, pbf->likes11[vi]);
              pbf->likei11 = vj;
            }
        }
      else
        {
          assert (0);
        }
    }
  else
    {
      assert (0);
    }

  return vi;
}

int 
bfmove_selectAnEdge (int ci, int li, im_bfupdate *pbf, im_bfupdate *qbf)
{
  /* pbf: sbfnew, qbf: sbfold */
  int pc;
  double ct;
  int ki;
  int mi;
  int pop;
  int i;
  int n;
  int vi;
  int sis;
  struct genealogy *G;
  struct edge *gtree; 
  int edge;
  int down;
  int downdown;
  int up;
  double sisupt;
  double sisdnt;
  double upt;
  double dnt;
  int ngnodes;
  int ngenes;
  char iscross;
  int ei;
  int bi;
  int nA;
  double dlikeA;
  int d;
  double u;
  int newA;
  double sisoldtl;
  double sisnewtl;
  double downtl;
  double bessiv;
  int iszero;

  G = &C[ci]->G[li];
  gtree = G->gtree;
  edge = pbf->m1->edgeid;
  down = gtree[edge].down;
  ct = qbf->m1->dnt;
  ki = findperiodoftime (ci, ct);
  pc = qbf->m1->fpop;

  ngnodes = L[li].numlines;
  ngenes = L[li].numgenes;

  /* Find all edges that are crossed by time [[ct]]. Exclude edges that are not
   * a part of a partial genealogy. */
  pbf->n[pc] = 0;
  if (qbf->is_finiterootbranch == 'T')
    {
      pbf->n[pc] = 1;
      pbf->l[pc][0] = G->root;
    }
  else
    {
      for (ei = 0; ei < ngnodes; ei++)
        {
          /* Exclude a gene at an active lineage. */
          if (gtree[ei].exist == 'F')
            {
              continue; 
            }
          assert (!(gtree[ei].pop < 0));
          assert (gtree[ei].pop < 2 * MAXPOPS - 1);
          if (ei < ngenes)
            {
              upt = 0.0;
            }
          else
            {
              up = gtree[ei].up[0];
              upt = gtree[up].time;
            }
          dnt = gtree[ei].time;
          iscross = 'F';
          if (upt < ct && ct < dnt)
            {
              pop = gtree[ei].pop;
              if (gtree[ei].mig[0].mt < 0.0)
                {
                  pop = saC.popndown[pop][ki];
                  if (pop == pc)
                    {
                      iscross = 'T';
                    }
                }
              else
                {
                  mi = 0;
                  while (gtree[ei].mig[mi].mt > -0.5)
                    {
                      dnt = gtree[ei].mig[mi].mt;
                      if (upt < ct && ct < dnt)
                        {
                          pop = saC.popndown[pop][ki];
                          if (pop == pc)
                            {
                              iscross = 'T';
                            }
                          break;
                        }
                      upt = dnt;
                      pop = gtree[ei].mig[mi].mp;
                      mi++;
                    }
                  if (iscross == 'F')
                    {
                      dnt = gtree[ei].time;
                      if (upt < ct && ct < dnt)
                        {
                          pop = saC.popndown[pop][ki];
                          if (pop == pc)
                            {
                              iscross = 'T';
                            }
                        }
                    }
                }
            }
          if (iscross == 'T')
            {
              n = pbf->n[pc];
              pbf->l[pc][n] = ei;
              pbf->n[pc]++;
            }
        }
    }

  if (L[li].model == INFINITESITES)
    {
      assert (0);
    }
  else if (L[li].model == HKY || L[li].model == STEPWISE)
    {
      assert (pbf->n[pc] > 0);
#ifdef NOTCONSIDERINGDATA
      vi = randposint (pbf->n[pc]);
      break;
#endif /* NOTCONSIDERINGDATA */
      if (pbf->n[pc] > 0)
        {
          /* Compute all likelihoods for possible full genealogies. */
          
          /* We attach an active lineage to one of lineages. Each of the lineages is
           * a candidate for lineage [[sis]]. We freely attach it, 
           * compute likelihood, and detach it. We keep doing this until we fill all
           * likelihoods for all full possible genealogies. */
          for (i = 0; i < pbf->n[pc]; i++)
            {
              sis = pbf->l[pc][i]; /* FIXME: l is not set yet. */
              /* Attach it. */ 
              /* Attaching may be done using a function. */
              downdown = gtree[sis].down;
              sisupt = IMA_genealogy_time (ci, li, sis);
              assert (sisupt < ct); 

              sisdnt = gtree[sis].time;
              if (G->root == sis)
                {
                  assert (pbf->n[pc] == 1);
                  /* sisdnt is not ... */
                  /* sisdnt = gtree[sis].time; */
                }
              else
                {
                  assert (ct < sisdnt);
                }

              if (G->root == sis)
                {
                  assert (downdown == -1);
                  G->root = down;
                  G->roottime = ct;
                }
              else
                {
                  if (gtree[downdown].up[IM_EDGE_LEFT] == sis)
                    {
                      gtree[downdown].up[IM_EDGE_LEFT] = down;
                    }
                  else
                    {
                      gtree[downdown].up[IM_EDGE_RIGHT] = down;
                    }
                }
              gtree[down].down = downdown;
              gtree[sis].down = down;
              gtree[edge].down = down;
              gtree[down].up[IM_EDGE_LEFT] = sis;
              gtree[down].up[IM_EDGE_RIGHT] = edge;
              gtree[edge].time = ct;
              gtree[sis].time = ct;
              gtree[down].time = sisdnt;
              sisoldtl = sisdnt - sisupt;
              sisnewtl = ct - sisupt;
              downtl = sisdnt - ct;

              /* Compute likelihood. */
              if (L[li].model == HKY)
                {
                  pbf->likes1[i] = likelihoodHKY (ci, li, 
                                                  G->uvals[0], G->kappaval, 
                                                  -1, -1, -1, -1);
                }
              else if (L[li].model == STEPWISE)
                {
                  /**********************************************************/
                  /*                  SOMEWHAT COMPLICATED                  */
                  /**********************************************************/
                  /* All unlinked loci! */
                  u = G->uvals[0];
                  pbf->likes1[i] = 0.0;
                  nA = L[li].maxA[0] - L[li].minA[0] + 1;
                  for (bi = 0; bi < nA; bi++)
                    {
                      newA = L[li].minA[0] + bi;
                      dlikeA = 0.0; /* FIXME: */
                      if (G->root != down)
                        {
                          dlikeA -= gtree[sis].dlikeA[0];
                        }
                      else
                        {
                          assert (downdown == -1);
                        }

                      d = gtree[sis].A[0] - newA;
                      bessiv = bessi (d, sisnewtl * u);
                      if (!(bessiv > 0.0))
                        {
                          iszero = 1;
                          assert (0);
                          /* Impossible situation: IM_BESSI_MIN; */
                        }
                      else
                        {
                          dlikeA += (-(sisnewtl * u) + log (bessiv));
                        }
                      d = gtree[edge].A[0] - newA;
                      bessiv = bessi (d, ct * u);
                      if (!(bessiv > 0.0))
                        {
                          iszero = 1;
                          assert (0);
                          /* Impossible situation: IM_BESSI_MIN; */
                        }
                      else
                        {
                          dlikeA += (-(ct * u) + log (bessiv));
                        }
 
                      if (!(downdown < 0))
                        {
                          d = newA - gtree[downdown].A[0];
                          bessiv = bessi (d, downtl * u);
                          if (!(bessiv > 0.0))
                            {
                              iszero = 1;
                              assert (0);
                              /* Impossible situation: IM_BESSI_MIN; */
                            }
                          else
                            {
                              dlikeA += (-(downtl * u) + log (bessiv));
                            }
                        }
                      pbf->likes11[i][bi] = dlikeA;
                      if (bi == 0)
                        {
                          pbf->likes1[i] = dlikeA;
                        }
                      else
                        {
                          LogSum2 (pbf->likes1[i], pbf->likes1[i], dlikeA);
                        }
                    }
                  /**********************************************************/
                  /*                  SOMEWHAT COMPLICATED                  */
                  /**********************************************************/
                }

              if (sis == qbf->sis1) /* qbf: sbfold, where is sbfold.sis1 set? */
                {
                  /* IMA_genealogy_detach (ci, li, &sbfold, z1): sets sis1. */
                  pbf->likei1 = i;
                  if (L[li].model == STEPWISE)
                    {
                      pbf->likei11 = gtree[down].A[0] - L[li].minA[0];
                    }
                }

              /* Detach it. */
              up = sis;
              if (G->root == down)
                {
                  assert (downdown == -1);
                  gtree[up].down = -1;
                  G->root = up;
                  gtree[up].time = sisdnt;
                }
              else
                {
                  assert (!(downdown < 0));
                  gtree[up].down = downdown;
                  if (gtree[downdown].up[IM_EDGE_LEFT] == down)
                    {
                      gtree[downdown].up[IM_EDGE_LEFT] = up;
                    }
                  else
                    {
                      gtree[downdown].up[IM_EDGE_RIGHT] = up;
                    }
                  gtree[up].time = gtree[down].time;
                }
            }

          /* We need to set likei11. */
          /* vi = ran_discreteloga (pbf->n[pc], pbf->likes1); */

                  if (pbf->n[pc] > 1)
                    {
                          normalizeLoga (pbf->n[pc], pbf->likes1);
                    }
          else if (pbf->n[pc] == 1)
                    {
              pbf->likes1[0] = 0.0;
              pbf->likei1 = 0;
            }
          else
            {
              assert (0);
            }

          if (L[li].model == STEPWISE)
                    {
              nA = L[li].maxA[0] - L[li].minA[0] + 1;
              assert (!(pbf->likei11 < 0));
              normalizeLoga (nA, pbf->likes11[pbf->likei1]);
            }
        }
      else if (pbf->n[pc] == 1)
        {
          assert (0);
          pbf->likes1[0] = 0.0;
          pbf->likei1 = 0;
        }
    }
  else
    {
      assert (0);
    }

  return vi;
}

void 
IMA_sbf_correctnl (int ci, int li, int sis, int *n, int **l,
                   int *seqzu, int *seqzl)
{
  int i;
  int p;
  int *seqz; 
  int *downseqz;
  struct genealogy *G;
  struct edge *gtree; 
  int si;
  int down;
  char is_sisInL;
  char remove_sisInL;
  int p_sisInL;
  int i_sisInL; 

  G = &C[ci]->G[li];
  gtree = G->gtree;

  down = gtree[sis].down;
  seqz = gtree[sis].seq;
  if (down == -1)
    {
      downseqz = malloc (L[li].numsites * sizeof (int));
      for (si = 0; si < L[li].numsites; si++)
        {
          if (G->mut[si] == sis)
            {
              if (seqz[si] == 0)
                {
                  downseqz[si] = 1;
                }
              else if (seqz[si] == 1)
                {
                  downseqz[si] = 0;
                }
              else
                {
                  assert (0);
                }
            }
          else
            {
              downseqz[si] = seqz[si];
            }
        }
    }
  else
    {
      downseqz = gtree[down].seq;
      for (si = 0; si < L[li].numsites; si++)
        {
          if (G->mut[si] == sis)
            {
              assert (seqz[si] != downseqz[si]);
            }
          else
            {
              assert (seqz[si] == downseqz[si]);
            }
        } 
    }
 
  is_sisInL = 'F';
  for (p = 0; p < numpopsizeparams; p++)
    {
      for (i = 0; i < n[p]; i++)
        {
          if (l[p][i] == sis)
            {
              assert (is_sisInL == 'F');
              is_sisInL = 'T';
              p_sisInL = p;
              i_sisInL = i;
            }
        }
    }
  if (is_sisInL == 'T')
    {
      remove_sisInL = 'F';
      for (si = 0; si < L[li].numsites; si++)
        {
          if (G->mut[si] == sis)
            {
              assert (seqz[si] != downseqz[si]);
              if (seqz[si] == seqzu[si])
                {
                  /* No Code! */
                }
              else
                {
                  if (seqzu[si] == seqzl[si])
                    {
                      /* No Code! */
                    }
                  else
                    {
                      remove_sisInL = 'T';
                      break;
                    }
                }
              assert (seqz[si] != downseqz[si]);
            }
          else
            {
              assert (seqz[si] == seqzu[si]);
              assert (seqz[si] == seqzl[si]);
              assert (seqz[si] == downseqz[si]);
            }
        } 
    }
  if (remove_sisInL == 'T')
    {
      for (i = i_sisInL; i < n[p_sisInL] - 1; i++)
        {
          l[p_sisInL][i] = l[p_sisInL][i + 1];
        }
      n[p_sisInL]--;
    }
  
  if (down == -1)
    {
      free (downseqz);
      downseqz = NULL;
    }
  return;
}
char
IMA_sbf_correctn3 (int ci, int li, int sis, 
                   int *seqzu, int *seqzl)
{
  char canzxy;
  int *seqz; 
  struct genealogy *G;
  struct edge *gtree; 
  int si;

  G = &C[ci]->G[li];
  gtree = G->gtree;
  canzxy = 'T';
  seqz = gtree[sis].seq;
  for (si = 0; si < L[li].numsites; si++)
    {
      if (G->mut[si] == sis)
        {
          if (seqz[si] == seqzu[si])
            {
              /* No Code! */
            }
          else
            {
              if (seqzu[si] == seqzl[si])
                {
                  /* No Code! */
                }
              else
                {
                  canzxy = 'F';
                  break;
                }
            }
        }
      else
        {
          assert (seqz[si] == seqzu[si]);
          if (seqz[si] == seqzl[si])
            {
              /* No Code! */
            }
          else
            {
              canzxy = 'F';
              break;
            }
        }
    } 

  return canzxy;
}

double 
bfmove_nextI (int ci, int li, 
              double lq, im_event *ek, 
              int np, double rates[], 
              im_ginfo *gk, 
              im_bfupdate *pbf,
              im_bfupdate *qbf,
              int *ki,
              int z1,
              char *is_coalesced1,
              int *pc1)
{
  int pi;
  int pj;
  int pk;
  int pl;
  int edge;
  int n1;
  int n2;
  struct genealogy *G;
  struct edge *gtree; 
  
  G = &C[ci]->G[li];
  gtree = G->gtree;

  switch (ek->type)
    {
    case 'a': /* Active lineage coalesces. */
      assert (np == 1 || np == 2);
      assert (rates[0] > 0.0);
      assert (ek->ei == z1);
      assert (*is_coalesced1 == 'F');
      lq += log (rates[0]);
      if (L[li].model == INFINITESITES)
        {
          assert (pbf->nz1[*pc1] > 0);
          lq -= log ((double) pbf->nz1[*pc1]);
        }
      else if (L[li].model == HKY)
        {
          assert (pbf->n[*pc1] > 0);
          lq += pbf->likes1[pbf->likei1]; /* This has been normalized. */
        }
      else if (L[li].model == STEPWISE)
        {
          assert (pbf->n[*pc1] > 0);
          lq += pbf->likes1[pbf->likei1]; /* This has been normalized. */
          lq += pbf->likes11[pbf->likei1][pbf->likei11]; /* This has been normalized. */
        }
      else
        {
          assert (0);
        }
      *is_coalesced1 = 'T';
      break;
    case 'p': /* Active lineage [[z1]] migrates.  */
      pj = ek->pj;
      assert (ek->ei == z1);
      assert (gk->mms[*pc1][pj] > 0.0);
      lq += log (gk->mms[*pc1][pj]);
      *pc1 = pj;
      break;
    case 'e':
      assert (pbf->is_finiterootbranch == 'T');
      break;
    case 'c':
      pi = ek->pj;
      edge = ek->ei;
      assert (*is_coalesced1 == 'F');
      if (*is_coalesced1 == 'F')
        {
          if (L[li].model == INFINITESITES)
            {
              /* NOTE: We do not update sbfold.l!!!           */
              /* We could have called IMA_sbf_c with n and l. */
              pbf->n[pi]--;
              IMA_sbf_c (ci, li, edge, pi, pbf->seqz1,
                         pbf->nz1, pbf->lz1, pbf->a);
            }
          else if (L[li].model == HKY || L[li].model == STEPWISE)
            {
              IMA_sbf_c (ci, li, edge, pi, NULL, 
                         pbf->n, pbf->l, pbf->a);
            }
        }
      break;
    case 'm':
      edge = ek->ei;
      pi = ek->pi;
      pj = ek->pj;
      assert (*is_coalesced1 == 'F');
      if (L[li].model == INFINITESITES)
        {
          pbf->n[pi]--;
          pbf->n[pj]++;
          IMA_sbf_m (ci, li, edge, pi, pj, 
                     pbf->nz1, pbf->lz1, pbf->a);
        }
      else if (L[li].model == HKY || L[li].model == STEPWISE)
        {
          IMA_sbf_m (ci, li, edge, pi, pj, 
                     pbf->n, pbf->l, pbf->a);
        }
      break;
    case 's':
      assert (ek->p == *ki);
      (*ki)++;
      pj = C[ci]->droppops[*ki][IM_EDGE_LEFT];
      pk = C[ci]->droppops[*ki][IM_EDGE_RIGHT];
      pl = C[ci]->addpop[*ki];
      pbf->n[pl] = pbf->n[pj] + pbf->n[pk];
      n1 = pbf->n[pj];
      n2 = pbf->n[pk];
      if (n1 > 0)
        {
          memcpy (pbf->l[pl], pbf->l[pj], n1 * sizeof (int));
        }
      if (n2 > 0)
        {
          memcpy (&pbf->l[pl][n1], pbf->l[pk], n2 * sizeof (int));
        }
      pbf->n[pj] = 0;
      pbf->n[pk] = 0;
      assert (*is_coalesced1 == 'F');
      *pc1 = saC.popndown[*pc1][*ki];
      assert (L[li].model == INFINITESITES 
              || L[li].model == HKY
              || L[li].model == STEPWISE);
      /* Change the number and sets as well. */
      n1 = pbf->nz1[pj];
      n2 = pbf->nz1[pk];
      pbf->nz1[pl] = n1 + n2;
      if (n1 > 0)
        {
          memcpy (pbf->lz1[pl], pbf->lz1[pj], n1 * sizeof (int));
        }
      if (n2 > 0)
        {
          memcpy (&pbf->lz1[pl][n1], pbf->lz1[pk], n2 * sizeof (int));
        }
      pbf->nz1[pj] = 0;
      pbf->nz1[pk] = 0;
      break;
    }
  return lq;
}
int 
bfmove_ratesI (double (*rates)[], 
               im_ginfo *gk, 
               im_ginfo *gkz1, 
               im_bfupdate *pbf,
               int ki,
               int pc1)
{
  int np;
  int li;

  np = 0;
  li = pbf->m1->li;

  if (L[li].model == INFINITESITES)
    {
      if (gkz1->thetas[pc1] == 0.0)
        {
          assert (pbf->nz1[pc1] == 0);
          (*rates)[np++] = 0.0;
        }
      else
        {
          assert (gkz1->thetas[pc1] > 0.0);
          (*rates)[np++] = 2 * pbf->nz1[pc1] / gkz1->thetas[pc1];
        }
    }
  else if (L[li].model == HKY || L[li].model == STEPWISE)
    {
      if (gk->thetas[pc1] == 0.0)
        {
          assert (pbf->n[pc1] == 0);
          (*rates)[np++] = 0.0;
        }
      else
        {
          assert (gk->thetas[pc1] > 0.0);
          (*rates)[np++] = 2 * pbf->n[pc1] / gk->thetas[pc1];
        }
    }
  else
    {
      assert (0);
    }
  if (ki < lastperiodnumber)
    {
      assert (gk->ms[ki][pc1] > 0.0);
      (*rates)[np++] = gk->ms[ki][pc1];
    }

  return np;
}

void
bfmove_ratesIII (double (*rates)[], im_ginfo *gk, im_bfupdate *pbf,
                 int ki,
                 int pc1, int pc3)
{
  if (ki < lastperiodnumber)
    {
      (*rates)[0] = gk->ms[ki][pc1];
      (*rates)[1] = gk->ms[ki][pc3];
    }
  else
    {
      (*rates)[0] = 0.0;
      (*rates)[1] = 0.0;
    }
  (*rates)[2] = 0.0; 
  assert (pbf->canz31 == 'T');
  if (pc1 == pc3) 
    {
      if (gk->thetas[pc3] == 0.0)
        {
          /* FIXME: is this a rightful max? */
          (*rates)[2] = 2 * (2 - 1) / itheta[pc3].pr.max; 
        }
      else
        {
          (*rates)[2] = 2 * (2 - 1) / gk->thetas[pc3];
        }
    }
  return; 
}
void 
IMA_genealogy_detach (int ci, int li, im_bfupdate *pbf,
                      int z1)
{
  struct genealogy *G;
  struct edge *gtree; 
  double dnt_genealogy;
  int si;
  int down1;
  int downdown1;
  double sisupt;
  double bessiv;
  int d;
  double sislen;

  
  G = &C[ci]->G[li];
  gtree = G->gtree;

  pbf->is_finiterootbranch = 'F';
  dnt_genealogy = IMA_genealogy_time (ci, li, G->root);
  assert (!gsl_fcmp(dnt_genealogy, G->roottime, 1e-7));
  assert (pbf->m1->edgeid == z1);
  down1 = gtree[z1].down;

  
  if (L[li].model == INFINITESITES)
    {
      memcpy (pbf->seqz1, gtree[down1].seq, L[li].numsites * sizeof (int));
    }
  else if (L[li].model == HKY || L[li].model == STEPWISE)
    {
      /* No code. */
    }
  else
    {
      assert (0);
    }

  pbf->sis1 = IMA_genealogy_sister (ci, li, z1);
  assert (!(pbf->sis1 < 0));

  sisupt = IMA_genealogy_time (ci, li, pbf->sis1);

  downdown1 = gtree[down1].down;
  IMA_genealogy_save_edge (ci, li, z1);
  IMA_genealogy_save_edge (ci, li, pbf->sis1);
  IMA_genealogy_save_edge (ci, li, down1);
  if (!(downdown1 < 0))
    {
      IMA_genealogy_save_edge (ci, li, downdown1);
    }
  pbf->is_finiterootbranch = 'F';
  if (downdown1 == -1)
    {
      assert (G->root == down1);
      G->root = pbf->sis1;
      gtree[pbf->sis1].down = -1;
      pbf->is_finiterootbranch = 'T';
      /**************************************************/
      /* We keep gtree[sis1].time                       */
      /* We keep migration events along the sister.     */
      /* We keep coalescent event at the original root  */
      /* by naming it ``end.''                          */
      /* We do not change roottime.                     */
      /**************************************************/
    }
  else
    {
      IMA_genealogy_join (gtree, pbf->sis1, down1, downdown1);
    }

  if (L[li].model == INFINITESITES)
    {
      /* Move mutations on down branch to sister. */
      for (si = 0; si < L[li].numsites; si++)
        {
          if (G->mut[si] == down1)
            {
              G->mut[si] = pbf->sis1;
            }
        }
    }
  else if (L[li].model == HKY)
    {
      /* No code. */
    }
  else if (L[li].model == STEPWISE)
    {
      pbf->A1 = gtree[down1].A[0];
/* 
      if (downdown1 == -1)
            {
                }
          else
            {
                  d = gtree[pbf->sis1].A[0] - gtree[downdown1].A[0];
                  sislen = gtree[pbf->sis1].time - sisupt;
                  bessiv = bessi (d, sislen * G->uvals[0]);
                  if (!(bessiv > 0.0))
                        {
                          assert (0); 
                        }
                  else
                        {
                      gtree[pbf->sis1].dlikeA[0] = (-(sislen * G->uvals[0]) + log (bessiv));
                        }
                }
*/
    }
  else
    {
      assert (0);
    }

  gtree[z1].exist = 'F';
  gtree[down1].exist = 'F';

  return;
}
int
IMA_genealogy_findIntSeq (int ci, int li)
{
  int si;
  int ei;
  int edge;
  int sis;
  int down;
  int ledge;
  int redge;
  int nmut;
  int ngenes;
  int ngnodes;
  char flag;
  struct genealogy *G;
  struct edge *gtree; 
  int accp;

  assert (L[li].model == INFINITESITES);
  G = &C[ci]->G[li];
  gtree = G->gtree;
  ngenes = L[li].numgenes;
  ngnodes = L[li].numlines;
  /* Reset seq and mut of a full genealogy */
  for (si = 0; si < L[li].numsites; si++)
    {
      for (ei = ngenes; ei < ngnodes; ei++)
        {
          gtree[ei].seq[si] = -1; 
        }
      G->mut[si] = -1;
    }
  for (si = 0; si < L[li].numsites; si++)
    {
      for (ei = 0; ei < ngenes; ei++)
        {
          edge = ei;
          flag = 'T';
          while (flag == 'T')
            {
              flag = 'F';
              down = gtree[edge].down;
              if (down == -1)
                {
                  assert (G->root == edge);
                  break;
                }
              if ((sis = gtree[down].up[IM_EDGE_LEFT]) == edge)
                {
                  sis = gtree[down].up[IM_EDGE_RIGHT];
                }
              assert (gtree[edge].seq[si] == 0
                      || gtree[edge].seq[si] == 1
                      || gtree[edge].seq[si] == 2);
              assert (!(gtree[down].seq[si] == -1 && gtree[edge].seq[si] == -1));
              if (gtree[down].seq[si] != gtree[edge].seq[si])
                {
                  if (gtree[edge].seq[si] == 2)
                    {
                      if (gtree[sis].seq[si] != -1)
                        {
                          gtree[down].seq[si] = gtree[sis].seq[si];
                        }
                      else
                        {
                          gtree[down].seq[si] = 2;
                        }
                    }
                  else if ((gtree[sis].seq[si] != -1 && gtree[sis].seq[si] != 2)
                           && (gtree[sis].seq[si] != gtree[edge].seq[si]))
                    {
                      gtree[down].seq[si] = 2;
                      /* assert (G->mut[si] == -1 || G->mut[si] == down); */
                      G->mut[si] = down;
                    }
                  else
                    {
                      gtree[down].seq[si] = gtree[edge].seq[si];
                    }
                  flag = 'T';
                }
              else if (gtree[down].seq[si] == 2)
                {
                  assert (gtree[edge].seq[si] == 2);
                  if (gtree[sis].seq[si] != -1 && gtree[sis].seq[si] != 2)
                    {
                      gtree[down].seq[si] = gtree[sis].seq[si];
                      flag = 'T';
                    }
                }
              edge = down;
            }
        }
    }
  /* Find a single edge with a mutation event at site [[si]]. */
  /* All edges except one must be labelled either 0 or 1. The only
   * edge is labelled 2. The edge is can be found at G->mut[si]. */
  accp = 1;
  for (si = 0; si < L[li].numsites; si++)
    {
      nmut = 0; 
      for (ei = ngenes; ei < ngnodes; ei++)
        {
          assert (gtree[ei].seq[si] == 0
                  || gtree[ei].seq[si] == 1
                  || gtree[ei].seq[si] == 2);
          if (gtree[ei].seq[si] == 2)
            {
              nmut++;
              ledge = gtree[ei].up[IM_EDGE_LEFT];
              redge = gtree[ei].up[IM_EDGE_RIGHT];
              /* assert (gtree[ledge].seq[si] != gtree[redge].seq[si]); */
              if (gtree[ledge].seq[si] == gtree[redge].seq[si])
                {
                  nmut = 999;
                  break;
                }
              if (ei == G->root)
                {
                  gtree[ei].seq[si] = 0;
                  if (gtree[ledge].seq[si] == 0)
                    {
                      G->mut[si] = redge;
                    }
                  else
                    {
                      G->mut[si] = ledge;
                    }
                }
              else
                {
                  down = gtree[ei].down;
                  assert (!(down < 0));
                  assert (down < ngnodes);
                  if (gtree[ledge].seq[si] == gtree[down].seq[si])
                    {
                      G->mut[si] = redge;
                      gtree[ei].seq[si] = gtree[ledge].seq[si];
                    }
                  else
                    {
                      G->mut[si] = ledge;
                      gtree[ei].seq[si] = gtree[redge].seq[si];
                    }
                }
            }
        }
      /* assert (nmut == 1); */
      if (nmut != 1)
        {
          accp = 0;
          break;
        }
    }
  return accp;
}

int 
IMA_sbf_c (int ci, int li, int edge, int pi, 
           int *seqz,
           int *n, int **l, int *a)
{
  int j;
  int ei;
  char can_coalesce;
  int si;
  int ni;
  int nl;
  struct genealogy *G;
  struct edge *gtree;
  int ledge;
  int redge;

  G = &C[ci]->G[li];
  gtree = G->gtree;
  ledge = gtree[edge].up[IM_EDGE_LEFT];
  redge = gtree[edge].up[IM_EDGE_RIGHT];

  /* Remove edges (both left and right) that coalesce together if there is. It
   * is possible that they are absent in [[l]]. */
  nl = n[pi];
  j = 0;
  for (ei = 0; ei < nl; ei++)
    {
      if (l[pi][ei] == ledge || l[pi][ei] == redge)
        {
          /* assert ledge could have coalesced with z. */
        }
      else
        {
          a[j] = l[pi][ei];
          j++;
        }
    }
  assert (j == nl || j == nl - 1 || j == nl - 2);
  if (j < nl)
    {
      if (j > 0)
        {
          memcpy (l[pi], a, j * sizeof (int));
        }
      n[pi] = j;
    }

  /* Add edge if it can coalesce with active lineage z. */
  can_coalesce = 'T';
  if (seqz != NULL)
    {
      for (si = 0; si < L[li].numsites; si++)
        {
          if (G->mut[si] == edge)
            {
              continue;
            }
          if (gtree[edge].seq[si] != seqz[si])
            {
              can_coalesce = 'F';
            }
        }
    }
  if (can_coalesce == 'T') /* edge and z can coalesce together? */
    {
      ni = n[pi];
      l[pi][ni] = edge;
      n[pi]++;
    }
  return 0;
} 

int 
IMA_sbf_m (int ci, int li, int edge, int pi, int pj,
           int *n, int **l, int *a)
{
  int j;
  int ei;
  int ni;
  int nl;

  assert (pi != pj);

  /* Remove if edge exists in pi. */
  nl = n[pi];
  j = 0;
  for (ei = 0; ei < nl; ei++)
    {
      if (l[pi][ei] == edge)
        {
          /* assert ledge could have coalesced with z1. */
        }
      else
        {
          a[j] = l[pi][ei];
          j++;
        }
    }
  assert (j == nl || j == nl - 1);
  if (j < nl)
    {
      /* Add the edge to population pj. */
      if (j > 0)
        {
          memcpy (l[pi], a, j * sizeof (int));
        }
      n[pi] = j;
      ni = n[pj];
      for (ei = 0; ei < ni; ei++)
        {
          assert (l[pj][ei] != edge);
        }
      l[pj][ni] = edge;
      n[pj]++;
    }
  else if (j == nl)
    {
      /* No Code! */
    }
  else
    {
      assert (0);
    }
  return 0;
}

int
IMA_sbf_nsum (im_bfupdate *pbf, int ci, int li, 
              int pc1, int pc3)
{
  int nsum;
  int pi;

  /* Make sure that we have one inactive lineage. */
  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      if (pi == pc3)
        {
          assert (pbf->n[pc3] == 1);
        }
      else
        {
          assert (pbf->n[pi] == 0);
        }
    }
  pbf->n[pc1]++;

  nsum = 0;
  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      nsum += pbf->n[pi];
    }
  assert (nsum == 2);

  return nsum;
}
void 
IMA_sbf_lineages (im_bfupdate *pbf, int ci, int li, int z1)
{
  int n;
  int p;
  int np;
  int ni;
  int si;
  int ei;
  int edge;
  struct edge *gtree;
  struct genealogy *G;
  int ngenes;
  char can_coalesce;

  for (p = 0; p < numpopsizeparams; p++)
    {
      pbf->n[p] = 0;
      pbf->nz1[p] = 0; /* Only for infinite-sites model */
    }
  ngenes = L[li].numgenes;
  G = &C[ci]->G[li];
  gtree = G->gtree;
  for (ei = 0; ei < ngenes; ei++)
    {
      /* Exclude a gene at an active lineage. */
      if (ei == z1) 
        {
          continue; 
        }
      assert (!(gtree[ei].pop < 0));
      assert (gtree[ei].pop < 2 * MAXPOPS - 1);
      p = gtree[ei].pop;
      /* Add gene [[ei]] to population [[p]]. */
      n = pbf->n[p];
      pbf->l[p][n] = ei;
      pbf->n[p]++;
    }

  /* Number of lineages in populations for active lineages to be able to
   * coalesce. */
  if (L[li].model == INFINITESITES)
    {
      for (p = 0; p < numpopsizeparams; p++)
        {
          np = pbf->n[p];
          for (ni = 0; ni < np; ni++)
            {
              edge = pbf->l[p][ni];
              can_coalesce = 'T';
              for (si = 0; si < L[li].numsites; si++)
                {
                  if (G->mut[si] == edge)
                    {
                      continue;
                    }
                  if (gtree[edge].seq[si] != pbf->seqz1[si])
                    {
                      can_coalesce = 'F';
                      break;
                    }
                }
              if (can_coalesce == 'T') /* If edge and z1 can coalesce together? */
                {
                  n = pbf->nz1[p];
                  pbf->lz1[p][n] = edge;
                  pbf->nz1[p]++;
                }
            }
        }
    }
  return;
}
void 
IMA_intervals_collectsbfnew (int ci, int li, im_edgemiginfo *newm3)
{
  double ct;
  int pc;

  im_event *ek;
  int mi;
  int mj;
  int ei;
  int si;
  int ki;
  int pi;
  int z3;
  struct edge *gtree;
  struct genealogy *G;
  double t;
  double roottime;
  char m3type;
  int down1;

  G = &C[ci]->G[li];
  gtree = C[ci]->G[li].gtree;
  if (sbfnew.sis1 == G->root)
    {
      sbfnew.is_finiterootbranch = 'T';
    }
  else
    {
      sbfnew.is_finiterootbranch = 'F';
    }
  
  down1 = gtree[sbfnew.m1->edgeid].down;
  if (L[li].model == INFINITESITES)
    {
      assert (L[li].numsites > 0); /* FIXME: Crash! */
      memcpy (sbfnew.seqz1, gtree[down1].seq, L[li].numsites * sizeof (int)); 
    }
  else if (L[li].model == HKY)
    {
      assert (L[li].numsites > 0); /* FIXME: Crash! */
      /* No code. */
    }
  else if (L[li].model == STEPWISE)
    {
      /* No code. */
    }
  else
    {
      assert (0);
    }

  ct = sbfnew.m1->dnt;
  pc = sbfnew.m1->fpop;

  /***************************************************************************/
  /*                             sbfold.m3                                   */
  /***************************************************************************/
  /* Fill sbfold.m3 when detaching an external root-connecting branch and    */
  /* attaching it to an inactive lineage.                                    */
  if (newm3 == NULL && sbfold.is_finiterootbranch == 'T')
    {
      ei = G->root;
      roottime = IMA_genealogy_time (ci, li, ei);

      sbfold.m3 = &saosms[2 * li];
      sbfold.m3->li = li;
      sbfold.m3->edgeid = ei;
      if (ct < roottime)
        {
          m3type = 'A';
          sbfold.m3->upt = roottime;
          sbfold.m3->pop = gtree[ei].pop;
        }
      else
        {
          m3type = 'B'; 
          sbfold.m3->upt = ct;
          sbfold.m3->pop = pc;
        }
      assert (ct < gtree[ei].time);
      assert (gtree[ei].time < TIMEMAX - 1.0);
      sbfold.m3->dnt = gtree[ei].time;

      mi = 0;
      mj = 0;
      while (gtree[ei].mig[mi].mt > -0.5)
        {
          if (gtree[ei].mig[mi].mt < ct)
            { 
              mi++;
              continue;
            }
          sbfold.m3->mig[mj] = gtree[ei].mig[mi];
          mi++;
          mj++;
        }
      sbfold.m3->mig[mj].mt = -1.0;
      sbfold.m3->mpall = mj;
      if (mj == 0)
        {
          pi = sbfold.m3->pop;
        }
      else
        {
          pi = sbfold.m3->mig[mj - 1].mp;
        }
      ki = findperiodoftime (ci, sbfold.m3->dnt);
      sbfold.m3->fpop = saC.popndown[pi][ki];
      if (ct < roottime)
        {
          /* No Code! */
        }
      else
        {
          assert (sbfnew.m1->dnt == sbfold.m3->upt);
        }
      assert (sbfnew.m1->dnt < sbfold.m3->dnt);

      /* Why are we doing this here? */
      /* Set sbfnew's canz31? */
      z3 = G->root;
      sbfnew.canz31 = 'T';
      if (L[li].model == INFINITESITES)
        {
          for (si = 0; si < L[li].numsites; si++)
            {
              if (G->mut[si] == z3)
                {
                  continue;
                }
              if (gtree[z3].seq[si] != sbfnew.seqz1[si])
                {
                  sbfnew.canz31 = 'F';
                }
            }
        }
      else if (L[li].model == HKY || L[li].model == STEPWISE)
        {
          /* No code. */
        }

      assert (sbfnew.canz31 == 'T');
    }
  else
    {
      sbfold.m3 = NULL;
    }

  /***************************************************************************/
  /*                                 sbfnew                                  */
  /***************************************************************************/

  /* Copy e1 from sbfold's to sbfnew's. */
  if (sbfnew.nmax1 < sbfold.n1)
    {
      sbfnew.nmax1 = sbfold.n1;
      sbfnew.e1 = realloc (sbfnew.e1, sbfnew.nmax1 * sizeof (im_event));
    }
  assert (sbfold.n1 > 0);
  memcpy (sbfnew.e1, sbfold.e1, sbfold.n1 * sizeof (im_event));
  sbfnew.n1 = sbfold.n1;

  /* Case I */
  if (newm3 != NULL)
    {
      /* sbfnew is lengthened by newm3. */
      /* if there is end event, we remove it */
      /* we can add ``end'' event. */
      if (sbfnew.e1[sbfnew.n1 - 1].type == 'e')
        {
          sbfnew.n1--;
        }
       
      if (sbfnew.nmax1 < sbfnew.n1 + newm3->mpall + numsplittimes + 1)
        {
          sbfnew.nmax1 = sbfnew.n1 + newm3->mpall + numsplittimes + 1;
          sbfnew.e1 = realloc (sbfnew.e1, sbfnew.nmax1 * sizeof (im_event));
        }
      /* ek = &sbfnew.e1[sbfnew.n1]; */
      ek = sbfnew.e1 + sbfnew.n1;
      mi = 0;
      while (newm3->mig[mi].mt > -0.5)
        {
          t = newm3->mig[mi].mt;
          ek->type = 'm';
          ek->t = t;
          ek->p = findperiodoftime (ci, t);
          if (mi == 0)
            {
              ek->pi = newm3->pop;
            }
          else
            {
              ek->pi = newm3->mig[mi - 1].mp;
            }
          assert (!(ek->pi < 0));
          assert (ek->pi < numpopsizeparams);
          ek->pi = saC.popndown[ek->pi][ek->p];
          ek->pj = newm3->mig[mi].mp;
          assert (!(ek->pi < 0));
          assert (ek->pi < numpopsizeparams);
          assert (!(ek->pj < 0));
          assert (ek->pj < numpopsizeparams);
          ek->ei = newm3->edgeid;
          ek++;
          mi++;
        }
      assert (newm3->mpall == mi); 
      /* end event */
      t = newm3->dnt;
      ek->type = 'e';
      ek->t = t;
      ek->p = findperiodoftime (ci, t);
      if (mi == 0)
        {
          ek->pi = newm3->pop;
        }
      else
        {
          ek->pi = newm3->mig[mi - 1].mp;
        }
      assert (!(ek->pi < 0));
      assert (ek->pi < numpopsizeparams);
      ek->pj = saC.popndown[ek->pi][ek->p];
      assert (!(ek->pj < 0));
      assert (ek->pj < numpopsizeparams);
      ek->ei = newm3->edgeid;
      sbfnew.n1 += newm3->mpall + 1;

      /* We need to sort the events. */
      ek = sbfnew.e1 + sbfnew.n1;
      for (ki = 0; ki < numsplittimes; ki++)
        {
          t = C[ci]->tvals[ki];
          if (newm3->upt < t && t < newm3->dnt)
            {
              ek->type = 's'; 
              ek->t = t;
              ek->p = ki;    /* Split time happens right before the period line. */
              ek->pi = ki;   /* from 0 (pi) to 1 (pj) for the first split time.  */
              ek->pj = ki + 1;
              ek->ei = -1;
              ek++;
              sbfnew.n1++;
            }
        }
      qsort (sbfnew.e1, (size_t) sbfnew.n1, sizeof (im_event), im_event_cmp);
    }
  /* Case II */
  else if (newm3 == NULL && sbfold.is_finiterootbranch == 'T')
    {
      assert (sbfold.m3 != NULL);
      /* sbfnew is shortened by oldm3. */
      /* if there is end event, we remove it */
      /* we may add ``end'' event. */
      if (sbfnew.e1[sbfnew.n1 - 1].type == 'e')
        {
          sbfnew.n1--;
        }
      /* Remove events below oldm3 */
      if (m3type == 'A')
        {
           while (sbfnew.e1[sbfnew.n1 - 1].t > sbfold.m3->upt) 
            {
              if (sbfnew.e1[sbfnew.n1 - 1].type == 'c'
                  && sbfnew.e1[sbfnew.n1 - 1].ei == sbfold.m3->edgeid)
                {
                  break;
                }
              sbfnew.n1--;
            }
        }
      else if (m3type == 'B')
        {
          while (sbfnew.e1[sbfnew.n1 - 1].t > sbfold.m3->upt) 
            {
              assert (!(sbfnew.e1[sbfnew.n1 - 1].type == 'c'
                        && sbfnew.e1[sbfnew.n1 - 1].ei == sbfold.m3->edgeid));
              sbfnew.n1--;
            }
          ek = sbfnew.e1 + sbfnew.n1;
          ek->ei = sbfold.m3->edgeid;
          ek->type = 'e';
          ek->t = sbfold.m3->upt;
          ek->p = findperiodoftime (ci, ek->t);
          mi = 0;
          t = gtree[ek->ei].mig[mi].mt;
          if (t > -0.5 && t < ek->t)
            {
              mi++;
              t = gtree[ek->ei].mig[mi].mt;
            }
          if (mi == 0)
            {
              ek->pi = gtree[ek->ei].pop;
            }
          else
            {
              ek->pi = gtree[ek->ei].mig[mi - 1].mp;
            }
          ek->pj = sbfold.m3->pop;
          sbfnew.n1++;
        }
      else
        {
          assert (0);
        }
      /* No need to sort */
    }
  return;
}

double 
bfq (im_bfupdate *pbf, im_bfupdate *qbf, 
     int ci, int li, 
     im_ginfo *gk,
     im_ginfo *gkz1)
{
  double lq;
  char is_coalesced1;
  char is_coalesced3;
  int pc1;
  double prevt;
  double dt;
  int ki;
  int ti;
  int ti3;
  int pj;
  int np;
  int mi;
  double totalrate;
  im_event *ek;
  double lasteventtime;
  double rates[2 * MAXPOPS];
  int nsum;
  double t;
  struct edge *gtree;
  struct genealogy *G;
  int ngenes;
  double newct;
  int pc3;
  int ri;
  int nr;
  int z1;
  int z3;
  double r;
  char eventAadded;

  G = &C[ci]->G[li];
  gtree = G->gtree;
  ngenes = L[li].numgenes;

  /***************************************************************************/
  /* Prepare [[e2]] and possibly [[e3]].                                     */
  /***************************************************************************/
  /* Type 'a': lineage 1 coalesces. 
   *      'p': lineage 1 migrates.
   *      'r': lineage 3 migrates. */
  lasteventtime = pbf->e1[pbf->n1 - 1].t;

  if (pbf->nmax2 < pbf->n1)
    {
      pbf->nmax2 = pbf->n1;
      pbf->e2 = realloc (pbf->e2, pbf->nmax2 * sizeof (im_event));
    }
  pbf->e2 = memcpy (pbf->e2, pbf->e1, pbf->n1 * sizeof (im_event));
  pbf->n2 = pbf->n1;
  if (pbf->nmax2 < pbf->n2 + qbf->m1->mpall + 1)
    {
      pbf->nmax2 = pbf->n2 + qbf->m1->mpall + 1;
      pbf->e2 = realloc (pbf->e2, pbf->nmax2 * sizeof (im_event));
    }
  ek = pbf->e2 + pbf->n2;
  for (mi = 0; mi < qbf->m1->mpall; mi++)
    {
      t = qbf->m1->mig[mi].mt;
      assert (t > -0.5);
      if (t > lasteventtime)
        {
          break;
        }
      ek->type = 'p';
      ek->t = t;
      ek->p = findperiodoftime (ci, t);
      if (mi == 0)
        {
          ek->pi = qbf->m1->pop;
        }
      else
        {
          ek->pi = qbf->m1->mig[mi - 1].mp;
        }
      ek->pi = saC.popndown[ek->pi][ek->p];
      ek->pj = qbf->m1->mig[mi].mp;
      ek->ei = qbf->m1->edgeid;
      ek++;
      pbf->n2++;
    }
  if (qbf->m1->dnt < lasteventtime)
    {
      /* Active lineage [[z1]] coalesces. */
      ek->type = 'a';
      ek->t = qbf->m1->dnt;
      ek->p = -1;
      ek->pi = -1;
      ek->pj = qbf->m1->fpop;
      ek->ei = qbf->m1->edgeid;
      ek++;
      pbf->n2++;
      eventAadded = 'y';
    }
  else
    {
      eventAadded = 'n';
    }
  qsort (pbf->e2, (size_t) pbf->n2, sizeof (im_event), im_event_cmp);

  /* old3 = two or three active lineages   */
  /* qbf->m1->mpall and a coalescent       */
  /* qbf->m3->mpall and a coalescent       */
  pbf->n3 = 0;
  if (qbf->m3 != NULL)
    {
      if (pbf->nmax3 < qbf->m1->mpall + 1)
        {
          pbf->nmax3 = qbf->m1->mpall + 1;
          pbf->e3 = realloc (pbf->e3, pbf->nmax3 * sizeof (im_event));
        }
      ek = pbf->e3;
      for (mi = 0; mi < qbf->m1->mpall; mi++)
        {
          t = qbf->m1->mig[mi].mt; 
          assert (t > -0.5);
          if (t < lasteventtime)
            {
              continue;
            }
          ek->type = 'p';
          ek->t = t;
          ek->p = findperiodoftime (ci, t);
          if (mi == 0)
            {
              ek->pi = qbf->m1->pop;
            }
          else
            {
              ek->pi = qbf->m1->mig[mi - 1].mp;
            }
          ek->pi = saC.popndown[ek->pi][ek->p];
          ek->pj = qbf->m1->mig[mi].mp;
          ek->ei = qbf->m1->edgeid;
          ek++;
          pbf->n3++;
        }
      assert (qbf->m1->dnt > lasteventtime);
      assert (eventAadded = 'n');
      ek->type = 'a';
      ek->t = qbf->m1->dnt;
      ek->p = -1;
      ek->pi = -1;
      ek->pj = qbf->m1->fpop;
      ek->ei = qbf->m1->edgeid;
      ek++;
      pbf->n3++;

      if (pbf->nmax3 < pbf->n3 + qbf->m3->mpall + 1)
        {
          pbf->nmax3 = pbf->n3 + qbf->m3->mpall + 1;
          pbf->e3 = realloc (pbf->e3, pbf->nmax3 * sizeof (im_event));
        }
      ek = pbf->e3 + pbf->n3;
      for (mi = 0; mi < qbf->m3->mpall; mi++)
        {
          t = qbf->m3->mig[mi].mt;
          assert (t > lasteventtime);
          assert (t > -0.5);
          ek->type = 'r';
          ek->t = t;
          ek->p = findperiodoftime (ci, t);
          if (mi == 0)
            {
              ek->pi = qbf->m3->pop;
            }
          else
            {
              ek->pi = qbf->m3->mig[mi - 1].mp;
            }
          ek->pi = saC.popndown[ek->pi][ek->p];
          ek->pj = qbf->m3->mig[mi].mp;
          ek->ei = qbf->m3->edgeid;
          ek++;
          pbf->n3++;
        }
      if (pbf->nmax3 < pbf->n3 + numsplittimes + 1)
        {
          pbf->nmax3 = pbf->n3 + numsplittimes + 1;
          pbf->e3 = realloc (pbf->e3, pbf->nmax3 * sizeof (im_event));
        }
      ek = pbf->e3 + pbf->n3;
      for (ki = 0; ki < numsplittimes; ki++)
        {
          t = C[ci]->tvals[ki];
          if (qbf->m3->upt < t && t < qbf->m3->dnt)
            {
              ek->type = 's'; 
              ek->t = t;
              ek->p = ki;        /* Split time happens right before the period line. */
              ek->pi = ki;       /* from 0 (pi) to 1 (pj) for the first split time.  */
              ek->pj = ki + 1;
              ek++;
              pbf->n3++;
            }
        }
      qsort (pbf->e3, (size_t) pbf->n3, sizeof (im_event), im_event_cmp);
      assert (pbf->e3[pbf->n3 - 1].type == 'a');
    }
  /***************************************************************************/
  /* Prepare [[e2]] and possibly [[e3]].                                     */
  /***************************************************************************/

  /* Add probabilities before the bottom of a partial genealogy */
  /* pbf:sbfold, and qbf:sbfnew.                                */
  z1 = qbf->m1->edgeid;
  assert (!(z1 < 0));
  is_coalesced1 = 'F';
  pc1 = qbf->m1->pop;
  assert (!(pc1 < 0));
  IMA_sbf_lineages (pbf, ci, li, z1);

  ki = 0;
  ti = 0;
  lq = 0.0; 
  prevt = 0.0;
  ek = pbf->e2;
  while (is_coalesced1 == 'F') 
    {
      dt = ek->t - prevt;
      np = bfmove_ratesI (&rates, gk, gkz1,
                          pbf, ki, 
                          pc1);
      totalrate = 0.0;
      for (pj = 0; pj < np; pj++)
        {
          assert (!(rates[pj] < 0.0));
          totalrate += rates[pj];
        }
      if (totalrate == 0.0)
        {
          /* assert (0); I do not know when this would happen. We have to fix
           * it. */
          assert (ek->type == 'c'
                  || ek->type == 'm'
                  || ek->type == 's'
                  || ek->type == 'e');
        }
      lq -= totalrate * dt;
      lq = bfmove_nextI (ci, li, lq, ek, np, rates, 
                         gk, pbf, qbf,
                         &ki,
                         z1,
                         &is_coalesced1,
                         &pc1);
      prevt = ek->t;
      ti++;
      if (ti < pbf->n2)
        {
          ek++;
        }
      else
        {
          break;
        }
    }

  /* Add probabilities after the bottom of a partial genealogy. */
  if (is_coalesced1 == 'F')
    {
      assert (pbf->m3 == NULL);
      assert (qbf->m3 != NULL);
      assert (pbf->n3 > 0);
      assert (pbf->n2 == ti);
      assert (pbf->e2[ti - 1].type != 's');
      assert (pbf->e2[ti - 1].type == 'e' || pbf->e2[ti - 1].type == 'c');
      assert (pbf->e2[ti - 1].t == prevt);
      assert (findperiodoftime (ci, prevt) == ki);

      is_coalesced3 = 'F';
      newct = prevt;
      pc3 = pbf->e2[ti - 1].pj;
      z3 = pbf->e2[ti - 1].ei;
      assert (z3 == qbf->m3->edgeid);
      assert (!(pc1 < 0));
      assert (!(pc3 < 0));
      assert (qbf->m3->upt == newct);
      assert (qbf->m3->li == li);
      assert (qbf->m3->edgeid == G->root);
      assert (qbf->m3->pop == pc3);

      /* Set canz31? */
      assert (is_coalesced1 == 'F');
      assert (pbf->canz31 == 'T');
      
      nsum = IMA_sbf_nsum (pbf, ci, li, pc1, pc3);
      assert (nsum == 2);

      ti3 = 0;
      ek = pbf->e3;
      while (nsum == 2)
        {
          dt = ek->t - prevt;

          bfmove_ratesIII (&rates, gk, pbf, ki, pc1, pc3);
          nr = 3;
          totalrate = 0.0;
          for (ri = 0; ri < nr; ri++)
            {
              totalrate += rates[ri];
            }
          assert (totalrate > 0.0);
          assert (dt > 0.0);
          lq -= totalrate * dt;
          switch (ek->type)
            {
            case 's': /* split event. */
              assert (ek->pi == ki);
              ki++;
              assert (ek->pj == ki);
              pc1 = saC.popndown[pc1][ki];
              pc3 = saC.popndown[pc3][ki];
              break;
            case 'p': /* one left active lineage migrates. */
              pj = ek->pj;
              assert (ek->ei == z1);      /* lineage 1 coalesces. */
              assert (rates[0] > 0.0);
              assert (gk->mms[pc1][pj] > 0.0);
              assert (gk->ms[ki][pc1] > 0.0);
              r = gk->mms[pc1][pj];
              lq += log (r);
              pc1 = pj;
              break;
            case 'r': /* Lineage 3 migrates. */
              pj = ek->pj;
              assert (ek->ei == z3);
              assert (rates[1] > 0.0); 
              assert (gk->mms[pc3][pj] > 0.0);
              assert (gk->ms[ki][pc3] > 0.0);
              lq += log (gk->mms[pc3][pj]);
              pc3 = pj;
              break;
            case 'a': /* Two lineages coalesce. */
              is_coalesced3 = 'T';
              assert (ek->ei == z1);
              assert (is_coalesced1 == 'F');
              is_coalesced1 = 'T';
              assert (rates[2] > 0.0); 
              r = rates[2];
              lq += log (r);
              lq = bfmove_selectA1 (ci, li, qbf, lq);
              nsum = 1;
              break;
            default:
              assert (0);
              break;
            }
          prevt = ek->t;
          ti3++;
          ek++;
        }
      assert (ti3 == pbf->n3);
    }

  return lq;
}

double 
IMA_choose_migtopop (int *topop, int ci, int li, int ki, int pi, im_ginfo *gk)
{
  double pr[MAXPOPS];
  int npminus1;
  int ni;
  int pj;
  int i;
  double total;

  total = 0.0;
  npminus1 = npops - ki - 1;
  for (ni = 0; ni < npminus1; ni++)
    {
      pj = saC.popnmig[ki][pi][ni];
      assert (gk->mms[pi][pj] > 0.0);
      pr[ni] = gk->mms[pi][pj];
      total += pr[ni];
    }
  assert (gk->ms[ki][pi] == total);

  i = ran_discretea (npminus1, pr);
  *topop = saC.popnmig[ki][pi][i];
  return pr[i] / total;
}

void
IMA_edgemiginfo_copymig (im_edgemiginfo *em, int ci, int li, int ei)
{
  int mi;
  int down;
  struct edge *gtree;

  gtree = C[ci]->G[li].gtree;
  down = gtree[ei].down;
  assert (!(down < 0));
  em->li = li;
  em->edgeid = ei;
  em->upt = IMA_genealogy_time (ci, li, ei);
  em->dnt = gtree[ei].time;
  em->pop = gtree[ei].pop;
  em->fpop = gtree[down].pop;
  mi = 0;
  while (gtree[ei].mig[mi].mt > -0.5)
    {
      em->mig[mi] = gtree[ei].mig[mi];
      mi++;
    }
  em->mig[mi].mt = -1;
  em->mpall = mi;

  return;
}
void
IMA_edgemiginfo_copypartmig (im_edgemiginfo *em, 
                             int ci, int li, int ei, double upt)
{
  int mi;
  int mj;
  int ki;
  int pi;
  int down;
  double ct;
  struct edge *gtree;

  gtree = C[ci]->G[li].gtree;
  down = gtree[ei].down;
  em->li = li;
  em->edgeid = ei;
  ct = IMA_genealogy_time (ci, li, ei);
  em->upt = upt;
  em->dnt = gtree[ei].time;
  assert (gtree[ei].time < TIMEMAX - 1.0);
  assert (ct < upt);
  assert (upt < gtree[ei].time);
  mi = 0;
  while (gtree[ei].mig[mi].mt > -0.5 && gtree[ei].mig[mi].mt < upt)
    {
      mi++;
    }
  if (mi == 0)
    {
      pi = gtree[ei].pop;
    }
  else
    {
      pi = gtree[ei].mig[mi - 1].mp;
    }
  ki = findperiodoftime (ci, upt);
  em->pop = saC.popndown[pi][ki];
  
  mj = 0;
  while (gtree[ei].mig[mi].mt > -0.5)
    {
      em->mig[mj] = gtree[ei].mig[mi];
      mi++;
      mj++;
    }
  em->mig[mj].mt = -1.0;
  em->mpall = mj;
  if (mj == 0)
    {
      pi = em->pop;
    }
  else
    {
      pi = em->mig[mj - 1].mp;
    }
  ki = findperiodoftime (ci, em->dnt);
  em->fpop = saC.popndown[pi][ki];
  return;
}

void
IMA_edgemiginfo_appendmig (im_edgemiginfo *em, int ci, int li, int ei)
{
  int mi;
  int mj;
  int down;
  struct edge *gtree;

  assert (!(em->edgeid < 0));
  gtree = C[ci]->G[li].gtree;
  down = gtree[ei].down;
  assert (!(down < 0));

  em->dnt = gtree[ei].time;
  em->fpop = gtree[down].pop;

  mi = 0;
  while (em->mig[mi].mt > -0.5)
    {
      mi++;
    }
  assert (em->mpall == mi);
  mj = 0;
  while (gtree[ei].mig[mj].mt > -0.5)
    {
      em->mig[mi] = gtree[ei].mig[mj];
      mi++;
      mj++;
    }
  em->mig[mi].mt = -1;
  em->mpall = mi;

  return;
}

void
IMA_genealogy_copymig (int ci, int li, int ei, im_edgemiginfo *em)
{
  int mi;
  struct edge *gtree;

  gtree = C[ci]->G[li].gtree;

  assert (em->li == li);
  assert (em->edgeid == ei);
  mi = 0;
  while (em->mig[mi].mt > -0.5)
    {
      checkmigt (mi, &gtree[ei]);
      gtree[ei].mig[mi] = em->mig[mi];
      mi++;
    }
  gtree[ei].mig[mi].mt = -1.0;
  assert (em->mpall == mi);

  return;
}

void
IMA_genealogy_appendmig (int ci, int li, int sis, im_edgemiginfo *em, 
                         int edge, int down)
{
  double sisupt;
  double sisdnt;
  int mj;
  int mi;
  struct genealogy *G;
  struct edge *gtree;

  assert (em->li == li);
  assert (em->edgeid == sis);

  G = &C[ci]->G[li];
  gtree = G->gtree;
  IMA_genealogy_save_edge (ci, li, sis);
  sisupt = IMA_genealogy_time (ci, li, sis);
  sisdnt = gtree[sis].time;
  assert (G->root == sis);
  G->root = down;
  assert (G->roottime < em->dnt);
  G->roottime = em->dnt;
  gtree[down].down = -1;
  gtree[sis].down = down;
  gtree[edge].down = down;
  gtree[down].up[IM_EDGE_LEFT] = sis;
  gtree[down].up[IM_EDGE_RIGHT] = edge;
  gtree[edge].time = em->dnt;
  gtree[sis].time = em->dnt;
  gtree[down].time = TIMEMAX; 
  /* A partial genealogy becomes a full genealogy although it may have another
   * active lineage that needs to be attached. */

  mi = 0;
  while (gtree[sis].mig[mi].mt > -0.5 && gtree[sis].mig[mi].mt < em->upt)
    {
      mi++;
    }
  mj = 0;
  while (em->mig[mj].mt > -0.5)
    {
      checkmigt (mi, &gtree[sis]);
      gtree[sis].mig[mi] = em->mig[mj];
      mi++;
      mj++;
    }
  gtree[sis].mig[mi].mt = -1.0;
  assert (em->mpall == mj);
  gtree[down].mig[0].mt = -1.0;
  gtree[down].pop = em->fpop;

  return;
}

void 
IMA_genealogy_join (im_edge *gtree, int up, int down, int downdown)
{
  int mi;
  int mj;

  gtree[up].down = downdown;
  if (gtree[downdown].up[IM_EDGE_LEFT] == down)
    {
      gtree[downdown].up[IM_EDGE_LEFT] = up;
    }
  else
    {
      gtree[downdown].up[IM_EDGE_RIGHT] = up;
    }

  mi = 0;
  while (gtree[up].mig[mi].mt > -0.5)
    {
      mi++;
    }
  mj = 0;
  while (gtree[down].mig[mj].mt > -0.5)
    {
      checkmig (mi + 1, &(gtree[up].mig), &(gtree[up].cmm));
      gtree[up].mig[mi] = gtree[down].mig[mj];
      mi++;
      mj++;
    }
  gtree[up].mig[mi].mt = -1;
  gtree[up].time = gtree[down].time;
  return;
}

void 
IMA_genealogy_absorbdown (im_edge *gtree, int up, int down)
{
  int mi;
  int mj;

  mi = 0;
  while (gtree[up].mig[mi].mt > -0.5)
    {
      mi++;
    }
  mj = 0;
  while (gtree[down].mig[mj].mt > -0.5)
    {
      checkmig (mi + 1, &(gtree[up].mig), &(gtree[up].cmm));
      gtree[up].mig[mi] = gtree[down].mig[mj];
      mi++;
      mj++;
    }
  gtree[up].mig[mi].mt = -1;
  gtree[up].time = gtree[down].time;
  return;
}

void 
IMA_genealogy_bisectsis (int ci, int li, 
                         int sis, int edge, int down, 
                         im_edgemiginfo *em, 
                         int *seqz)
{
  double sisupt;
  double sisdnt;
  double ct;
  int nowpop;
  int downdown;
  int mj;
  int mi;
  int ki;
  int si;
  struct genealogy *G;
  struct edge *gtree;

  /* I do not like that I have to keep calling these two lines. I wanted to
   * avoid using them. It turned out that I have to change some basic structures
   * */
  G = &C[ci]->G[li];
  gtree = G->gtree;
  downdown = gtree[sis].down;
  IMA_genealogy_save_edge (ci, li, sis);
  if (!(downdown < 0))
    {
      IMA_genealogy_save_edge (ci, li, downdown);
    }

  sisupt = IMA_genealogy_time (ci, li, sis);
  sisdnt = gtree[sis].time;
  ct = em->dnt;
  assert (sisupt < ct && ct < sisdnt);
  if (G->root == sis)
  {
    assert (downdown == -1);
    G->root = down;
    G->roottime = ct;
  }
  else
  {
    if (gtree[downdown].up[IM_EDGE_LEFT] == sis)
      {
        gtree[downdown].up[IM_EDGE_LEFT] = down;
      }
    else
      {
        gtree[downdown].up[IM_EDGE_RIGHT] = down;
      }
  }
  gtree[down].down = downdown;

  gtree[sis].down = down;
  gtree[edge].down = down;
  gtree[down].up[IM_EDGE_LEFT] = sis;
  gtree[down].up[IM_EDGE_RIGHT] = edge;

  /* Population label and migration events */
  gtree[edge].time = ct;
  gtree[sis].time = ct;
  gtree[down].time = sisdnt;
  mi = 0;
  while (gtree[sis].mig[mi].mt > -0.5 
         && gtree[sis].mig[mi].mt < ct)
    {
      mi++;
    }
  if (mi > 0)
    {
      nowpop = gtree[sis].mig[mi - 1].mp;
    }
  else
    {
      nowpop = gtree[sis].pop;
    }

  ki = findperiodoftime (ci, ct);
  nowpop = saC.popndown[nowpop][ki];
  assert (em->fpop == nowpop);
  gtree[down].pop = em->fpop;

  mj = 0;
  while (gtree[sis].mig[mj + mi].mt > -0.5)
    {
      checkmigt (mj, &gtree[down]);
      gtree[down].mig[mj] = gtree[sis].mig[mj + mi];
      assert (nowpop != gtree[down].mig[mj].mp);
      nowpop = gtree[down].mig[mj].mp;
      mj++;
    }
  gtree[sis].mig[mi].mt = -1.0;
  gtree[down].mig[mj].mt = -1.0;

  if (L[li].model == INFINITESITES)
    {
      /* Allocate proper things at the internal node sequence. */
      for (si = 0; si < L[li].numsites; si++)
        {
          if (G->mut[si] == sis)
            {
              if (gtree[sis].seq[si] == seqz[si])
                {
                  G->mut[si] = down;
                }
            }
        }
    }
  else if (L[li].model == HKY || L[li].model == STEPWISE)
    {
      /* No code. */
    }
  else
    {
      assert (0);
    }
  return;
}


int
IMA_genealogy_lineages (int ci, int li, int pi, double td, int *seq) 
{
  int n;
  struct genealogy *G;
  struct edge *gtree;
  int ei;
  int mi;
  int pc;
  int ngnodes;
  double upt;
  double dnt;
  double bt;
  double et;
  int ki;                /* period */
  char is_exist;
  int down;

  G = &C[ci]->G[li];
  gtree = G->gtree;
  ki = findperiodoftime (ci, td);
  ngnodes = L[li].numlines;
  n = 0;
  for (ei = 0; ei < ngnodes; ei++)
    {
      if (gtree[ei].exist == 'F')
        {
          continue;
        }
      is_exist = 'F';
      upt = IMA_genealogy_time (ci, li, ei);
      dnt = gtree[ei].time;
      pc = gtree[ei].pop;
      if (upt < td && td < dnt)
        {
          bt = upt;
          if (gtree[ei].cmm > 0)
            {
              mi = 0;
              while (is_exist == 'F' && gtree[ei].mig[mi].mt > -0.5)
                {
                  et = gtree[ei].mig[mi].mt;
                  if (bt < td && td < et)
                    {
                      pc = saC.popndown[pc][ki];
                      if (pc == pi)
                        {
                          is_exist = 'T';
                        }
                    }
                  pc = gtree[ei].mig[mi].mp;
                  mi++;
                  bt = et;
                } 
            }
          if (is_exist == 'F')
            {
              et = dnt;
              if (bt < td && td < et)
                {
                  pc = saC.popndown[pc][ki];
                  if (pc == pi)
                    {
                      is_exist = 'T';
                    }
                }
            }
        }
      if (is_exist == 'T')
        {
          if (seq == NULL)
            {
              slineages[n] = ei;
              n++;
            }
          else
            {
              /* Infinte-sites data check! */
              if (L[li].model == INFINITESITES)
                {
                  down = gtree[ei].down;
                  if (memcmp (seq, gtree[ei].seq, L[li].numsites * sizeof (int)) == 0)
                    {
                      slineages[n] = ei;
                      n++;
                    }
                  else
                    {
                      if (down == -1) /* We should find its bottom sequence. */
                        {
                          /* FIXME: */
                          assert (0);
                          if (memcmp (seq, gtree[down].seq, L[li].numsites * sizeof (int)) == 0)
                            {
                              slineages[n] = ei;
                              n++;
                            }
                        }
                      else
                        {
                          if (memcmp (seq, gtree[down].seq, L[li].numsites * sizeof (int)) == 0)
                            {
                              slineages[n] = ei;
                              n++;
                            }
                        }
                    }
                }
            }
        }
    }

  return n;
}

int
IMA_genealogy_migrates (int ci, int ki, int pi)
{
  int n;
  int pop;
  int npopi;
  int i;

  npopi = npops - ki;
  assert (npopi > 1);
  n = 0;
  for (i = 0; i < npopi; i++)
    {
      pop = C[ci]->plist[ki][i];
      if (pop == pi)
        {
          continue;
        }
      else
        {
          smigrates[n] = pop;
          n++;
        }
    }
  return n;
}

int
IMA_genealogy_samplepop (int ci, int li, int z, double td)
{
  struct genealogy *G;
  struct edge *gtree;
  int ei;
  int mi;
  int pc;
  int ngnodes;
  double upt;
  double dnt;
  double bt;
  double et;
  int ki;                /* period */
  int pi;
  int n;
  char is_exist;

  for (pi = 0; pi < numtreepops; pi++)
    {
      snsample_pop[pi] = 0;
    }

  G = &C[ci]->G[li];
  gtree = G->gtree;
  ki = findperiodoftime (ci, td);
  ngnodes = L[li].numlines;
  n = 0;
  for (ei = 0; ei < ngnodes; ei++)
    {
      if (ei == z)
        {
          pc = gtree[ei].pop;
          snsample_pop[pc]++;
          continue;
        }
      is_exist = 'F';
      upt = IMA_genealogy_time (ci, li, ei);
      dnt = gtree[ei].time;
      pc = gtree[ei].pop;
      if (upt < td && td < dnt)
        {
          bt = upt;
          if (gtree[ei].cmm > 0)
            {
              mi = 0;
              while (is_exist == 'F' && gtree[ei].mig[mi].mt > -0.5)
                {
                  et = gtree[ei].mig[mi].mt;
                  if (bt < td && td < et)
                    {
                      pc = saC.popndown[pc][ki];
                      snsample_pop[pc]++;
                      is_exist = 'T';
                    }
                  pc = gtree[ei].mig[mi].mp;
                  mi++;
                } 
            }
          else
            {
              et = dnt;
              if (bt < td && td < et)
                {
                  pc = saC.popndown[pc][ki];
                  snsample_pop[pc]++;
                  is_exist = 'T';
                }
            }
        }
    }

  return n;
}
int
IMA_intervals_collect (int ci, int li, double upt, double dnt)
{
  struct genealogy *G;
  struct edge *gtree;
  im_event *ek;
  int enmax;
  int ki;
  int mi;
  int ei;
  int a;
  int ngnodes;
  int ngenes;
  double t;

  G = &C[ci]->G[li];
  gtree = G->gtree;
  ngnodes = L[li].numlines;
  ngenes = L[li].numgenes;

  a = 0;
  ek = sbfold.e1;
  enmax = sbfold.nmax1;
  for (ei = ngenes; ei < ngnodes; ei++)
    {
      if (gtree[ei].exist == 'F')
        {
          continue;
        }
      t = IMA_genealogy_time (ci, li, ei);
      assert (upt < t);
      ek->type = 'c'; 
      ek->t = t;
      ek->p = -1;
      ek->pi = -1;
      ek->pj = gtree[ei].pop;
      ek->ei = ei;
      ek++;
      a++;
    }
  for (ki = 0; ki < numsplittimes; ki++)
    {
      t = C[ci]->tvals[ki];
      if (upt < t && t < dnt)
        {
          ek->type = 's'; 
          ek->t = t;
          ek->p = ki;    /* Split time happens right before the period line. */
          ek->pi = ki;   /* from 0 (pi) to 1 (pj) for the first split time.  */
          ek->pj = ki + 1;
          ek->ei = -1;
          ek++;
          a++;
        }
    }
  for (ei = 0; ei < ngnodes; ei++)
    {
      if (gtree[ei].exist == 'F')
        {
          continue;
        }
      mi = 0;
      assert (gtree[ei].cmm > 0);
      t = gtree[ei].mig[mi].mt;
      while (t > -0.5)
        {
          assert (upt < t);
          if (!(a < enmax - 1))
            {
              enmax += 5;
              sbfold.e1 = realloc (sbfold.e1, enmax * sizeof (im_event));
              sbfold.nmax1 = enmax;

              /* IMPORTANT! The pointer must repoint to a right place. */
              ek = &sbfold.e1[a]; 
            }

          ek->type = 'm';
          ek->t = t;
          ek->p = findperiodoftime (ci, t);
          if (mi == 0)
            {
              ek->pi = gtree[ei].pop;
            }
          else
            {
              ek->pi = gtree[ei].mig[mi - 1].mp;
            }
          ek->pi = saC.popndown[ek->pi][ek->p];
          ek->pj = gtree[ei].mig[mi].mp;
          ek->ei = ei;
          ek++;
          a++;

          mi++;
          t = gtree[ei].mig[mi].mt;
        }
      if (ei == G->root)
        {
          if (sbfold.is_finiterootbranch == 'T')
            {
              t = gtree[ei].time;
              assert (t < TIMEMAX - 1.0);
              if (!(a < enmax - 1))
                {
                  enmax += 5;
                  sbfold.e1 = realloc (sbfold.e1, enmax * sizeof (im_event));
                  sbfold.nmax1 = enmax;

                  /* IMPORTANT! The pointer must repoint to a right place. */
                  ek = &sbfold.e1[a];
                }

              ek->type = 'e';
              ek->t = t;
              ek->p = findperiodoftime (ci, t);
              if (mi == 0)
                {
                  ek->pi = gtree[ei].pop;
                }
              else
                {
                  ek->pi = gtree[ei].mig[mi - 1].mp;
                }
              ek->pi = saC.popndown[ek->pi][ek->p];
              ek->pj = ek->pi;
              ek->ei = ei;
              ek++;
              a++;
            }
        }
    }

  assert (a > 1);
  sbfold.n1 = a;
  qsort (sbfold.e1, (size_t) a, sizeof (im_event), im_event_cmp);
  assert (upt < sbfold.e1[0].t);
  
  return a;
}
int
IMA_genealogy_derivetheta (char w, int ci, int li, int z1)
{
  im_event *ek;
  double prevt;
  double tau;
  double h2term;
  int pi;
  int pj;
  int ki;
  int ai;
  int ei;
  int pii;
  double fm1;
  double fm1z1;
  int ngenes;
  int npminus1; 
  int ninterval;
  struct edge *gtree;
  int ni;
  int pk;
  int pl;
  int ti;
  im_ginfo *gk;
  im_ginfo *gkz1;
  im_bfupdate *pbf;
  int n1;
  int n2;
  int edge;

  assert (assignmentoptions[POPULATIONASSIGNMENT] == 1);
  gtree = C[ci]->G[li].gtree;

  if (w == 'o')
    {
      pbf = &sbfold;
      gk = &saGiold[li];
      gkz1 = &saGioldz1[li];
    }
  else if (w == 'n')
    {
      pbf = &sbfnew;
      gk = &saGinew[li];
      gkz1 = &saGinewz1[li];
    }
  else
    {
      assert (0);
    }

  ek = pbf->e1;
  ninterval = pbf->n1;
  ngenes = L[li].numgenes;
  IMA_sbf_lineages (pbf, ci, li, z1);
  IMA_memory_resetginfo (gk);
  IMA_memory_resetginfo (gkz1);

  h2term = 1 / (2 * L[li].hval);
  prevt = 0.0;
  ki = 0;
  npminus1 = npops - ki - 1;
  for (ai = 0; ai < ninterval; ai++)
    {
      tau = ek->t - prevt;
      assert (tau > 0.0);
      /* For all pop'n during period ki */
      for (pii = 0; pii < npops - ki; pii++) 
        {
          pi = C[ci]->plist[ki][pii];

          /* Jody's code does not have plus one or +1. When we have only single
           * lineage, we should have positve coalescent rate. 1 * (1 - 1) is
           * zero, and we would have 0 rate of coalescent for a single lineage.
           * That may be why we add one to the number of lineages. The idea is
           * that we should consider the active lineage. */
          gk->fc[pi] += snnminus1[pbf->n[pi] + 1] * tau * h2term; 
          gkz1->fc[pi] += snnminus1[pbf->nz1[pi] + 1] * tau * h2term; 
          assert (gk->fc[pi] < DBL_MAX);
          assert (gkz1->fc[pi] < DBL_MAX);
          if (modeloptions[NOMIGRATION] == 0 && ki < lastperiodnumber)
            {
              fm1 = (pbf->n[pi] + 1) * tau;
              fm1z1 = (pbf->nz1[pi] + 1) * tau;
              for (ni = 0; ni < npminus1; ni++)
                {
                  pj = saC.popnmig[ki][pi][ni];
                  gk->fm[pi][pj] += fm1;
                  gkz1->fm[pi][pj] += fm1z1;
                }
            }
        }

      switch (ek->type)
        {
        case 'c':
          pi = ek->pj;
          gk->cc[pi]++;
          gkz1->cc[pi]++;
          edge = ek->ei;
          if (L[li].model == INFINITESITES)
            {
              /* NOTE: We do not update sbfold.l!!! */
              pbf->n[pi]--;
              IMA_sbf_c (ci, li, edge, pi, pbf->seqz1,
                         pbf->nz1, pbf->lz1, pbf->a);
            }
          else if (L[li].model == HKY || L[li].model == STEPWISE)
            {
              /* NOTE: We do not update sbfold.nz1 .lz1 !!! */
              IMA_sbf_c (ci, li, edge, pi, NULL, 
                         pbf->n, pbf->l, pbf->a);
            }
          break;
        case 'm':
          edge = ek->ei;
          pi = ek->pi;
          pj = ek->pj;
          gk->cm[pi][pj]++;
          gkz1->cm[pi][pj]++;
          if (L[li].model == INFINITESITES)
            {
              pbf->n[pi]--;
              pbf->n[pj]++;
              IMA_sbf_m (ci, li, edge, pi, pj, 
                         pbf->nz1, pbf->lz1, pbf->a);
            }
          else if (L[li].model == HKY || L[li].model == STEPWISE)
            {
              IMA_sbf_m (ci, li, edge, pi, pj, 
                         pbf->n, pbf->l, pbf->a);
            }
          break;
        case 's':
          assert (ek->p == ki);
          ki++;
          pj = C[ci]->droppops[ki][IM_EDGE_LEFT];
          pk = C[ci]->droppops[ki][IM_EDGE_RIGHT];
          pl = C[ci]->addpop[ki];
          pbf->n[pl] = pbf->n[pj] + pbf->n[pk];
          n1 = pbf->n[pj];
          n2 = pbf->n[pk];
          if (n1 > 0)
            {
              memcpy (pbf->l[pl], pbf->l[pj], n1 * sizeof (int));
            }
          if (n2 > 0)
            {
              memcpy (&pbf->l[pl][n1], pbf->l[pk], n2 * sizeof (int));
            }
          pbf->n[pj] = 0;
          pbf->n[pk] = 0;
          npminus1--;
          assert (!(npminus1 < 0));
          if (!(z1 < 0))
            {
              /* *pc1 = saC.popndown[*pc1][ki]; */
            }
          assert (L[li].model == INFINITESITES 
                  || L[li].model == HKY
                  || L[li].model == STEPWISE);
          if (!(z1 < 0))
            {
              /* Change the number and sets as well. */
              n1 = pbf->nz1[pj];
              n2 = pbf->nz1[pk];
              pbf->nz1[pl] = n1 + n2;
              if (n1 > 0)
                {
                  memcpy (pbf->lz1[pl], pbf->lz1[pj], n1 * sizeof (int));
                }
              if (n2 > 0)
                {
                  memcpy (&pbf->lz1[pl][n1], pbf->lz1[pk], n2 * sizeof (int));
                }
              pbf->nz1[pj] = 0;
              pbf->nz1[pk] = 0;
            }
          break;
        case 'e':
          assert (ai == ninterval - 1);
          break;
        }
      prevt = ek->t;
      ek++;
    }

  IMA_ginfo_computetheta (gk);
  IMA_ginfo_computetheta (gkz1);
  
  return 0;
}

int
IMA_ginfo_computetheta (im_ginfo *gk)
{  
  int pi;
  int pj;
  int ti;
  int ni;
  int npminus1;
  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      if (gk->cc[pi] > 0)         /* Case (A) */
        {
          gk->thetas[pi] = 2 * gk->fc[pi] / (gk->cc[pi] + 1);
        }
      else
        {
          if (gk->fc[pi] > 0.0)   /* Case (B) */
            {
              assert (gk->cc[pi] == 0);
              gk->thetas[pi] = 2 * gk->fc[pi] / (1);
              /* gk->thetas[pi] = itheta[pi].pr.max; FIXME: is this a rightful max? */
            }
          else                    /* Case (C) */
            {
              gk->thetas[pi] = 0.0;
            }
        }
      for (pj = 0; pj < numpopsizeparams; pj++)
        {
          if (pi == pj)
            {
              continue;
            }
          if (gk->cm[pi][pj] > 0) /* Case (D) */
            {
              assert (gk->fm[pi][pj] > 0.0);
              gk->mms[pi][pj] = (gk->cm[pi][pj] + 1) / gk->fm[pi][pj];
            }
          else                   /* Case (E) and (F) */
            {
              assert (gk->cm[pi][pj] == 0);
              if (gk->fm[pi][pj] > 0.0)
                {
                  gk->mms[pi][pj] = (1) / gk->fm[pi][pj];
                }
              else
                {
                  gk->mms[pi][pj] = BFMPRIORMIN;
                }
            }
        }
    }

  /* We need number of periods where migrations are allowed. */
  /* npops - 1 for tree model and 1 for tree model.          */
  /* lastperiodnumber is for this                            */
  for (ti = 0; ti < lastperiodnumber; ti++)
    {
      for (pi = 0; pi < numpopsizeparams; pi++)
        {
          if (saC.popnmig[ti][pi] != NULL)
            {
              gk->ms[ti][pi] = 0.0;
              npminus1 = npops - ti - 1;
              for (ni = 0; ni < npminus1; ni++)
                {
                  pj = saC.popnmig[ti][pi][ni];
                  gk->ms[ti][pi] += gk->mms[pi][pj];
                }
            }
        }
    }
  return 0;
}
void 
IMA_memory_initsaGi ()
{
  IMA_memory_initginfo (&saGiold);
  IMA_memory_initginfo (&saGioldz1);
  IMA_memory_initginfo (&saGioldz2);
  IMA_memory_initginfo (&saGinew);
  IMA_memory_initginfo (&saGinewz1);
  IMA_memory_initginfo (&saGinewz2);
  return;
}

void 
IMA_memory_freesaGi ()
{
  IMA_memory_freeginfo (&saGiold);
  IMA_memory_freeginfo (&saGioldz1);
  IMA_memory_freeginfo (&saGioldz2);
  IMA_memory_freeginfo (&saGinew);
  IMA_memory_freeginfo (&saGinewz1);
  IMA_memory_freeginfo (&saGinewz2);
  return;
}

void 
IMA_memory_initginfo (im_ginfo **g)
{
  int li;
  int pi;

  (*g) = malloc (nloci * sizeof (im_ginfo));
  for (li = 0; li < nloci; li++)
    {
      (*g)[li].cc = malloc (numpopsizeparams * sizeof (int));
      (*g)[li].fc = malloc (numpopsizeparams * sizeof (double));
      (*g)[li].cm = malloc (numpopsizeparams * sizeof (int *));
      (*g)[li].fm = malloc (numpopsizeparams * sizeof (double *));
      (*g)[li].mms = malloc (numpopsizeparams * sizeof (double *));
      for (pi = 0; pi < numpopsizeparams; pi++)
        {
          (*g)[li].cm[pi] = malloc (numpopsizeparams * sizeof (int));
          (*g)[li].fm[pi] = malloc (numpopsizeparams * sizeof (double));
          (*g)[li].mms[pi] = malloc (numpopsizeparams * sizeof (double));
        }
      (*g)[li].ms = malloc (lastperiodnumber * sizeof (double *));
      for (pi = 0; pi < lastperiodnumber; pi++)
        {
          (*g)[li].ms[pi] = malloc (numpopsizeparams * sizeof (double));
        }
      (*g)[li].thetas = malloc (numpopsizeparams * sizeof (double));
    }
  return;
}

void
IMA_memory_freeginfo (im_ginfo **g)
{
  int li;
  int pi;

  for (li = 0; li < nloci; li++)
    {
      free ((*g)[li].cc);
      (*g)[li].cc = NULL;
      free ((*g)[li].fc);
      (*g)[li].fc = NULL;
      for (pi = 0; pi < numpopsizeparams ; pi++)
        {
          free ((*g)[li].cm[pi]);
          (*g)[li].cm[pi] = NULL;
          free ((*g)[li].fm[pi]);
          (*g)[li].fm[pi] = NULL;
          free ((*g)[li].mms[pi]);
          (*g)[li].mms[pi] = NULL;
        }
      free ((*g)[li].cm);
      (*g)[li].cm = NULL;
      free ((*g)[li].fm);
      (*g)[li].fm = NULL;
      free ((*g)[li].thetas);
      (*g)[li].thetas = NULL;
      free ((*g)[li].mms);
      (*g)[li].mms = NULL;
      for (pi = 0; pi < lastperiodnumber; pi++)
        {
          free ((*g)[li].ms[pi]);
          (*g)[li].ms[pi] = NULL;
        }
      free ((*g)[li].ms);
      (*g)[li].ms = NULL;
    }
  free (*g);
  (*g) = NULL;
  return;
}

void 
IMA_memory_resetginfo (im_ginfo *g)
{
  int ki;
  int pi;
  int pj;

  g->cc = memset (g->cc, 0, numpopsizeparams * sizeof (int));
  g->fc = memset (g->fc, 0, numpopsizeparams * sizeof (double));
  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      g->cm[pi] = memset (g->cm[pi], 0, numpopsizeparams * sizeof (int));
      g->fm[pi] = memset (g->fm[pi], 0, numpopsizeparams * sizeof (double));
      g->mms[pi] = memset (g->mms[pi], 0, numpopsizeparams * sizeof (double));
    }
  for (ki = 0; ki < lastperiodnumber; ki++)
    {
      g->ms[ki] = memset (g->ms[ki], 0, numpopsizeparams * sizeof (double));
    }
  g->thetas = memset (g->thetas, 0, numpopsizeparams * sizeof (double));

  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      assert (g->cc[pi] == 0);
      assert (g->fc[pi] == 0.0);
      for (pj = 0; pj < numpopsizeparams; pj++)
        {
          assert (g->cm[pi][pj] == 0);
          assert (g->fm[pi][pj] == 0.0);
          assert (g->mms[pi][pj] == 0.0);
        }
      assert (g->thetas[pi] == 0.0);
    }
  for (ki = 0; ki < lastperiodnumber; ki++)
    {
      for (pj = 0; pj < numpopsizeparams; pj++)
        {
          assert (g->ms[ki][pj] == 0.0);
        }
    }

  return;
}

int
ran_discrete (int n, ...)
{
  register int i;
  va_list ap;
  double s;
  double a;
  double l[10];
  double u;
  int v;

  assert (n > 1);
  assert (n <= 10);
  va_start (ap, n);
  a = va_arg (ap, double);
  s = a;
  l[0] = a;
  for (i = 1; i < n; i++)
    {
      a = va_arg (ap, double);
      l[i] = a;
      s += a;
    }
  va_end (ap);

  u = s * uniform ();

  s = 0.0;
  for (i = 0; i < n; i++)
    {
      s += l[i];
      if (u < s)
        {
          v = i;
          break;
        }
    }
  return v;
}

int
ran_discretea (int n, double pr[])
{
  register int i;
  int v;
  double s;
  double u;

  s = 0;
  for (i = 0; i < n; i++)
    {
      s += pr[i]; 
    }
  u = s * uniform ();
  s = 0.0;
  for (i = 0; i < n; i++)
    {
      s += pr[i];
      if (u < s)
        {
          v = i;
          break;
        }
    }
  return v;

}

int
ran_discreteloga (int n, double pr[])
{
  register int i;
  int v;
  double s;
  double t;
  double u;

  t = pr[0];
  for (i = 1; i < n; i++)
    {
      LogSum2 (t, t, pr[i]);
    }
  u = log (uniform ()) + t;

  s = pr[0];
  if (u < s)
    {
      v = 0;
    }
  else
    {
      for (i = 1; i < n; i++)
        {
          LogSum2 (s, s, pr[i]);
          if (u < s)
            {
              v = i;
              break;
            }
        }
    }
  return v;
}

int
normalizeLoga (int n, double pr[])
{
  register int i;
  double t;

  t = pr[0];
  for (i = 1; i < n; i++)
    {
      LogSum2 (t, t, pr[i]);
    }
  
  for (i = 0; i < n; i++)
    {
      pr[i] -= t;
    }

  return 0;
}

int
IMA_create_intervals ()
{
  int n;
  int li;
  int ai;
  int pi;
  int largestnA;

  n = 0;
  largestnA = -1;
  for (li = 0; li < nloci; li++)
    {
      if (n < L[li].numlines)
        {
          n = L[li].numlines;
        }
      if (largestnA < L[li].maxA[0] - L[li].minA[0] + 1)
        {
          largestnA = L[li].maxA[0] - L[li].minA[0] + 1;
        }
    }
  sbfold.nmax1 = 2 * n;
  sbfold.nmax2 = 2 * n;
  sbfold.nmax3 = 2 * n;
  sbfnew.nmax1 = 2 * n;
  sbfnew.nmax2 = 2 * n;
  sbfnew.nmax3 = 2 * n;

  sbfold.e1 = (im_event *) malloc (sbfold.nmax1 * sizeof (im_event));
  sbfold.e2 = (im_event *) malloc (sbfold.nmax2 * sizeof (im_event));
  sbfold.e3 = (im_event *) malloc (sbfold.nmax3 * sizeof (im_event));
  sbfold.nz1 = (int *) malloc (numpopsizeparams * sizeof (int));
  /* sbfold.nz2 = (int *) malloc (numpopsizeparams * sizeof (int)); */
  sbfold.n = (int *) malloc (numpopsizeparams * sizeof (int));
  sbfold.lz1 = (int **) malloc (numpopsizeparams * sizeof (int *));
  /* sbfold.lz2 = (int **) malloc (numpopsizeparams * sizeof (int *)); */
  sbfold.l = (int **) malloc (numpopsizeparams * sizeof (int *));
  sbfold.a = (int *) malloc (gi_largestngenes * sizeof (int));
  sbfold.likes1 = (double *) malloc (gi_largestngenes * sizeof (double));
  sbfold.likes11 = (double **) malloc (gi_largestngenes * sizeof (double *));
  for (ai = 0; ai < gi_largestngenes; ai++)
    {
      sbfold.likes11[ai] = (double *) malloc (largestnA * sizeof (double));
    }

  sbfold.seqz1 = (int *) malloc (gi_largestnumsites * sizeof (int));
  /* sbfold.seqz2 = (int *) malloc (gi_largestnumsites * sizeof (int)); */
  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      assert (gi_largestngenes > 0);
      sbfold.lz1[pi] = (int *) malloc (gi_largestngenes * sizeof (int));
      /* sbfold.lz2[pi] = (int *) malloc (gi_largestngenes * sizeof (int)); */
      sbfold.l[pi] = (int *) malloc (gi_largestngenes * sizeof (int));
    }
  sbfnew.e1 = (im_event *) malloc (sbfnew.nmax1 * sizeof (im_event));
  sbfnew.e2 = (im_event *) malloc (sbfnew.nmax2 * sizeof (im_event));
  sbfnew.e3 = (im_event *) malloc (sbfnew.nmax3 * sizeof (im_event));
  sbfnew.nz1 = (int *) malloc (numpopsizeparams * sizeof (int));
  /* sbfnew.nz2 = (int *) malloc (numpopsizeparams * sizeof (int)); */
  sbfnew.n = (int *) malloc (numpopsizeparams * sizeof (int));
  sbfnew.lz1 = (int **) malloc (numpopsizeparams * sizeof (int *));
  /* sbfnew.lz2 = (int **) malloc (numpopsizeparams * sizeof (int *)); */
  sbfnew.l = (int **) malloc (numpopsizeparams * sizeof (int *));
  sbfnew.a = (int *) malloc (gi_largestngenes * sizeof (int));
  sbfnew.likes1 = (double *) malloc (gi_largestngenes * sizeof (double));
  sbfnew.likes11 = (double **) malloc (gi_largestngenes * sizeof (double *));
  for (ai = 0; ai < gi_largestngenes; ai++)
    {
      sbfnew.likes11[ai] = (double *) malloc (largestnA * sizeof (double));
    }

  sbfnew.seqz1 = (int *) malloc (gi_largestnumsites * sizeof (int));
  /* sbfnew.seqz2 = (int *) malloc (gi_largestnumsites * sizeof (int)); */
  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      assert (gi_largestngenes > 0);
      sbfnew.lz1[pi] = (int *) malloc (gi_largestngenes * sizeof (int));
      /* sbfnew.lz2[pi] = (int *) malloc (gi_largestngenes * sizeof (int)); */
      sbfnew.l[pi] = (int *) malloc (gi_largestngenes * sizeof (int));
    }

  slineages = (int *) malloc (n * sizeof (int));
  smigrates = (int *) malloc (npops * sizeof (int));

  return 0;
}

int
IMA_delete_intervals ()
{
  int pi, ai;

  free (sbfold.e1);
  sbfold.e1 = NULL;
  sbfold.n1 = 0;
  sbfold.nmax1 = 0;
  free (sbfold.e2);
  sbfold.e2 = NULL;
  sbfold.n2 = 0;
  sbfold.nmax2 = 0;
  free (sbfold.e3);
  sbfold.e3 = NULL;
  sbfold.n3 = 0;
  sbfold.nmax3 = 0;
  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      free (sbfold.lz1[pi]);
      sbfold.lz1[pi] = NULL;
      /* free (sbfold.lz2[pi]); */
      /* sbfold.lz2[pi] = NULL; */
      free (sbfold.l[pi]);
      sbfold.l[pi] = NULL;
    }
  free (sbfold.lz1);
  sbfold.lz1 = NULL;
  /* free (sbfold.lz2); */
  /* sbfold.lz2 = NULL; */
  free (sbfold.l);
  sbfold.l = NULL;
  free (sbfold.nz1);
  sbfold.nz1 = NULL;
  /* free (sbfold.nz2); */
  /* sbfold.nz2 = NULL; */
  free (sbfold.n);
  sbfold.n = NULL;
  free (sbfold.a);
  sbfold.a = NULL;
  free (sbfold.likes1);
  sbfold.likes1 = NULL;
  for (ai = 0; ai < gi_largestngenes; ai++)
    {
      free (sbfold.likes11[ai]);
      sbfold.likes11[ai] = NULL;
    }
  free (sbfold.likes11);
  sbfold.likes11 = NULL;

  XFREE (sbfold.seqz1);
  /* XFREE (sbfold.seqz2); */

  free (sbfnew.e1);
  sbfnew.e1 = NULL;
  sbfnew.n1 = 0;
  sbfnew.nmax1 = 0;
  free (sbfnew.e2);
  sbfnew.e2 = NULL;
  sbfnew.n2 = 0;
  sbfnew.nmax2 = 0;
  free (sbfnew.e3);
  sbfnew.e3 = NULL;
  sbfnew.n3 = 0;
  sbfnew.nmax3 = 0;
  for (pi = 0; pi < numpopsizeparams; pi++)
    {
      free (sbfnew.lz1[pi]);
      sbfnew.lz1[pi] = NULL;
      /* free (sbfnew.lz2[pi]); */
      /* sbfnew.lz2[pi] = NULL; */
      free (sbfnew.l[pi]);
      sbfnew.l[pi] = NULL;
    }
  free (sbfnew.lz1);
  sbfnew.lz1 = NULL;
  /* free (sbfnew.lz2); */
  /* sbfnew.lz2 = NULL; */
  free (sbfnew.l);
  sbfnew.l = NULL;
  free (sbfnew.nz1);
  sbfnew.nz1 = NULL;
  /* free (sbfnew.nz2); */
  /* sbfnew.nz2 = NULL; */
  free (sbfnew.n);
  sbfnew.n = NULL;
  free (sbfnew.a);
  sbfnew.a = NULL;
  free (sbfnew.likes1);
  sbfnew.likes1 = NULL;
  for (ai = 0; ai < gi_largestngenes; ai++)
    {
      free (sbfnew.likes11[ai]);
      sbfnew.likes11[ai] = NULL;
    }
  free (sbfnew.likes11);
  sbfnew.likes11 = NULL;

  XFREE (sbfnew.seqz1);
  /* XFREE (sbfnew.seqz2); */

  free (slineages);
  slineages = NULL;

  free (smigrates);
  smigrates = NULL;
  return 0;
}

int 
im_event_cmp (const void *a, const void *b)
{
  const im_event *A = (const im_event *) a;
  const im_event *B = (const im_event *) b;
  if (A->t > B->t)
    return 1;
  else if (A->t < B->t)
    return -1;
  else
    return 0;
}

double
likelihoodJC (int ci, int li, double u)
{
  double l_k[4];
  double v;
  int si;
  im_edge *gtree; 

  gtree = C[ci]->G[li].gtree;
  v = 0.0;
  for (si = 0; si < L[li].numsites; si++)
    {
      computeLk (&l_k, u, li, gtree, C[ci]->G[li].root, si);
      v += logsum (4, l_k[0], l_k[1], l_k[2], l_k[3]);
      v -= M_LN2 * 2;
    }

  return v;
}

/* Function: computeLk
 * l_k: return value for the conditional likelihood
 * u  : branch length scaler
 * li : locus index
 * t  : tree
 * k  : node index
 * si : site index */
int
computeLk (double (*l_k)[], double u, int li, im_edge *t, int k, int si)
{
  int s;
  int x;
  int y;
  int left;
  int right;
  int base;
  int basel;
  int baser;
  double l_kl[4] = { 0.0, 0.0, 0.0, 0.0};
  double l_kr[4] = { 0.0, 0.0, 0.0, 0.0};
  double l_kleft;
  double l_kright;
  double A;
  double tl;
  double tr;

  left = t[k].up[0];
  right = t[k].up[1];
  base = -1;
  if (left < 0 && right < 0)
    {
      assert (!(k < 0));
      assert (k < L[li].numgenes);
      base = L[li].seq[k][si];
    }
  else
    {
      assert (!(left < 0 || right < 0));
      basel = computeLk (&l_kl, u, li, t, left, si);
      baser = computeLk (&l_kr, u, li, t, right, si);
      tl = u * IMA_edge_length (t, left);
      tr = u * IMA_edge_length (t, right);

      for (s = 0; s < 4; s++)
        {
          l_kleft = 0.0;
          if (basel == -1)
            {
              for (x = 0; x < 4; x++)
                {
                  A = PijJC (x, s, tl) + l_kl[x];
                  LogSum2 (l_kleft, l_kleft, A);
                }
            }
          else
            {
              assert (basel == 0 || basel == 1 || basel == 2 || basel == 3);
              l_kleft = PijJC (basel, s, tl);
            }

          l_kright = 0.0;
          if (baser == -1)
            {
              for (y = 0; y < 4; y++)
                {
                  A = PijJC (y, s, tr) + l_kr[y];
                  LogSum2 (l_kright, l_kright, A);
                }
            }
          else
            {
              assert (baser == 0 || baser == 1 || baser == 2 || baser == 3);
              l_kright = PijJC (baser, s, tr);
            }
          (*l_k)[s] = l_kleft + l_kright;
        }
    }

  return base; 
}

double
PijJC (int i, int j, double t)
{
  double v;
  double A;

  if (i == j)
    {
      A = -4.0 * t / 3.0;
      LogDiff (v, 0.0, A);
      v -= M_LN2 * 2;
    }
  else
    {
      A = -4.0 * t / 3.0 + log (3.0);
      LogSum2 (v, 0.0, A);
      v -= M_LN2 * 2;
    }
  return v;
}

double
IMA_edge_length (im_edge *t, int ei)
{
  double upt;
  int up0;

  up0 = t[ei].up[0];
  if (up0 < 0)
    {
      upt = 0.0;
    }
  else
    {
      upt = t[up0].time;
    }
  return t[ei].time - upt;
}

double
likelihoodSWP (int ci, int li, int ai, double u)
{
  double *l_k;
  double v;
  int si;
  im_edge *gtree; 
  int nallele;

  nallele = L[li].maxA[ai] - L[li].minA[ai];
  assert (nallele > 1);
  assert (nallele < 10);
  l_k = (double *) malloc (nallele * sizeof (double));
  gtree = C[ci]->G[li].gtree;
  computeLkSWP (l_k, u, li, ai, gtree, C[ci]->G[li].root, nallele);
  v = logsuma (nallele, l_k);

  free (l_k);
  l_k = NULL;

  return v;
}

/* Function: computeLkSWP
 * l_k: return value for the conditional likelihood
 * u  : branch length scaler
 * li : locus index
 * t  : tree
 * k  : node index
 * si : site index */
int
computeLkSWP (double *l_k, double u, int li, int ai, im_edge *t, int k, int nallele)
{
  int s;
  int x;
  int y;
  int left;
  int right;
  int base;
  int basel;
  int baser;
  double *l_kl;
  double *l_kr;
  double l_kleft;
  double l_kright;
  double A;
  double tl;
  double tr;

  l_kl = (double *) malloc (nallele * sizeof (double));
  l_kr = (double *) malloc (nallele * sizeof (double));
  memset (l_kl, 0, nallele * sizeof (double));
  memset (l_kr, 0, nallele * sizeof (double));

  left = t[k].up[0];
  right = t[k].up[1];
  base = -1;
  if (left < 0 && right < 0)
    {
      assert (!(k < 0));
      assert (k < L[li].numgenes);
      base = t[k].A[ai];
    }
  else
    {
      assert (!(left < 0 || right < 0));
      basel = computeLkSWP (l_kl, u, li, ai, t, left, nallele);
      baser = computeLkSWP (l_kr, u, li, ai, t, right, nallele);
      tl = u * IMA_edge_length (t, left);
      tr = u * IMA_edge_length (t, right);

      for (s = 0; s < nallele; s++)
        {
          l_kleft = 0.0;
          if (basel == -1)
            {
              for (x = 0; x < nallele; x++)
                {
                  A = PijSWP (x, s, tl) + l_kl[x];
                  LogSum2 (l_kleft, l_kleft, A);
                }
            }
          else
            {
              assert (basel == 0 || basel == 1 || basel == 2 || basel == 3);
              l_kleft = PijSWP (basel, s, tl);
            }

          l_kright = 0.0;
          if (baser == -1)
            {
              for (y = 0; y < nallele; y++)
                {
                  A = PijSWP (y, s, tr) + l_kr[y];
                  LogSum2 (l_kright, l_kright, A);
                }
            }
          else
            {
              assert (baser == 0 || baser == 1 || baser == 2 || baser == 3);
              l_kright = PijSWP (baser, s, tr);
            }
          l_k[s] = l_kleft + l_kright;
        }
    }
  free (l_kl);
  l_kl = NULL;
  free (l_kr);
  l_kr = NULL;

  return base; 
}

double
PijSWP (int i, int j, double t)
{
  double bessiv;
  int d;

  d = abs (i - j);
  bessiv = bessi (d, t);
  return bessiv;
}


void 
IMA_structurama_uncertainty (char *fname)
{
  int d;
  FILE *fp;
  int lenpname;
  char *pname;
  int **m;
  int nrow;
  int ncol;
  int i;
  int *B;
  int total_distance;

  double *au;

  lenpname = 20 + strlen (fname);

  pname = malloc (lenpname * sizeof (char));
  sprintf (pname, "%s.in.p", fname);
  IMA_io_readp (pname, &m, &nrow, &ncol);

  B = malloc (ncol * sizeof (int));
  au = malloc (ncol * sizeof (double));
  memset (au, 0, ncol * sizeof (double));
  sprintf (pname, "%s.in.sum_assignments", fname);
  IMA_io_readsumassignment (pname, ncol, B);
  d = IMA_gusfield_distanceRemoved (B, m, nrow, ncol, au);

  sprintf (pname, "%s.uct", fname);
  fp = fopen (pname, "w");

  fprintf (fp, "t = {");
  for (i = 0; i < ncol; i++)
    {
      fprintf (fp, "{%d, %lf}", B[i], au[i]);
      if (i < ncol - 1)
        {
          fprintf (fp, ", ");
        }
    }
  fprintf (fp, "}\n");
  fclose (fp);
  fp = NULL;
  free (au);
  au = NULL;
  free (B);
  B = NULL;
  for (i = 0; i < nrow; i++)
    {
      free (m[i]);
      m[i] = NULL;
    }
  free (m);
  m = NULL;
  free (pname);
  pname = NULL;


  return;
}

void 
IMA_structurama_summary (char *fname)
{
  FILE *fp;
  int lenpname;
  char *pname;
  int **m;
  int nrow;
  int ncol;
  int i;
  int *A;
  int *B;
  int total_distance;
  double mean_distance;
  double normalized_distance;
  double variance_distance;
  double squared_distance;
  double meansquared_distance;
  int b_distance; /* distance between true and estimate */

  lenpname = 20 + strlen (fname);

  pname = malloc (lenpname * sizeof (char));
  sprintf (pname, "%s.in.p", fname);

  IMA_io_readp (pname, &m, &nrow, &ncol);
  if (strlen (trueassignment) == ncol)
    {
      A = malloc (ncol * sizeof (int));
      IMA_convertStr2Asn (ncol, A, trueassignment);

      total_distance = IMA_gusfield_distance (A, m, nrow, ncol); 
      mean_distance = ((double) total_distance) / ((double) nrow);
      normalized_distance = mean_distance / ncol;
      variance_distance = IMA_gusfield_varianceDistance (A, m, nrow, ncol, mean_distance); 
      squared_distance = IMA_gusfield_squaredDistance (A, m, nrow, ncol); 
      meansquared_distance = ((double) squared_distance) / ((double) nrow);
      sprintf (pname, "%s.sum", fname);
      fp = fopen (pname, "w");
      fprintf (fp, "%d\n", total_distance);
      fprintf (fp, "%g\n", mean_distance);
      fprintf (fp, "%g\n", normalized_distance);
      fprintf (fp, "%g\n", variance_distance);
      fprintf (fp, "%g\n", squared_distance);
      fprintf (fp, "%g\n", meansquared_distance);

      B = malloc (ncol * sizeof (int));
      sprintf (pname, "%s.in.sum_assignments", fname);
      IMA_io_readsumassignment (pname, ncol, B);
      b_distance = IMA_gdistance (ncol, A, B);
      fprintf (fp, "%d\n", b_distance);
      free (A);
      A = NULL;
      free (B);
      B = NULL;

      fclose (fp);
      fp = NULL;
    }
  free (A);
  A = NULL;
  for (i = 0; i < nrow; i++)
    {
      free (m[i]);
      m[i] = NULL;
    }
  free (m);
  m = NULL;
  free (pname);
  pname = NULL;
  return;
}

void 
IMA_io_readsumassignment (char *fname, int n, int *A)
{
  FILE *fp;
  char c;
  int i, d;
  
  fp = fopen (fname, "r");
  for (i = 0; i < 10; i++)
    {
      skip_a_line (fp);
    }
  for (i = 0; i < n; i++)
    {
      c = 'x';
      while (c != ')')
        {
          c = fgetc (fp);
        }
      d = read_int (fp);
      A[i] = d;
    }
  fclose (fp);
  fp = NULL;
 
  return;
}

int
IMA_convertStr2Asn (int ncol, int *a, char *s)
{
  int len;
  int i;
  char str[2];

  len = strlen (s);
  assert (len == ncol);

  str[1] = '\0';
  for (i = 0; i < len; i++)
    {
      str[0] = s[i];
      a[i] = atoi (str);
    }
  return 1;
}

/*
static int
skip_a_line (FILE *fp)
{
  char c;
  c = 'x';
  while (c != '\n')
    {
      c = fgetc(fp);
    }
  return 0;
}
*/

int
IMA_io_readp (char *fname, int ***m, int *mr, int *mc)
{
  FILE *fp;
  char c;
  int ncol, nrow;
  int i, j;
  double lnL;
  int gi;
  
  fp = fopen (fname, "r");
  
  /* Count the number of columns. */
  c = fgetc(fp);
  ncol = 1;
  while (c != '\n')
    {
      if (c == '\t')
        {
          ncol++;
        }
      c = fgetc(fp);
    }

/* DOS?
  ncol -= 3;
  nrow = -2;
*/
  ncol -= 2;
  nrow = -1;

  while (!feof (fp))
    {
      if (c == '\n')
        {
          nrow++;
        }
      c = fgetc(fp);      
    }
  
  (*m) = malloc (nrow * sizeof (int *));
  for (i = 0; i < nrow; i++)
    {
      (*m)[i] = malloc (ncol * sizeof (int));
    }
  
  /* Reposition pointer to the beginning of the file. */
  fseek (fp, 0L, SEEK_SET);

/* DOS? Remove this line. */
  skip_a_line (fp);

  for (i = 0; i < nrow; i++)
    {
/* DOS? Remove the skip_a_line before for statement.
      skip_a_line (fp);
*/
      gi = read_int (fp);
      lnL = read_double (fp);
      for (j = 0; j < ncol; j++)
        {
          gi = read_int (fp);
          (*m)[i][j] = gi;
        }
    }
  
  fclose (fp);
  fp = NULL;
  
  *mr = nrow;
  *mc = ncol;
  fprintf (stdout, "ncol: %5d\nnrow: %5d\n", ncol, nrow);
  return 0;
}

int
IMA_gusfield_distanceRemoved (int *A, int **m, int nrow, int ncol, double *r)
{
  int d;
  int i, j;
  char *c;
  
  c = malloc (ncol * sizeof (char));
  d = 0;
  for (i = 0; i < nrow; i++)
    {
      d += IMA_gdistanceRemoved (ncol, A, m[i], c);
      for (j = 0; j < ncol; j++)
        {
          if (c[j] == 'T')
            {
              r[j] += 1.0;
            }
        }
    }
  for (j = 0; j < ncol; j++)
    {
      r[j] /= ((double) nrow);
    }
  
  free (c);
  c = NULL;
  return d;
}

int
IMA_gusfield_distance (int *A, int **m, int nrow, int ncol)
{
  int d;
  int i;
  
  d = 0;
  for (i = 0; i < nrow; i++)
    {
      d += IMA_gdistance (ncol, A, m[i]);
    }
  return d;
}

int
IMA_gusfield_squaredDistance (int *A, int **m, int nrow, int ncol)
{
  int d, d1;
  int i;
  
  d = 0;
  for (i = 0; i < nrow; i++)
    {
      d1 = IMA_gdistance (ncol, A, m[i]);
      d1 *= d1;
      d += d1;
    }
  return d;
}


double
IMA_gusfield_varianceDistance (int *A, int **m, int nrow, int ncol, double mean)
{
  double v;
  double xi;
  int i;
  
  v = 0.0;
  for (i = 0; i < nrow; i++)
    {
      xi = (double) IMA_gdistance (ncol, A, m[i]);
      v += ((xi - mean) * (xi - mean));
    }
  v /= ((double) (nrow - 1.0)); 
  return v;
}

int
IMA_gdistanceRemoved (int n, int *Ai, int *Bi, char *Ci)
{
  int rk; /* rank for the best permutation */
  int idummy;
  int i, j, k;
  int d, d1, d2;
  int r, s;
  UByteP *ar;
  UByteP *bs;
  UByteP c;
  int sP;
  int g;
  int *pi;
  int npi; /* number of permutations */
  int *A;
  int *B;
  /* Find r and s. */
  r = IMA_max_intarray (n, Ai);
  s = IMA_max_intarray (n, Bi);
  if (r > s)
    {
      A = Bi;
      B = Ai;
      i = s;
      s = r;
      r = i;
    }
  else
    {
      A = Ai;
      B = Bi;
    }
  assert (!(r > s));
  sP = (n - 1) / 8 + 1;
  BitPowerSetNew (ar, sP, r, i);
  BitPowerSetNew (bs, sP, s, i);
  BitSetNew (c, sP);
  /* pi is initialized to be the first of permutations. */
  pi = malloc ((s + 1) * sizeof (int));
  for (i = 0; i < s + 1; i++)
    {
      pi[i] = i;
    }
  
  /* Set a and b. */
  for (i = 0; i < n; i++)
    {
      g = A[i];
      assert (!(g > r));
      BitTrue (ar[g - 1], i);
      g = B[i];
      assert (!(g > s));
      BitTrue (bs[g - 1], i);
    }
  
  d = n;
  npi = IMA_factorial(s);
  i = 0;
  rk = i;
  while (i++ < npi)
    {
      d2 = 0;
      for (j = 1; j < s + 1; j++)
        {
          k = pi[j] - 1;
          if (j < r + 1)
            {
              BitDifference (c, ar[j - 1], bs[k], sP, idummy);
              BitNTrues (d1, c, sP, idummy);
              d2 += d1;              
            }
        }
      if (d2 < d)
        {
          rk = i;
          d = d2;
        }
      IMA_permLexSuccessor(s, pi);
    }
  
  assert (rk > 0);
  IMA_permLexUnrank(pi, s, rk - 1);
  memset (Ci, 'F', n * sizeof (char));
  for (j = 1; j < s + 1; j++)
    {
      k = pi[j] - 1;
      if (j < r + 1)
        {
          BitDifference (c, ar[j - 1], bs[k], sP, idummy);
          for (i = 0; i < n; i++)
            {
              if (BitIsTrue (c, i))
                {
                  Ci[i] = 'T';
                }
            }          
        }          
    }  
  
  free (pi);
  pi = NULL;
  BitSetDelete (c);
  BitPowerSetDelete (ar, r, i);
  BitPowerSetDelete (bs, s, i);
  return d;
}

/* Return a swap order not the distance between two assignment. */
int
IMA_gdistanceOrder (int n, int *Ai, int *Bi, int m, int *od)
{
  int idummy;
  int i, j, k;
  int d, d1, d2;
  int r, s;
  UByteP *ar;
  UByteP *bs;
  UByteP c;
  int sP;
  int g;
  int *pi;
  int npi; /* number of permutations */
  int *A;
  int *B;
  int v;
  
  /* Find r and s. */
  r = IMA_max_intarray (n, Ai);
  s = IMA_max_intarray (n, Bi);
  if (r > s)
    {
      A = Bi;
      B = Ai;
      i = s;
      s = r;
      r = i;
    }
  else
    {
      A = Ai;
      B = Bi;
    }
  assert (!(r > s));
  /* r is not greater than s. */
  assert (m == s);

  sP = (n - 1) / 8 + 1;
  BitPowerSetNew (ar, sP, r, i);
  BitPowerSetNew (bs, sP, s, i);
  BitSetNew (c, sP);
  /* pi is initialized to be the first of permutations. */
  pi = malloc ((s + 1) * sizeof (int));
  for (i = 0; i < s + 1; i++)
    {
      pi[i] = i;
    }
  
  /* Set a and b. */
  for (i = 0; i < n; i++)
    {
      g = A[i];
      assert (!(g > r));
      BitTrue (ar[g - 1], i);
      g = B[i];
      assert (!(g > s));
      BitTrue (bs[g - 1], i);
    }
  
  d = n;
  npi = IMA_factorial(s);
  assert (npi > 0);
  i = 0;
  v = 0;
  while (i++ < npi)
    {
      d2 = 0;
      for (j = 1; j < s + 1; j++)
        {
          k = pi[j] - 1;
          if (j < r + 1)
            {
              BitDifference (c, ar[j - 1], bs[k], sP, idummy);
              BitNTrues (d1, c, sP, idummy);
              d2 += d1;              
            }
        }

      if (d2 < d)
        {
          d = d2;
          v = i;
        }
      IMA_permLexSuccessor(s, pi);
    }

  IMA_permLexRank (v, pi);
  d2 = 0;
  for (j = 1; j < s + 1; j++)
    {
      k = pi[j] - 1;
      od[j - 1] = k;
      if (j < r + 1)
        {
          BitDifference (c, ar[j - 1], bs[k], sP, idummy);
          BitNTrues (d1, c, sP, idummy);
          d2 += d1;              
        }
    }
  assert (d == d2);
  
  free (pi);
  pi = NULL;
  BitSetDelete (c);
  BitPowerSetDelete (ar, r, i);
  BitPowerSetDelete (bs, s, i);
  return v;
}

int
IMA_gdistance (int n, int *Ai, int *Bi)
{
  int idummy;
  int i, j, k;
  int d, d1, d2;
  int r, s;
  UByteP *ar;
  UByteP *bs;
  UByteP c;
  int sP;
  int g;
  int *pi;
  int npi; /* number of permutations */
  int *A;
  int *B;
  
  /* Find r and s. */
  r = IMA_max_intarray (n, Ai);
  s = IMA_max_intarray (n, Bi);
  if (r > s)
    {
      A = Bi;
      B = Ai;
      i = s;
      s = r;
      r = i;
    }
  else
    {
      A = Ai;
      B = Bi;
    }
  assert (!(r > s));
  sP = (n - 1) / 8 + 1;
  BitPowerSetNew (ar, sP, r, i);
  BitPowerSetNew (bs, sP, s, i);
  BitSetNew (c, sP);
  /* pi is initialized to be the first of permutations. */
  pi = malloc ((s + 1) * sizeof (int));
  for (i = 0; i < s + 1; i++)
    {
      pi[i] = i;
    }
  
  /* Set a and b. */
  for (i = 0; i < n; i++)
    {
      g = A[i];
      assert (!(g > r));
      BitTrue (ar[g - 1], i);
      g = B[i];
      assert (!(g > s));
      BitTrue (bs[g - 1], i);
    }
  
  d = n;
  npi = IMA_factorial(s);
  assert (npi > 0);
  i = 0;
  while (i++ < npi)
    {
      d2 = 0;
      for (j = 1; j < s + 1; j++)
        {
          k = pi[j] - 1;
          if (j < r + 1)
            {
              BitDifference (c, ar[j - 1], bs[k], sP, idummy);
              BitNTrues (d1, c, sP, idummy);
              d2 += d1;              
            }
        }
      if (d2 < d)
        {
          d = d2;
        }
      IMA_permLexSuccessor(s, pi);
    }
  
  free (pi);
  pi = NULL;
  BitSetDelete (c);
  BitPowerSetDelete (ar, r, i);
  BitPowerSetDelete (bs, s, i);
  return d;
}

int
IMA_gdistanceSimple (int n, int *Ai, int *Bi)
{
  int d;
  int i;

  assert (n > 0);
  d = 0;
  for (i = 0; i < n; i++)
    {
      if (Ai[i] != Bi[i])
        {
          d++;
        }
    }
  return d;
}

int
IMA_max_intarray (int n, int *A)
{
  int i;
  int m;
  m = A[0];
  for (i = 1; i < n; i++)
    {
      if (m < A[i])
        {
          m = A[i];
        }
    }
  return m;
}
int
IMA_permFprintf (FILE *fp, int n, int *p)
{
  int i;
  
  for (i = 1; i < n + 1; i++)
    {
      if (i == n)
        {
          fprintf (fp, "%2d\n", p[i]);
        }
      else
        {
          fprintf (fp, "%2d ", p[i]);          
        }
    }
  return 0;
}

/* Algorithm 2.14 of Kreher and Stinson */
int
IMA_permLexSuccessor (int n, int *p)
{
  int i, j, t, h;
  int *rho;
  
  rho = malloc ((n + 1) * sizeof (int));
  p[0] = 0;
  i = n - 1;
  while (p[i + 1] < p[i])
    {
      i = i - 1;
    }
  if (i == 0)
    {
      free (rho);
      rho = NULL;
      return 1; /* undefined */
    }
  j = n;
  while (p[j] < p[i])
    {
      j = j - 1;
    }
  t = p[j];
  p[j] = p[i];
  p[i] = t;
  for (h = i + 1; h < n + 1; h++)
    {
      rho[h] = p[h];
    }
  for (h = i + 1; h < n + 1; h++)
    {
      p[h] = rho[n + i + 1 - h];
    }
  free (rho);
  rho = NULL;
  return 0;
}

int
IMA_permLexRank (int n, int *p)
{
  int r, i, j;
  int *rho;
  
  rho = malloc ((n + 1) * sizeof (int));
  r = 0;
  memcpy (rho, p, (n + 1) * sizeof (int));
  for (j = 1; j < n + 1; j++)
    {
      r = r + (rho[j] - 1) * IMA_factorial (n - j);
      for (i = j + 1; i < n + 1; i++)
        {
          if (rho[i] > rho[j])
            {
              rho[i] = rho[i] - 1;
            }
        }
    }
  free (rho);
  rho = NULL;
  return r;
}

int
IMA_permLexUnrank (int *p, int n, int r)
{
  int i, j, d;
  
  p[n] = 1;
  for (j = 1; j < n; j++)
    {
      d = (r % IMA_factorial(j + 1)) / IMA_factorial(j);
      r = r - d * IMA_factorial(j);
      p[n - j] = d + 1;
      for (i = n - j + 1; i < n + 1; i++)
        {
          if (p[i] > d)
            {
              p[i] = p[i] + 1;
            }
        }
    }
  return 0;
}

int
IMA_factorial (int n)
{
  assert (!(n < 0));
  assert (n < 13);
  if (n == 0)
    {
      return 1;
    }
  if (n == 1)
    {
      return 1;
    }
  if (n == 2)
    {
      return 2;
    }
  else
    {
      return n * IMA_factorial (n - 1);
    }
}

void
IMA_savelocusinfo (int ci, int li)
{
  int ai;
  struct genealogy *G = NULL;
  G = &(C[ci]->G[li]);
  saT[li].savedroot = G->root;
  saT[li].savedmignum = G->mignum;
  saT[li].savedroottime = G->roottime;
  saT[li].savedlength = G->length;
  saT[li].savedtlength = G->tlength;
  saT[li].savedpdg = G->pdg;        
  saT[li].savedplg = G->plg;
  for (ai = 0; ai < L[li].nlinked; ai++)
    {
      saT[li].savedpdg_a[ai] = G->pdg_a[ai];
    }
  return;
}

void
IMA_restorelocusinfo (int ci, int li)
{
  int ai;
  struct genealogy *G = NULL;
  G = &(C[ci]->G[li]);
  G->root = saT[li].savedroot;
  G->mignum = saT[li].savedmignum;
  G->roottime = saT[li].savedroottime;
  G->length = saT[li].savedlength;
  G->tlength = saT[li].savedtlength;
  G->pdg = saT[li].savedpdg;
  G->plg = saT[li].savedplg;
  for (ai = 0; ai < L[li].nlinked; ai++)
    {
      G->pdg_a[ai] = saT[li].savedpdg_a[ai];
    }
  return;
}

void 
IMA_genealogy_save_reset (int li)
{
  saT[li].savedi = 0;
  return;
}

void 
IMA_genealogy_save_init ()
{
  int li;
  int ngnodes;
  struct genealogy *G;

  saT = NULL; 

  init_genealogy_weights (&saC.savedallgweight); 
  init_probcalc (&saC.savedallpcalc);
  saT = malloc (nloci * sizeof (im_savedlocus));
  if (saT == NULL)
    {
      IM_err (IMERR_MEM, "saT");
    }
  for (li = 0; li < nloci; li++)
    {
      G = &C[0]->G[li];
      ngnodes = L[li].numlines;
      saT[li].nlinked = L[li].nlinked;
      saT[li].model = L[li].model;
      saT[li].gtree = NULL;
      saT[li].gtree = malloc (ngnodes * sizeof (im_savedlocus));
      if (saT[li].gtree == NULL)
        {
          IM_err (IMERR_MEM, "sat[%d].gtree", li);
        }
      memset (saT[li].gtree, 0, ngnodes * sizeof (im_savedlocus));
      saT[li].saved = NULL;
      saT[li].savedi = 0;
      saT[li].savedn = 0;
      init_genealogy_weights (&saT[li].savedgweight); 
      saT[li].savedpdg_a = NULL;
      saT[li].savedpdg_a = malloc (L[li].nlinked * sizeof (double));
      if (saT[li].savedpdg_a == NULL)
        {
          IM_err (IMERR_MEM, "saT[%d] savedpdg_a", li);
        }
    }

  return;
}

void 
IMA_genealogy_save_free ()
{
  struct edge *saved = NULL; /* genealogy of locus li */
  int savedi;
  int savedn;
  int si;
  int li;
  int ei;
  int ngnodes;
  im_savedlocus *G;

  for (li = 0; li < nloci; li++)
    {
      G = &saT[li];
      saved = G->saved;
      savedi = G->savedi;
      savedn = G->savedn;

      for (si = 0; si < savedn; si++)
        {
          free (saved[si].mig); 
          saved[si].mig = NULL;
          free (saved[si].A);
          saved[si].A = NULL;
          free (saved[si].dlikeA);
          saved[si].dlikeA = NULL;
        }
      free (saved);
      saved = NULL;
      free (G->savedpdg_a);
      G->savedpdg_a = NULL;

      /* member gtree of im_savedlocus */
      ngnodes = L[li].numlines;
      for (ei = 0; ei < ngnodes; ei++)
        {
          if (G->gtree[ei].cmm > 0)
            {
              free (G->gtree[ei].savedmig);
              G->gtree[ei].savedmig = NULL;
              G->gtree[ei].cmm = 0;
            }
        }
      free (G->gtree);
      G->gtree = NULL;
    }
  free (saT);
  saT = NULL;

  return;
}

int
IMA_genealogy_save_edge (int ci, int li, int ei)
{
  int si;
  int ai;
  int mi;
  int savedi;
  int savedn;
  im_savedlocus *G = NULL;
  struct edge *gtree = NULL;
  struct edge *saved = NULL;

  G = &saT[li];
  gtree = C[ci]->G[li].gtree;
  assert (!(ei < 0));
  assert (ei < L[li].numlines);
  if (!(G->savedi < G->savedn))
    {
      /* allocate more memory of saved */
      G->savedn++;
      G->saved = realloc (G->saved, G->savedn * sizeof(struct edge));
      G->saved[G->savedi].cmm = 0;
      G->saved[G->savedi].mig = NULL;

      G->saved[G->savedi].A = malloc (L[li].nlinked * sizeof (int));
      G->saved[G->savedi].dlikeA = malloc (L[li].nlinked * sizeof (double));

    }
  gtree = C[ci]->G[li].gtree;
  saved = G->saved;
  savedi = G->savedi;
  savedn = G->savedn;

  /* Should we check this loop first?          */
  /* We may use bitset to check the existence. */
  for (si = 0; si < savedi; si++)
    {
      if (saved[si].ei == ei)
        {
          return si;
        }
    }

  saved[savedi].ei = ei;
  saved[savedi].up[0] = gtree[ei].up[0]; 
  saved[savedi].up[1] = gtree[ei].up[1]; 
  saved[savedi].down = gtree[ei].down;
  saved[savedi].time = gtree[ei].time;
  saved[savedi].pop = gtree[ei].pop;

  if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
    {       
      for (ai = (L[li].model == JOINT_IS_SW); 
           ai < L[li].nlinked;
           ai++)
        {
          saved[savedi].A[ai] = gtree[ei].A[ai];
          saved[savedi].dlikeA[ai] = gtree[ei].dlikeA[ai];
        }
    }

  if (modeloptions[NOMIGRATION] == 0)
    {
      if (saved[savedi].cmm < gtree[ei].cmm)
        {
          saved[savedi].cmm = gtree[ei].cmm;
          saved[savedi].mig = realloc (saved[savedi].mig, saved[savedi].cmm *
                                       sizeof (struct migstruct));
        }
      for (mi = 0; mi < gtree[ei].cmm; mi++)
        {
          saved[savedi].mig[mi].mt = gtree[ei].mig[mi].mt;
          saved[savedi].mig[mi].mp = gtree[ei].mig[mi].mp;
          if (gtree[ei].mig[mi].mt < -0.5)
            {
              break;
            }
        }
    }

  G->savedi++;
  return savedi; 
}

void 
IMA_genealogy_restore_edges (int ci, int li)
{
  im_savedlocus *G = NULL;
  struct edge *gtree = NULL; /* genealogy of locus li */
  struct edge *saved = NULL; /* genealogy of locus li */
  int savedi;
  int savedn;
  int mi;
  int si;
  int ei;
  int ai;

  G = &saT[li];
  gtree = C[ci]->G[li].gtree;
  saved = G->saved;
  savedi = G->savedi;
  savedn = G->savedn;

  for (si = 0; si < savedi; si++)
    {
      ei = saved[si].ei;
      assert (!(ei < 0));
      assert (ei < L[li].numlines);
      gtree[ei].up[0] = saved[si].up[0]; 
      gtree[ei].up[1] = saved[si].up[1]; 
      gtree[ei].down = saved[si].down;
      gtree[ei].time = saved[si].time;
      gtree[ei].pop = saved[si].pop;

      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        {       
          for (ai = (L[li].model == JOINT_IS_SW); 
               ai < L[li].nlinked;
               ai++)
            {
              gtree[ei].A[ai] = saved[si].A[ai];
              gtree[ei].dlikeA[ai] = saved[si].dlikeA[ai];
            }
        }

      if (modeloptions[NOMIGRATION] == 0)
        {
          for (mi = 0; mi < saved[si].cmm; mi++)
            {
              assert (mi < gtree[ei].cmm);
              gtree[ei].mig[mi].mt = saved[si].mig[mi].mt;
              gtree[ei].mig[mi].mp = saved[si].mig[mi].mp;
              if (gtree[ei].mig[mi].mt < -0.5)
                {
                  break;
                }
            }
        }
    }
  return;
}

void 
IMA_genealogy_joinsisdown (int ci, int li, int ei)
{
  int i;
  int gparent;
  int gsister;
  int ggrandparent;
  im_savedlocus *S = NULL;
  struct edge *gtree = NULL;
  struct genealogy *G = NULL;    /* locus li */
  int tmrcachange;
  int rootmove;
  int mrootdrop;
  double lmrootdrop;
  /* extend sis, and XFREE up the down edge */
  int sis;
  int j, ai, downdown, down;
  double uptime;

  S = &saT[li];
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  gparent = gtree[ei].down;
  gsister = IMA_genealogy_sister (ci, li, ei);
  assert (!(ei < 0));
  assert (!(gparent < 0));
  assert (!(gsister < 0));
  ggrandparent = gtree[gparent].down;

  /* save four edges and join gsister down to grandparent */
  i = IMA_genealogy_save_edge (ci, li, ei);
  S->copyedge[IM_COPYEDGE_EDGE] = i;
  i = IMA_genealogy_save_edge (ci, li, gparent);
  S->copyedge[IM_COPYEDGE_DOWNEDGE] = i;
  i = IMA_genealogy_save_edge (ci, li, gsister);
  S->copyedge[IM_COPYEDGE_SISEDGE] = i;
  
  if (ggrandparent < 0)
    {
      G->root = gsister;
    }
  else
    {
      IMA_genealogy_save_edge (ci, li, ggrandparent);
    }

  /*********************************************************************** 
  IMA_genealogy_join_edges (ci, li, gsister, gparent, ggrandparent,
                            &tmrcachange,
                            &rootmove,
                            &mrootdrop,
                            &lmrootdrop);
  ***********************************************************************/ 
  sis = gsister;
  gtree = C[ci]->G[li].gtree;
  down = gtree[sis].down;
  i = 0;
  while (gtree[sis].mig[i].mt > -0.5)
    i++;
  j = -1;

  do
    {
      j++;
      checkmigt (i + 1, &gtree[sis]);
      copymig (&gtree[sis].mig[i], &gtree[down].mig[j]);
      i++;
    } while (gtree[down].mig[j].mt > -0.5);
  gtree[sis].time = gtree[down].time;

  /* set the up to which sis now connects */
  gtree[sis].down = gtree[down].down;
  downdown = gtree[sis].down;
  if (downdown != -1)
    {
      rootmove = 0;
      if (gtree[downdown].up[0] == down)
        gtree[downdown].up[0] = sis;

      else
        gtree[downdown].up[1] = sis;

      mrootdrop = 0;
      lmrootdrop = 0.0;
    }
  else
    {
      rootmove = 1;
      tmrcachange += 1;
      /* figure out total time and number of migrants being dropped */
      i = -1;
      do
        {
          i++;
        }
      while (gtree[sis].mig[i].mt > -0.5);
      mrootdrop = i;
      if (sis < L[li].numgenes)
        uptime = 0;
      else
        uptime = gtree[gtree[sis].up[0]].time;

      if (uptime < C[ci]->tvals[lastperiodnumber - 1])
        {
          if (C[ci]->G[li].roottime < C[ci]->tvals[lastperiodnumber - 1])
            lmrootdrop = C[ci]->G[li].roottime - uptime;
          else
            lmrootdrop = C[ci]->tvals[lastperiodnumber - 1] - uptime;
        }
      else
        {
          lmrootdrop = 0.0;
        }
      C[ci]->G[li].root = sis;
      gtree[sis].down = -1;
      gtree[sis].time = TIMEMAX;
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
        for (ai = (L[li].model == JOINT_IS_SW);
             ai < L[li].nlinked; ai++)
          gtree[sis].dlikeA[ai] = 0;
      gtree[sis].mig[0].mt = -1;
      C[ci]->G[li].roottime = uptime;
    }


  /* mrootdrop and lmrootdrop DELETE THESE? */
  S->rootmove = rootmove; 
  S->mrootdrop = mrootdrop;
  S->lmrootdrop = lmrootdrop;

  return;
}

void 
IMA_genealogy_splitsisdown (int ci, int li, int ei, int gparent, int gsister)
{
  /* split newsis into two parts, and make a new down edge out of the lower part */
  int i, j, downdown, nowpop;
  double curt;
  struct edge *gtree; 
  int slidingedge;
  int down;
  int newsis;
  struct genealogy *G = NULL;    /* locus li */
  int ggrandparent;

  G = &(C[ci]->G[li]);
  gtree = G->gtree;

  /* save two edges and split gsister down to grandparent */
  ggrandparent = gtree[gsister].down;
  IMA_genealogy_save_edge (ci, li, gsister);
  if (!(ggrandparent < 0))
    {
      IMA_genealogy_save_edge (ci, li, ggrandparent);
    }

  /* IMA_genealogy_split_edges (ci, li, ei, gsister, gparent, ggrandparent); */
  /* IMA_genealogy_split_edges (int ci, int li, int ei1, int ei2, int ei3, int ei4) */
  /* ei1: ei
   * ei2: gsister,
   * ei3: gparent,
   * ei4: ggrandparent
   */

  slidingedge = ei;
  down = gparent;
  newsis = gsister;
  /* the following four lines used to be all of the function 
   * slidingedge = ei1;
   * down = ei3;
   * newsis = ei2;
   * splitsisdown (ci, li, slidingedge, down, newsis); */

  /* down time set */
  gtree = C[ci]->G[li].gtree;
  curt = gtree[slidingedge].time;
  gtree[down].time = gtree[newsis].time;
  gtree[newsis].time = curt;

  /* set the up of the edge to which down now connects, depends on whether newsis is the root */
  downdown = gtree[newsis].down;
  if (downdown != -1)
  {
    if (gtree[downdown].up[0] == newsis)
      gtree[downdown].up[0] = down;
    else
      gtree[downdown].up[1] = down;
  }
  else
  {
    /* newsis is the current root so the root must move down */
    C[ci]->G[li].root = down;
    C[ci]->G[li].roottime = curt;

    /* C[ci]->G[li].rootmove = 1; FIXME: where is this used? */

    if (C[ci]->G[li].roottime > TIMEMAX)
      {
        IM_err (IMERR_MEM, "Error 72");
      }
    gtree[down].mig[0].mt = -1;
    if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      for (i = (L[li].model == JOINT_IS_SW); i < L[li].nlinked; i++)
        gtree[down].dlikeA[i] = 0;
  }
  gtree[down].down = downdown;

  /* divide the migration events along newsis into upper part for newsis and lower part for down */
  /* this might have bugs setting the population of gtree[down] */
  /* if downdown is negative, i and j must be 0 */
  i = 0;
  while (gtree[newsis].mig[i].mt > -0.5 && gtree[newsis].mig[i].mt < curt)
    i++;
  if (i > 0)
    nowpop = gtree[newsis].mig[i - 1].mp;
  else
    nowpop = gtree[newsis].pop;

  j = findperiodoftime (ci, curt);
  gtree[down].pop = saC.popndown[nowpop][j];
  /********************************
   * j = findperiod (ci, curt);
   * while (C[ci]->poptree[nowpop].e <= j && C[ci]->poptree[nowpop].e != -1)
   *   nowpop = C[ci]->poptree[nowpop].down;
   * gtree[down].pop = nowpop;
  ********************************/
  j = findperiod (ci, curt);
  while (C[ci]->poptree[nowpop].e <= j && C[ci]->poptree[nowpop].e != -1)
    nowpop = C[ci]->poptree[nowpop].down;
  assert (gtree[down].pop == nowpop);
  /* CAN BE DELETED ABOVE 4 LINES */
  nowpop = gtree[down].pop;

  j = 0;
  if (downdown != -1)
  {
    while (gtree[newsis].mig[j + i].mt > -0.5)
    {
      checkmigt (j, &gtree[down]);
      copymig (&gtree[down].mig[j], &gtree[newsis].mig[j + i]);
      /* checkmig (j, &(gtree[down].mig), &(gtree[down].cmm)); */
      /* gtree[down].mig[j] = gtree[newsis].mig[j + i]; */
      assert (nowpop != gtree[down].mig[j].mp);
      nowpop = gtree[down].mig[j].mp;
      j++;
    }
  }
  if (downdown == -1)
    {
      assert (i == 0 && j == 0);
    }
  gtree[newsis].mig[i].mt = -1;
  gtree[down].mig[j].mt = -1;

  gtree[newsis].down = gtree[slidingedge].down = down;
  gtree[down].up[0] = newsis;
  gtree[down].up[1] = slidingedge;





  return;
}

void
IMA_convert_IM2Structurama (char *structurama)
{
  int ii;
  int li;
  int di;
  char *indname;
  int igenotype1;
  int igenotype2;
  FILE *inStructurama;
  inStructurama = NULL;

  inStructurama = fopen (structurama, "w");
  fprintf (inStructurama, "#NEXUS\n\n");
  fprintf (inStructurama, "[ Generated by IMamp ]\n\n");
  fprintf (inStructurama, "begin data;\n");

  fprintf (inStructurama, "\tdimensions nind=%d nloci=%d;\n", snind, nloci);
  fprintf (inStructurama, "\tinfo\n");

  for (ii = 0; ii < snind; ii++)
    {
      indname = IMA_getIndividualName (ii);
      fprintf (inStructurama, "\t%s\t", indname);
      for (li = 0; li < nloci; li++)
        {
          di = IMA_getGenotype (ii, li, &igenotype1, &igenotype2);
          switch (di)
            {
            case 0:
              fprintf (inStructurama, "( ? , ? ) ");
              break;
            case 1: 
              fprintf (inStructurama, "( %d , ? ) ", igenotype1);
              break;
            case 2: 
              fprintf (inStructurama, "( %d , %d ) ", igenotype1, igenotype2);
              break;
            default:
              assert (0);
            }
        }
      fprintf (inStructurama, ",\n");
    }
  fprintf (inStructurama, "\t;\n");
  fprintf (inStructurama, "end;\n\n");

  fclose (inStructurama);
  inStructurama = NULL;

  return;
}

char*
IMA_getIndividualName (int iind)
{
  int ci;
  int li;
  int ei;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL; /* genealogy of locus li */

  ci = 0;
  for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      for (ei = 0; ei < L[li].numgenes; ei++)
        {
          if (gtree[ei].i == iind)
            {
              return L[li].gNames[ei];
            }
        }
    }
  return NULL; 
}

int 
IMA_getGenotype (int ii, int li, int *ig1, int *ig2)
{
  int v;
  int ei;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL; /* genealogy of locus li */

  G = &(C[sci]->G[li]);
  gtree = G->gtree;
  *ig1 = -1;
  *ig2 = -1;
  for (ei = 0; ei < L[li].numgenes; ei++)
    {
      if (gtree[ei].i == ii)
        {
          if (L[li].model == STEPWISE)
            {
              *ig1 = gtree[ei].A[0];
              if (!(L[li].pairs[ei] < 0))
                {
                  *ig2 = gtree[L[li].pairs[ei]].A[0];
                }
            }
          else
            {
              *ig1 = 0;
              if (!(L[li].pairs[ei] < 0))
                {
                  *ig2 = 1;
                }
            }
          break;
        }
    }
  if (*ig1 == -1 && *ig2 == -1)
    {
      v = 0;
    }
  else if (*ig1 != -1 && *ig2 == -1)
    {
      v = 1;
    }
  else if (*ig1 != -1 && *ig2 != -1)
    {
      v = 2;
    }
  else
    {
      assert (0);
    }
  return v;
}

void
IMA_set_popntree_nosplit ()
{
  int ti;
  int ni;
  int pi;
  int i;
  int j;
  int k;
  int npnodes;
  int biti;
  int bitj;
  int n;
  int ip;
  int nperiods;
  int isExist;

  assert (!(npops < 0)); 
  assert (npops < MAXPOPS); 


  npnodes = npops + 1;
  nperiods = 2;


  ssizePi = ((2 * npops - 1) - 1) / 8 + 1;

  assert (lastperiodnumber == 1);

  /* Ps */
  BitPowerSetNew (saC.Ps, ssizePi, nperiods, biti);
  for (ti = 0; ti < lastperiodnumber; ti++)
    {
      for (ni = 0; ni < npops; ni++)
        {
          assert (!(C[sci]->plist[ti][ni] < 0));
          BitTrue (saC.Ps[ti], C[sci]->plist[ti][ni]);
        }
    }
  BitTrue (saC.Ps[1], C[sci]->plist[1][0]);

  /* popndown */
  saC.popndown = malloc ((npops + 1) * sizeof(int *));
  for (pi = 0; pi < npops; pi++)
    {
      saC.popndown[pi] = malloc (nperiods * sizeof(int));
      saC.popndown[pi][0] = pi;
      saC.popndown[pi][1] = npops;
    }
  saC.popndown[npops] = malloc (nperiods * sizeof(int));
  saC.popndown[npops][0] = -1;
  saC.popndown[npops][1] = npops;

  /**********************************************************************/
  /*                           IMA_popnlist_set                         */
  /**********************************************************************/
  assert (!(npops < 0));
  assert (npops < MAXPOPS);

  /* popnlist */
  saC.popnlist = malloc ((npops + 1) * sizeof(gsl_block_uint **));
  saC.popnlist[npops] = malloc (nperiods * sizeof(gsl_block_uint *)); 
  saC.popnlist[npops][1] = gsl_block_uint_alloc (1);
  saC.popnlist[npops][1]->size = 1;
  saC.popnlist[npops][1]->data[0] = npops;
  saC.popnlist[npops][0] = gsl_block_uint_alloc (npops);
  saC.popnlist[npops][0]->size = npops;
  for (i = 0; i < npops; i++)
    {
      saC.popnlist[i] = malloc (nperiods * sizeof(gsl_block_uint *)); 
      saC.popnlist[i][1] = NULL;
      saC.popnlist[i][0] = gsl_block_uint_alloc (1);
      saC.popnlist[i][0]->size = 1;
      saC.popnlist[i][0]->data[0] = i;
      saC.popnlist[npops][0]->data[i] = i;
    }

  /* popnlistB: bit set version of popnlist */
  BitSetPowerSetNew(saC.popnlistB, ssizePi, npnodes, nperiods, biti, bitj);
  for (i = 0; i < npnodes; i++)
    {
      for (j = 0; j < nperiods; j++)
        {
          if (saC.popnlist[i][j] != NULL)
            {
              n = saC.popnlist[i][j]->size;
              for (k = 0; k < n; k++)
                {
                  ip = saC.popnlist[i][j]->data[k];
                  BitTrue (saC.popnlistB[i][j], ip);
                }
            }
        }
    }

  /**********************************************************************/
  /*                        popnmig nosplit                             */
  /**********************************************************************/
  /* numsplittimes: npops - 1 for tree model, and 0 for island model.   */
  /* We could have a variable for number of periods.                    */
  /* numpopsizeparams: 2*npops - 1 for tree model, and npops for island */
  saC.popnmig = malloc ((numsplittimes + 1) * sizeof (int **));
  for (ti = 0; ti < numsplittimes + 1; ti++)
    {
      saC.popnmig[ti] = malloc (numpopsizeparams * sizeof (int *));
      for (pi = 0; pi < numpopsizeparams; pi++)
        {
          /* #_populations to move for period ti, population pi */
          isExist = 0;
          for (ni = 0; ni < npops - ti; ni++)
            {
              if (C[sci]->plist[ti][ni] == pi)
                {
                  isExist = 1;
                }
            }
          if (isExist == 1)
            {
              n = npops - ti - 1;
            }
          else
            {
              n = 0;
            }

          if (n > 0)
            {
              saC.popnmig[ti][pi] = malloc (n * sizeof (int));
              i = 0;
              for (ni = 0; ni < npops - ti; ni++)
                {
                  if (C[sci]->plist[ti][ni] != pi)
                    {
                      saC.popnmig[ti][pi][i] = C[sci]->plist[ti][ni];
                      i++;
                    }
                }
            }
          else
            {
              saC.popnmig[ti][pi] = NULL;
            }
        }
    }

  return; 
}

void
IMA_set_popntree ()
{
  int ti;
  int ni;
  int pi;
  int prevti;
  int i;
  int j;
  int k;
  int isExist;
  int ipExist;
  int sizeExist;
  int npnodes = 2 * npops - 1;
  struct popedge *pnode;
  int biti;
  int bitj;
  int n;
  int ip;

  assert (!(npops < 0)); 
  assert (npops < MAXPOPS); 

  ssizePi = ((2 * npops - 1) - 1) / 8 + 1;

  BitPowerSetNew (saC.Ps, ssizePi, npops, biti);
  for (ti = 0; ti < npops; ti++)
    {
      for (ni = 0; ni < npops - ti; ni++)
        {
          assert (!(C[sci]->plist[ti][ni] < 0));
          BitTrue (saC.Ps[ti], C[sci]->plist[ti][ni]);
        }
    }
  saC.popndown = malloc ((2 * npops - 1) * sizeof(int *));
  for (pi = 0; pi < 2 * npops - 1; pi++)
    {
      saC.popndown[pi] = malloc (npops * sizeof(int));
      for (ti = 0; ti < npops; ti++)
        {
          saC.popndown[pi][ti] = -1;
        }
    }
  for (prevti = npops - 1; prevti > 0; prevti--)
    {
      ti = prevti - 1; 
      saC.popndown[C[sci]->addpop[prevti]][prevti] = C[sci]->addpop[prevti];
      saC.popndown[C[sci]->droppops[prevti][0]][ti] = C[sci]->droppops[prevti][0];
      saC.popndown[C[sci]->droppops[prevti][1]][ti] = C[sci]->droppops[prevti][1];
      for (i = prevti; i < npops; i++)
        {
          saC.popndown[C[sci]->droppops[prevti][0]][i] = 
            saC.popndown[C[sci]->addpop[prevti]][i];
          saC.popndown[C[sci]->droppops[prevti][1]][i] = 
            saC.popndown[C[sci]->addpop[prevti]][i];
        }
      for (i = 0; i < 2 * npops - 1; i++)
        {
          if (i == C[sci]->droppops[prevti][0] || i == C[sci]->droppops[prevti][1])
            {
              continue;
            }
          if (BitIsTrue (saC.Ps[ti], i))
            {
              saC.popndown[i][ti] = saC.popndown[i][prevti];
            }
        }
    }

  /**********************************************************************/
  /*                           IMA_popnlist_set                         */
  /**********************************************************************/
  assert (!(npops < 0));
  assert (npops < MAXPOPS);

  /* popnlist: bit set */
  saC.popnlist = malloc (npnodes * sizeof(gsl_block_uint **));
  for (i = 0; i < npnodes; i++)
    {
      saC.popnlist[i] = malloc (npops * sizeof(gsl_block_uint *)); 
      for (j = 0; j < npops; j++)
        {
          saC.popnlist[i][j] = NULL;
        }
    }
  for (i = 0; i < npnodes; i++)
    {
      pnode = &(C[sci]->poptree[i]); 
      if (pnode->e == -1)
        {
          k = pnode->b + 1;
        }
      else
        {
          k = pnode->e;
        }
      for (j = pnode->b; j < k; j++)
        {
          saC.popnlist[i][j] = gsl_block_uint_alloc (1);
          saC.popnlist[i][j]->data[0] = i;
        }
      for (j = pnode->b - 1; !(j < 0); j--)
        {
          isExist = 0; 
          ipExist = -1;
          sizeExist = saC.popnlist[i][j + 1]->size;
          for (k = 0; k < sizeExist; k++)
            {
              if (saC.popnlist[i][j + 1]->data[k] == C[sci]->addpop[j + 1])
                {
                  isExist = 1;
                  ipExist = k;
                }
            }
          if (isExist == 1)
            {
              saC.popnlist[i][j] = gsl_block_uint_alloc (sizeExist + 1);
            }
          else
            {
              saC.popnlist[i][j] = gsl_block_uint_alloc (sizeExist);
            }
          for (k = 0; k < sizeExist; k++)
            {
              saC.popnlist[i][j]->data[k] = saC.popnlist[i][j + 1]->data[k];
            }
          if (isExist == 1)
            {
              saC.popnlist[i][j]->data[sizeExist] = C[sci]->droppops[j + 1][0];
              saC.popnlist[i][j]->data[ipExist] = C[sci]->droppops[j + 1][1];
            }
        }
    }

  /* popnlistB: bit set version of popnlist */
  BitSetPowerSetNew(saC.popnlistB, ssizePi, npnodes, npops, biti, bitj);
  for (i = 0; i < npnodes; i++)
    {
      for (j = 0; j < npops; j++)
        {
          if (saC.popnlist[i][j] != NULL)
            {
              n = saC.popnlist[i][j]->size;
              for (k = 0; k < n; k++)
                {
                  ip = saC.popnlist[i][j]->data[k];
                  BitTrue (saC.popnlistB[i][j], ip);
                }
            }
        }
    }

  /**********************************************************************/
  /*                       popnmig tree model                           */
  /**********************************************************************/
  /* numsplittimes: npops - 1 for tree model, and 0 for island model.   */
  /* We could have a variable for number of periods.                    */
  /* numpopsizeparams: 2*npops - 1 for tree model, and npops for island */
  saC.popnmig = malloc ((numsplittimes + 1) * sizeof (int **));
  for (ti = 0; ti < lastperiodnumber; ti++)
    {
      saC.popnmig[ti] = malloc (numpopsizeparams * sizeof (int *));
      for (pi = 0; pi < numpopsizeparams; pi++)
        {
          /* #_populations to move for period ti, population pi */
          isExist = 0;
          for (ni = 0; ni < npops - ti; ni++)
            {
              if (C[sci]->plist[ti][ni] == pi)
                {
                  isExist = 1;
                }
            }
          if (isExist == 1)
            {
              n = npops - ti - 1;
            }
          else
            {
              n = 0;
            }

          if (n > 0)
            {
              saC.popnmig[ti][pi] = malloc (n * sizeof (int));
              i = 0;
              for (ni = 0; ni < npops - ti; ni++)
                {
                  if (C[sci]->plist[ti][ni] != pi)
                    {
                      saC.popnmig[ti][pi][i] = C[sci]->plist[ti][ni];
                      i++;
                    }
                }
            }
          else
            {
              saC.popnmig[ti][pi] = NULL;
            }
        }
    }

  return; 
}

void
IMA_unset_popntree ()
{
  int pi;
  int i;
  int j;
  int biti;
  int bitj;
  int npnodes; 
  int nperiods;
  int ti;

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 1)
    {
      npnodes = npops + 1;
      nperiods = 2;
    }
  else
    {
      npnodes = 2 * npops - 1;
      nperiods = npops;
    }

  assert (!(npops < 0)); 
  assert (npops < MAXPOPS); 

  BitPowerSetDelete (saC.Ps, nperiods, biti);
  for (pi = 0; pi < npnodes; pi++)
    {
      free (saC.popndown[pi]);
      saC.popndown[pi] = NULL;
    }
  free (saC.popndown);
  saC.popndown = NULL;

  
  /**********************************************************************/
  /*                           IMA_popnlist_set                         */
  /**********************************************************************/
  for (i = 0; i < npnodes; i++)
    {
      for (j = 0; j < nperiods; j++)
        {
          if (saC.popnlist[i][j] != NULL)
            {
              gsl_block_uint_free (saC.popnlist[i][j]);
              saC.popnlist[i][j] = NULL;
            }
        }
      free (saC.popnlist[i]);
      saC.popnlist[i] = NULL;
    }
  free (saC.popnlist);
  saC.popnlist = NULL;
  BitSetPowerSetDelete(saC.popnlistB, npnodes, nperiods, biti, bitj);

  /**********************************************************************/
  /*                               popnmig                              */
  /**********************************************************************/
  for (ti = 0; ti < lastperiodnumber; ti++)
    {
      for (pi = 0; pi < numpopsizeparams; pi++)
        {
          if (saC.popnmig[ti][pi] != NULL)
            {
              free (saC.popnmig[ti][pi]);
              saC.popnmig[ti][pi] = NULL;
            }
        }
      free (saC.popnmig[ti]);
      saC.popnmig[ti] = NULL;
    }
  free (saC.popnmig);
  saC.popnmig = NULL;

  return;
}

#ifndef NDEBUG
void
assertgenealogyloc (int ci, int li)
{
  if (modeloptions[NOMIGRATION] == 1)
    {
      /* for not individual? 
       * if (assignmentoptions[POPULATIONASSIGNMENTINDIVIDUAL] == 0) 
          assertgenealogyloclabel (ci, li);
          assertgenealogylocdiploid (ci, li);
       */
      if (assignmentoptions[POPULATIONASSIGNMENTLOCAL] == 0)
        {
          assertgenealogyindividual (ci);
        }
    }
  else
    {
      /* if (assignmentoptions[POPULATIONASSIGNMENTINDIVIDUAL] == 0)
          assertgenealogyloclabel (ci, li);
          assertgenealogylocdiploid (ci, li);
       */ 
      if (assignmentoptions[POPULATIONASSIGNMENTLOCAL] == 0)
        {
          assertgenealogyindividual (ci);
        }
    }
  return;
}

void
assertgenealogy (int ci)
{
  int li;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;

  assertgenealogyultrametric (ci);
  assertgenealogylabel (ci);
  assertgenealogybranch (ci);

  for (li = 0; li < nloci; li++)
    {
      assertgenealogylocdiploid (ci, li);
    }

  if (assignmentoptions[POPULATIONASSIGNMENTLOCAL] == 0)
    {
      assertgenealogyindividual (ci);
    }

  for (li = 0; li < nloci; li++)
    {
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      if (G->roottime > TIMEMAX / 10)
        {
          assert (0);
        }
    }

  return;
}

void
assertgenealogyindividual (int ci)
{
  int ii;
  int li;
  int ji;
  int ei;
  int pop1;
  int pop2;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;

  for (ii = 0; ii < snind; ii++)
    {
      ji = 0;
      li = saC.indlist[ii][ji].li;
      ei = saC.indlist[ii][ji].ei;
      G = &(C[ci]->G[li]);
      gtree = G->gtree;
      assert (gtree[ei].up[0] < 0 && gtree[ei].up[1] < 0);
      pop1 = gtree[ei].pop;
      for (ji = 1; ji < sngenes_ind[ii]; ji++)
        {
          li = saC.indlist[ii][ji].li;
          ei = saC.indlist[ii][ji].ei;
          G = &(C[ci]->G[li]);
          gtree = G->gtree;
          assert (gtree[ei].up[0] < 0 && gtree[ei].up[1] < 0);
          pop2 = gtree[ei].pop;
          assert (pop1 == pop2);
        }
    }
  return;
}

void
assertgenealogylocbranch (int ci, int li)
{
  
}

void 
assertgenealogybranch (int ci)
{
  int li;
  for (li = 0; li < nloci; li++)
    {
      assertgenealogylocbranch (ci, li);
    }
  return;
}

void
assertgenealogyedgelabel (int ci, int li, int ei)
{
  /* Check all labels of events along a branch! */
  int left;
  int right;
  int gparent;
  int nmigs;
  int mi;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;
  double ktime;
  double kprevtime;
  int kperiod;
  int kprevperiod;
  int pop;
  int prevpop;
  int ngenes;
  int ngnodes;
 
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  left = gtree[ei].up[0];
  right = gtree[ei].up[1];
  gparent = gtree[ei].down;
  if (modeloptions[NOMIGRATION] == 1)
    {
      nmigs = 0;
    }
  else
    {
      nmigs = IMA_genealogy_nmigration (ci, li, ei);
    }
  assert (!(nmigs < 0));

  if (gparent < 0)
    {
      assert (nmigs == 0);
      kprevperiod = -1;
      prevpop = -1;
    }
  else
    {
      ngenes = L[li].numgenes;
      ngnodes = 2 * ngenes - 1;
      assert (gparent < ngnodes);
      assert (!(gtree[gparent].pop < 0));
      kprevtime = gtree[ei].time;
      kprevperiod = findperiodoftime (ci, kprevtime);
      prevpop = gtree[gparent].pop; 
    }
  for (mi = 0; mi < nmigs; mi++)
    {
      ktime = gtree[ei].mig[nmigs - 1 - mi].mt;
      assert (!(ktime < 0.0));
      kperiod = findperiodoftime (ci, ktime);
      assert (!(prevpop < 0));
      pop = gtree[ei].mig[nmigs - 1 - mi].mp;
      IMA_check_popn_up (kprevperiod, prevpop, kperiod, pop);
      prevpop = pop;
      if (mi == nmigs - 1)
        {
          pop = gtree[ei].pop;
        }
      else
        {
          pop = gtree[ei].mig[nmigs - 2 - mi].mp;
        }
      IMA_check_popn_move (prevpop, kperiod, pop);
      assert (!(pop < 0));
      kprevperiod = kperiod;
      prevpop = pop;
    }
  ktime = IMA_genealogy_time (ci, li, ei);
  kperiod = findperiodoftime (ci, ktime);
  pop = gtree[ei].pop;
  IMA_check_popn_up (kprevperiod, prevpop, kperiod, pop);

  if (left < 0 && right < 0)
    {
      /* No Code */
    }
  else
    {
      assertgenealogyedgelabel (ci, li, left);
      assertgenealogyedgelabel (ci, li, right);
    }
  return;
}

void
assertgenealogylocdiploid (int ci, int li)
{
  int ei;
  int ei2;
  int ngenes;
  struct genealogy *G;
  struct edge *gtree;

  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  ngenes = L[li].numgenes;
  for (ei = 0; ei < ngenes; ei++)
    {
      ei2 = L[li].pairs[ei]; 
      if (ei2 < 0)
        {
          continue;
        }
      else
        {
          assert (gtree[ei].pop == gtree[ei2].pop);
        }
    }
  return;
}

void
assertgenealogylabel (int ci)
{
  int li;

  for (li = 0; li < nloci; li++)
    {
      assertgenealogyloclabel (ci, li);
    }
  return;
}

void
assertgenealogyloclabel (int ci, int li)
{
  struct genealogy *G = NULL;

  G = &(C[ci]->G[li]);
  assertgenealogyedgelabel (ci, li, G->root);

  return;
}

/*
 * New functions
 *
 * assertgenealogy: 
 * We check whether a genealogy is valid. With no migration
 * assumption all internal nodes with the same period are labelled the same
 * population. Tip nodes are labelled the same population as that of their
 * parent nodes if parent nodes are in period 0.
 * Diploid genes must be labelled the same population.
 */

void
assertgenealogyultrametric (int ci)
{
  int li;

  for (li = 0; li < nloci; li++)
    {
      assertgenealogylocultrametric (ci, li);
    }
  return;
}

void
assertgenealogylocultrametric (int ci, int li)
{
  int ngenes;
  int ei;
  int gparent;
  double l;
  double len;
  double lentoroot;
  struct genealogy *G = NULL;
  struct edge *gtree = NULL;
  G = &(C[ci]->G[li]);
  gtree = G->gtree;
  ngenes = L[li].numgenes;


  for (ei = 0; ei < ngenes; ei++)
    {
      gparent = gtree[ei].down;
      len = IMAedge_length (ci, li, ei);
      assert (len > 0.00000000000000001);
      while (gparent != G->root)
        {
          l = IMAedge_length (ci, li, gparent);
          assert (l > 0.00000000000000001);
          len += l;
          gparent = gtree[gparent].down;
        }
      if (ei == 0)
        {
          lentoroot = len;
        }
      assert (fabs(len - lentoroot) < 0.0001);
    }
  return;
}

void
IMA_check_popn_up (int kprevperiod, int prevpop, int kperiod, int pop)
{
  if (!(kprevperiod < 0))
    {
      assert (!(kperiod > kprevperiod));
    }
  assert (!(kperiod < 0));
  assert (kperiod < npops);
  assert (kprevperiod < npops);

  if (kprevperiod < 0) /* a root's population */
    {
      assert (BitIsTrue (saC.Ps[kperiod], pop));
    }
  else /* a non-root's population */
    {
      assert (BitIsTrue (saC.popnlistB[prevpop][kperiod], pop));
    }
  return;
}

void
IMA_check_popn_move (int notpop, int kperiod, int pop)
{
  int ti;
 
  /* [[notpop]] must exist in period [[kperiod]]. */ 
  assert (!(kperiod < 0));
  assert (kperiod < npops);
  assert (BitIsTrue (saC.Ps[kperiod], notpop));
  assert (notpop != pop);

  /* pop does not belong to notpop's descendants. */
  for (ti = kperiod; !(ti < 0); ti--)
    {
      if (BitIsTrue (saC.popnlistB[notpop][ti], pop))
        {
          assert (0);
        }
    }
  
  return;
}
#endif /* NDEBUG */


void
IMA_rearrange_readseq (int li)
{
  /* This function may rearrange the read input sequences.
   * It does nothing, and I leave this function for possible future usage. */
  return;
}

/*
 * Gene names and individuals
 */
void
IMA_ind_init ()
{
  int li;
  int ei;
  int k;
  int ngenes;
  sci = 0;

  saC.genelist = malloc (nloci * sizeof(gsl_block_uint **));
  for (li = 0; li < nloci; li++)
    {
      ngenes = L[li].numgenes;
      saC.genelist[li] = malloc (ngenes * sizeof(gsl_block_uint *));
      for (ei = 0; ei < ngenes; ei++)
        {
          saC.genelist[li][ei] = malloc (nloci * sizeof(gsl_block_uint));
          for (k = 0; k < nloci; k++)
            {
              saC.genelist[li][ei][k].size = 0;
              saC.genelist[li][ei][k].data = NULL;
            }
        }
    }

  sngenes_ind = NULL;

  return; 
}

void
IMA_ind_fin ()
{
  int ii;
  int li;
  int ei;
  int li2;
  int ngenes;
  for (li = 0; li < nloci; li++)
    {
      ngenes = L[li].numgenes;
      for (ei = 0; ei < ngenes; ei++)
        {
          for (li2 = 0; li2 < nloci; li2++)
            {
              if (saC.genelist[li][ei][li2].size > 0)
                {
                  assert (saC.genelist[li][ei][li2].data != NULL);
                  free (saC.genelist[li][ei][li2].data);
                  saC.genelist[li][ei][li2].data = NULL;
                }
            }
          free (saC.genelist[li][ei]);
          saC.genelist[li][ei] = NULL;
        }
      free (saC.genelist[li]);
      saC.genelist[li] = NULL;
    }
  free (saC.genelist);
  saC.genelist = NULL;

  free (sngenes_ind);
  sngenes_ind = NULL;

  for (ii = 0; ii < snind; ii++)
    {
      free (saC.indlist[ii]);
      saC.indlist[ii] = NULL;
    }
  free (saC.indlist);
  saC.indlist = NULL;

  for (ii = 0; ii < snind; ii++)
    {
      free (sind2gi[ii]);
      sind2gi[ii] = NULL;
    } 
  free (sind2gi);
  sind2gi = NULL;

  return; 
}

void
IMA_ind_edge_init_i ()
{
  int ci;
  int li;
  int ei;
  int ngnodes;
  struct genealogy *G = NULL; /* locus li */
  struct edge *gtree = NULL;  /* genealogy of locus li */

  /* initialize */
  for (ci = 0; ci < numchains; ci++)
    {
      for (li = 0; li < nloci; li++)
        {
          G = &C[ci]->G[li];
          gtree = G->gtree;
          ngnodes = L[li].numlines;
          for (ei = 0; ei < ngnodes; ei++)
            {
              gtree[ei].i = -1;
            }
        }
    }
  return;
}

void
IMA_ind_edge_copy_i ()
{
  int ci;
  int li;
  int ei;
  int ngenes;
  struct genealogy *G = NULL;     /* locus li */
  struct edge *gtree = NULL;  /* genealogy of locus li */
  struct genealogy *G2 = NULL;    /* locus li */
  struct edge *gtree2 = NULL; /* genealogy of locus li */

  /* initialize */
  for (ci = 1; ci < numchains; ci++)
    {
      for (li = 0; li < nloci; li++)
        {
          G = &(C[0]->G[li]);
          gtree = G->gtree;
          G2 = &(C[ci]->G[li]);
          gtree2 = G2->gtree;
          ngenes = L[li].numgenes;
          for (ei = 0; ei < ngenes; ei++)
            {
              assert (!(gtree[ei].i < 0));
              assert (gtree2[ei].i == -1);
              gtree2[ei].i = gtree[ei].i;
            }
        }
    }
  return;
}

void
IMA_ind_set ()
{
  int li;
  int li2;
  int ei;
  int ei2;
  int ii;
  int i;
  int j;
  int n;
  int ngenes1;
  int ngenes2;
  char *query;
  int di;
  int *nis;
  int *nii;
  int ngenes;
  struct genealogy *G = NULL;     /* locus li */
  struct edge *gtree = NULL;  /* genealogy of locus li */
  struct genealogy *G2 = NULL;    /* locus li */
  struct edge *gtree2 = NULL; /* genealogy of locus li */
  
/*
  if (assignmentoptions[POPULATIONASSIGNMENTIGNOREDIPLOID] == 1)
    {
      for (li = 0; li < nloci; li++)
        {
          ngenes1 = L[li].numgenes;
          for (ei = 0; ei < ngenes1; ei++)
            {
              query = L[li].gNames[ei];
              for (j = 0; j < ngenes1; j++)
                {
                  if (ei == j)
                    continue;
                  if (!strcmp(query, L[li].gNames[j]))
                    {
                      strcat (L[li].gNames[j], "@");
                      len = strlen (L[li].gNames[j]);
                    }
                }
            }
        }
    }
*/
  for (li = 0; li < nloci; li++)
    {
      ngenes1 = L[li].numgenes;
      for (ei = 0; ei < ngenes1; ei++)
        {
          query = L[li].gNames[ei];
          for (i = 0; i < nloci; i++)
            {
              ngenes2 = L[i].numgenes;
              for (j = 0; j < ngenes2; j++)
                {
                  if (li == i && ei == j)
                    continue;
                  if (!strcmp(query, L[i].gNames[j]))
                    {
                      assert (!(saC.genelist[li][ei][i].size < 0));
                      assert (!(saC.genelist[li][ei][i].size > 2));
                      if (saC.genelist[li][ei][i].size == 0)
                        {
                          assert (saC.genelist[li][ei][i].data == NULL);
                          saC.genelist[li][ei][i].data 
                            = malloc (sizeof(unsigned int));
                          saC.genelist[li][ei][i].size = 1;
                          saC.genelist[li][ei][i].data[0] = j; 
                        }
                      else if (saC.genelist[li][ei][i].size == 1)
                        {
                          saC.genelist[li][ei][i].data 
                            = realloc (saC.genelist[li][ei][i].data, 
                                       2 * sizeof(unsigned int));
                          saC.genelist[li][ei][i].size = 2;
                          saC.genelist[li][ei][i].data[1] = j; 
                        }
                    }
                }
            }
        }
    }

  /* set all i's to -1 */
  IMA_ind_edge_init_i ();

  /* set i's */
  di = 0; 
  sngenes = 0;
  for (li = 0; li < nloci; li++)
    {
      G = &C[sci]->G[li];
      gtree = G->gtree;
      ngenes1 = L[li].numgenes;
      sngenes += ngenes1;
      for (ei = 0; ei < ngenes1; ei++)
        {
          if (gtree[ei].i == -1)
            {
              gtree[ei].i = di;
              for (li2 = 0; li2 < nloci; li2++)
                {
                  n = saC.genelist[li][ei][li2].size;
                  for (i = 0; i < n; i++)
                    {
                      G2 = &(C[sci]->G[li2]);
                      gtree2 = G2->gtree;
                      ei2 = saC.genelist[li][ei][li2].data[i];
                      gtree2[ei2].i = di;
                    }
                }
              di++;
            }
        }
    }
  /* copy all i's of a chain to all the rest of chains  */
  IMA_ind_edge_copy_i ();

  /* total number of individuals */
  snind = di;
  /* list of genes for an individual */
  saC.indlist = malloc (snind * sizeof(im_node *));
  nis = malloc (snind * sizeof(int));
  nii = malloc (snind * sizeof(int));
  sngenes_ind = malloc (snind * sizeof(int)); 
  for (i = 0; i < snind; i++)
    {
      n = IMA_ind_n (i); 
      nis[i] = n;
      nii[i] = 0;
      /* saC.indlist[i] = malloc (n * sizeof(im_node)); */
      saC.indlist[i] = malloc (nloci * 2 * sizeof(im_node));
    }
  for (li = 0; li < nloci; li++)
    {
      G = &(C[sci]->G[li]);
      gtree = G->gtree;
      ngenes1 = L[li].numgenes;
      for (ei = 0; ei < ngenes1; ei++)
        {
          i = gtree[ei].i;
          j = nii[i];
          saC.indlist[i][j].li = li;
          saC.indlist[i][j].ei = ei;
          nii[i]++;
        }
    }
  for (i = 0; i < snind; i++)
    {
      assert (nis[i] == nii[i]);
      sngenes_ind[i] = nis[i];
    }
  sind2gi = malloc (snind * sizeof (int *));
  for (ii = 0; ii < snind; ii++)
    {
      sind2gi[ii] = malloc (nloci * sizeof (int));
      for (li = 0; li < nloci; li++)
        {
          G = &(C[sci]->G[li]);
          gtree = G->gtree;
          ngenes = L[li].numgenes;
          for (ei = 0; ei < ngenes; ei++)
            {
              if (gtree[ei].i == ii)
                {
                  sind2gi[ii][li] = ei;
                  break;
                }
            }
          if (ei == ngenes)
            {
              sind2gi[ii][li] = -1;
            }
        }
    }
 
  free (nis);
  nis = NULL;
  free (nii);
  nii = NULL;

  /* just print indlist */
  for (i = 0; i < snind; i++)
    {
      printf ("ind [%3d]\t", i);
      for (j = 0; j < sngenes_ind[i]; j++)
        {
          li = saC.indlist[i][j].li;
          ei = saC.indlist[i][j].ei;
          printf ("[%d](%d-%d)\t", j, li, ei);
        }
      printf ("\n");
    }
 
  return;
}

int
IMA_ind_n (int k)
{
  struct genealogy *G = NULL;     /* locus li */
  struct edge *gtree = NULL;  /* genealogy of locus li */
  int li;
  int ei;
  int n;
  int i;
  int ngenes;

  n = 0;
  for (li = 0; li < nloci; li++)
    {
      G = &(C[sci]->G[li]);
      gtree = G->gtree;
      ngenes = L[li].numgenes;
      for (ei = 0; ei < ngenes; ei++)
        {
          i = gtree[ei].i;
          if (i == k)
            {
              n++;
            }
        }
    }
  return n;
}

int
IMA_ind_find (int **v, int l1, int l2, ...)
{
  int i;
  int j;
  int n;
  int ei;
  va_list ap;

  n = 0;
  va_start  (ap, l2);
  ei = va_arg (ap, int);
  while (ei != -1)
    {
      n += saC.genelist[l1][ei][l2].size;
      ei = va_arg (ap, int);
    }
  va_end (ap);
  
  if (n > 0)
    {
      *v = malloc (n * sizeof(int));
      va_start  (ap, l2);
      for (i = 0; i < n; i++)
        {
          ei = va_arg (ap, int);
          for (j = 0; j < saC.genelist[l1][ei][l2].size; j++)
            {
              (*v)[i] = saC.genelist[l1][ei][l2].data[j];
              i++;
            }
        }
      va_end (ap);
    }
  else
    {
      *v = NULL;
    }
  return n;
}


int 
IMA_rgs_convert (int *a, int n)
{
  int max;
  int i;
  int *m;
  int e;

  m = malloc (n * sizeof (int));
  memset (m, -1, n * sizeof (int)); 

  max = 1;
  for (i = 0; i < n; i++)
    {
      e = a[i];
      assert (e < n);
      assert (!(e < 0));
      if (m[e] == -1)
        {
          m[e] = max;
          max++;
        }
      a[i] = m[e];
    }
  free (m);
  m = NULL;

  return max;

}

int 
IMA_rgf_setSnasn (int nind, int np)
{
  int v;
  int i;

  if (nind < np)
    {
      IM_err (IMERR_ASN, "There are more populations than individuals");
    }
  v = 0;
  for (i = 0 ; i < np; i++)
    {
      v += IMA_stirlings2 (nind, i + 1);
    }
  return v;
}

void 
IMA_rgf_tickSaasn (int ci)
{
  int li;
  int ii;
  int ei;
  int *a = NULL;
  int *pa;
  int r;

  assert (ci == 0);

  a = malloc (snind * sizeof (int));

  /* how do I know individuals of a locus? */
  pa = a;
  for (ii = 0; ii < snind; ii++)
    {
      li = saC.indlist[ii][0].li;
      ei = saC.indlist[ii][0].ei;
      *pa = C[ci]->G[li].gtree[ei].pop;
      pa++;
    }

  IMA_rgs_convert (a, snind);

  /* Add the number of assignments. */
  r = IMA_rankrgf2 (snind, npops, a);
  assert (r < snasn);
  saasn[r]++;

  free (a);
  a = NULL;
 
  return;
}

void 
IMA_rgf_saveSaasn ()
{
  int i;
  FILE *tfile;

  tfile = fopen ("a.txt", "w");
  for (i = 0; i < snasn; i++)
    {
      fprintf (tfile, "%d\n", saasn[i]);
    }
  fclose (tfile);
  tfile = NULL;
  return;
}

int
IMA_fprintrgf (FILE *fp, int *f, int m)
{
  int i;
  fprintf (fp, "{");
  for (i = 0; i < m; i++)
    {
      if (i < m - 1)
        {
          fprintf (fp, "%2d ", f[i]);                
        }
      else
        {
          fprintf (fp, "%2d}", f[i]);
        }
    }
  fprintf (fp, "\n");
  return 0;
}

int
IMA_generatergf (int m, int n)
{
  int *f;
  int *fm;
  int i;
  int j;
  char done;
  int maxn;
  int rank;
  
  rank = 0;
  f = malloc (m * sizeof (int));
  fm = malloc (m * sizeof (int));
  for (i = 0; i < m; i++)
    {
      f[i] = 1;
      fm[i] = 2;
    }
  done = 'N';
  while (done == 'N')
    {
      fprintf (stdout, "[%3d] ", rank);
      rank++;
      IMA_fprintrgf (stdout, f, m);          
      j = m;
      do
        {
          j = j - 1;
        } while (f[j] == fm[j]);
      if (j > 0)
        {
          f[j] = f[j] + 1;
          for (i = j + 1; i < m; i++)
            {
              f[i] = 1;
              if (f[j] == fm[j])
                {
                  if (fm[j] < n - 1)
                    {
                      maxn = fm[j];
                    }
                  else
                    {
                      maxn = n - 1;
                    }
                  fm[i] = maxn + 1;
                }
              else
                {
                  fm[i] = fm[j];
                }
            }
        }
      else
        {
          done = 'Y';
        }
    }
  
  free (f);
  free (fm);
  f = NULL;
  fm = NULL;
  return 0;
}

int**
IMA_generalizedrgf (int m)
{
  return IMA_generalizedrgf2 (m, m);
}

int**
IMA_generalizedrgf2 (int m, int n)
{
  int j;
  int i;
  int **d;
  int M;
  int N;
  
  M = m + 1;
  N = n + 1;
  d = malloc (M * sizeof (int *));
  for (i = 0; i < M; i++)
    {
      d[i] = malloc (M * sizeof (int));
      memset (d[i], 0, M * sizeof (int));
    }
  for (j = 0; j < N; j++)
    {
      d[0][j] = 1;
    }
  for (i = 1; i < M; i++)
    {
      for (j = 0; j < M - i; j++)
        {
          d[i][j] = j * d[i-1][j] + d[i-1][j+1];
        }
    }
  return d;
}

int
IMA_freejrgf (int ***d, int m)
{
  int i;
  int M;
  
  M = m + 1;
  for (i = 0; i < M; i++)
    {
      free ((*d)[i]);
      (*d)[i] = NULL;
    }
  free (*d);
  *d = NULL;
  return 0;
}

int
IMA_fprintjrgf (FILE *fp, int **d, int m)
{
  int i;
  int j;
  int M;
  
  M = m + 1;

  fprintf (fp, "  i | ");
  for (j = 0; j < M; j++)
    {
      fprintf (fp, "%4d ", j);
    }
  fprintf (fp, "\n");
  fprintf (fp, "-----");
  for (j = 0; j < M; j++)
    {
      fprintf (fp, "-----");
    }  
  fprintf (fp, "\n");
  
  for (i = 0; i < M; i++)
    {
      fprintf (fp, "%3d | ", i);
      for (j = 0; j < M - i; j++)
        {
          fprintf (fp, "%4d ", d[i][j]);
        }
      fprintf (fp, "\n");
    } 
  return 0;
}

int
IMA_rankrgf (int m, int *f)
{
  return IMA_rankrgf2 (m, m, f);
}

int
IMA_rankrgf2 (int m, int n, int *f)
{
  int r;
  int j;
  int i;
  int **d;
  int M;
  
  for (i = 0; i < m; i++)
    {
      if (f[i] > n)
        {
          printf ("INVALID RGFs - There are elements larger than %d\n", n);
          return -1;
        }
    }
  M = m + 1;
  d = IMA_generalizedrgf2 (m, n);
  r = 0;
  j = 1;
  for (i = 2; i < M; i++)
    {
      r = r + (f[i-1] - 1) * d[m-i][j];
      if (j < f[i-1])
        {
          j = f[i-1];
        }
    }
  IMA_freejrgf (&d, m);  
  return r;
}

int
IMA_unrankrgf (int *f, int r, int m)
{  
  IMA_unrankrgf2 (f, r, m, m);
  return 0;
}

int
IMA_unrankrgf2 (int *f, int r, int m, int n)
{
  int j;
  int i;
  int **d;
  int M;
  
  M = m + 1;
  d = IMA_generalizedrgf2 (m, n);
  j = 1;
  for (i = 2; i < M; i++)
    {
      if (j * d[m-i][j] <= r)
        {
          f[i-1] = j + 1;
          r = r - j * d[m-i][j];
          j = j + 1;
        }
      else
        {
          f[i-1] = r / d[m-i][j] + 1;
          r = r % d[m-i][j];
        }
    }
  IMA_freejrgf (&d, m);
  return 0;
}

int
IMA_stirlings2 (int m, int n)
{
  int i;
  int j;
  int M;
  int N;
  int min;
  int **s;
  int v;
  
  M = m + 1;
  N = n + 1;
  
  s = malloc (M * sizeof (int *));
  for (i = 0; i < M; i++)
    {
      s[i] = malloc ((M + 1) * sizeof (int));
    }
  
  s[0][0] = 1;
  for (i = 1; i < M; i++)
    {
      s[i][0] = 0;
    }
  for (i = 0; i < M; i++)
    {
      s[i][i+1] = 0;
    }
  for (i = 1; i < M; i++)
    {
      if (i < n)
        {
          min = i;
        }
      else
        {
          min = n;
        }
      for (j = 1; j < min + 1; j++)
        {
          s[i][j] = j * s[i-1][j] + s[i-1][j-1];
        }
    }
  v = s[m][n];
  
  for (i = 0; i < M; i++)
    {
      free (s[i]);
      s[i] = NULL;
    }
  free (s);
  s = NULL;
  
  return v;
}

int
IMA_maxrgf (int *f, int m)
{
  int i;
  int maxn;
  maxn = 0;
  for (i = 0; i < m; i++)
    {
      if (maxn < f[i])
        {
          maxn = f[i];
        }
    }
  return maxn;
}


void 
IMA_align_genealogy (char *fname, int n)
{
  int i;
  int j;
  int k;
  int l;
  int nsispop;
  int **sispop;
  int *od;
  int *od1;
  char *pname;
  int len;
  int *tasn;
  int *easn;
  int nrow;
  int ncol;
  int **m;
  int ncombsis;
  char cb[5];
  int d;
  int d2;
  unsigned char uc;
  int v;

  if (1)
    {
      IMA_structurama_run (fname);
    }
  else
    {
      IMA_structurama_localrun (fname);
    }
  if (trueassignment != NULL)
    {
      IMA_structurama_summary (fname);
    }
  IMA_structurama_uncertainty (fname);

  /* Align genealogies by swapping elements in [[gsampinf]]. */
  len = strlen (fname);
  pname = (char *) malloc ((len + 30) * sizeof (char)); /* .in.sum_assignments */

  /* Find a point estimate of assignment parameter. */
  tasn = malloc (snind * sizeof (int));
  sprintf (pname, "%s.in.sum_assignments", fname);
  IMA_io_readsumassignment (pname, snind, tasn);
  for (i = 0; i < snind; i++)
    {
      assert (tasn[i] > 0); 
      tasn[i] = tasn[i] - 1;
    }

  /* Read all posterior samples of assignment parameter. */
  sprintf (pname, "%s.asn", fname);
  IMA_io_readp (pname, &m, &nrow, &ncol);
  /* This may be different when there are multiple MCMC ouput files. */
  assert (ncol == snind);
  assert (nrow == n);  

  if (assignmentoptions[POPULATIONASSIGNMENTINFINITE] == 0)
    {
      assert (npops == 2);
      easn = malloc (snind * sizeof (int));
      nsispop = IMA_ptree_nsispop (0);
      if (nsispop < 1 || nsispop > 5)
        {
          IM_err (IMERR_ASN, 
                  " number (%d) of present-day sister populations is in between 1 and 5",
                  nsispop);
        }
      sispop = (int **) malloc (nsispop * sizeof (int *));
      for (j = 0; j < nsispop; j++)
        {
          sispop[j] = (int *) malloc (2 * sizeof (int));
        }
      IMA_ptree_sispop (0, sispop);
      ncombsis = 1 << nsispop;

      /* Compare two assignments. */
      for (i = 0; i < n; i++)
        {
          /* 2^^#sisters many iterations */
          d = snind;
          for (j = 0; j < ncombsis; j++)
            {
              uc = (unsigned char) j;
              for (k = 0; k < 5; k++)
                {
                  if (uc & 1<<k)
                    {
                      cb[k] = '1';
                    }
                  else
                    {
                      cb[k] = '0';
                    }
                }

              /* Find the number of individuals to be removed. */
              memcpy (easn, m[i], snind * sizeof (int));
              for (k = 0; k < 5; k++)
                {
                  if (cb[k] == '1')
                    {
                      for (l = 0; l < snind; l++)
                        {
                          if (easn[l] == sispop[k][0]) 
                            {
                              easn[l] = sispop[k][1];
                            }
                          else if (easn[l] == sispop[k][1])
                            {
                              easn[l] = sispop[k][0];
                            }
                        }
                      
                    }
                }
              d2 = IMA_gdistanceSimple (snind, tasn, easn);

              if (d > d2)
                {
                  d = d2;
                  v = j;
                }
            }

          /* v-th permutation. */
          uc = (unsigned char) v;
          for (k = 0; k < 5; k++)
            {
              if (uc & 1<<k)
                {
                  cb[k] = '1';
                }
              else
                {
                  cb[k] = '0';
                }
            }
          for (j = 0; j < nsispop; j++)
            {
              if (cb[j] == '1')
                {
                  IMA_gsampinf_swap (i, sispop[j][0], sispop[j][1]);
                }
            }
        }
      XFREE (easn);
      for (i = 0; i < nsispop; i++)
        {
          XFREE (sispop[i]);
        }
      XFREE (sispop);
    }
  else
    {
      od = (int *) malloc (npops * sizeof (int));
      od1 = (int *) malloc (npops * sizeof (int));
      for (i = 0; i < n; i++)
        {
          for (j = 0; j < npops; j++)
            {
              od1[j] = j;
            }
          IMA_gdistanceOrder (snind, tasn, m[i], npops, od);
          for (j = 0; j < npops; j++)
            {
              if (od[j] != od1[j])
                {
                  /* Check if the swap is allowed under a tree model. */
                  IMA_gsampinf_swap (i, od[j], od1[j]);
                  l = -1;
                  for (k = j + 1; k < npops; k++)
                    {
                      if (od1[k] == od[j])
                        {
                          l = k;
                          break;
                        }
                    }
                  assert (l > 0); 
                  SWAP (od1[j], od1[l]);
                }
            }
          for (j = 0; j < npops; j++)
            {
              assert (od[j] == od1[j]);
            }
        }
      XFREE (od);
      XFREE (od1);
    }

  XFREE (pname);
  XFREE (tasn);
  for (i = 0; i < nrow; i++)
    {
      XFREE (m[i]);
    }
  XFREE (m);

  return; 
}

void 
IMA_gsampinf_swap (int gi, int pi, int pj)
{
  int i, j;
  int mi;
  char migstr[10]; /* m0>1 or q0 */

  /* All vlaues in gsampinf associated labels pi and pj are swapped. */
  assert (npops == 2);
  assert (pi < npops && pj < npops);
  assert (!(pi < 0) && !(pj < 0));

  i = pi;
  j = pj;
  SWAP (gsampinf[gi][gsamp_ccp + i], gsampinf[gi][gsamp_ccp + j]); 
  SWAP (gsampinf[gi][gsamp_fcp + i], gsampinf[gi][gsamp_fcp + j]); 
  SWAP (gsampinf[gi][gsamp_hccp + i], gsampinf[gi][gsamp_hccp + j]); 
  SWAP (gsampinf[gi][gsamp_qip + i], gsampinf[gi][gsamp_qip + j]); 

  i = nummigrateparams;
  sprintf (migstr, "m%d>%d", pi, pj);
  assert (strlen (migstr) == 4);
  for (mi = 0; mi < nummigrateparams; mi++)
    {
      if (strcmp(migstr, imig[mi].str) == 0)
        {
          i = mi;
        }
    }
  assert (i < nummigrateparams);

  j = nummigrateparams;
  sprintf (migstr, "m%d>%d", pj, pi);
  assert (strlen (migstr) == 4);
  for (mi = 0; mi < nummigrateparams; mi++)
    {
      if (strcmp(migstr, imig[mi].str) == 0)
        {
          j = mi;
        }
    }
  assert (j < nummigrateparams);
  SWAP (gsampinf[gi][gsamp_mcp + i], gsampinf[gi][gsamp_mcp + j]); 
  SWAP (gsampinf[gi][gsamp_fmp + i], gsampinf[gi][gsamp_fmp + j]); 
  SWAP (gsampinf[gi][gsamp_mip + i], gsampinf[gi][gsamp_mip + j]); 

  return;
}

int 
IMA_ptree_nsispop (int ci)
{
  int i;
  int v;
  struct popedge *ptree;

  ptree = C[ci]->poptree;
  v = 0;
  for (i = npops; i < 2 * npops - 1; i++)
    {
      if (ptree[i].up[0] < npops && ptree[i].up[1] < npops) 
        {
          v++;
        }
    }
  return v;
}

int 
IMA_ptree_sispop (int ci, int **s)
{
  int i;
  int v;
  struct popedge *ptree;

  ptree = C[ci]->poptree;
  v = 0;
  for (i = npops; i < 2 * npops - 1; i++)
    {
      if (ptree[i].up[0] < npops && ptree[i].up[1] < npops) 
        {
          s[v][0] = ptree[i].up[0];
          s[v][1] = ptree[i].up[1];
          v++;
        }
    }
  return v;
}



void 
IMA_output_structurama_bat (char *fname)
{
  FILE *fp;
  char *batname;
  char *n;
  int lenbatname;

  fp = NULL;
  batname = NULL;

  lenbatname = strlen (fname) + 4 + 1;
  batname = malloc (lenbatname * sizeof (char));
  sprintf (batname, "%s.bat", fname);
  fp = fopen (batname, "w");
  n = strrchr (fname, '/');
  if (n == NULL)
    {
      /* fprintf (fp, "execute %s.in\n", fname); */
    }
  else
    {
      /* fprintf (fp, "execute %s.in\n", n + 1); */
    }
  fprintf (fp, "execute 1.in\n");
  fprintf (fp, "sum\n");
  fprintf (fp, "quit\n");

  free (batname);
  fclose (fp);
  fp = NULL;
  batname = NULL;
  return;
}

void 
IMA_structurama_localrun (char *fname)
{
  int lencmd;
  char *cmd;

  /* Send input files for STRUCTURAMA. */
  /* Make a system command for sending two input files for STRUCTURAMA. */
  lencmd = 65 + strlen (fname);

  cmd = malloc (lencmd * sizeof (char));
  sprintf (cmd, "cp %s.in 1.in", fname);
  system (cmd);
  sprintf (cmd, "cp %s.in.p 1.in.p", fname);
  system (cmd);
  sprintf (cmd, "cp %s.bat 1.bat", fname);
  system (cmd);

  /* Run STRUCTURAMA. */
  sprintf (cmd, "structurama_1.0 < 1.bat >> log.txt");
  system (cmd);

  sprintf (cmd, "mv 1.in.sum_assignments %s.in.sum_assignments", fname);
  system (cmd);
  sprintf (cmd, "mv 1.in.sum_dist.nex %s.in.sum_dist.nex", fname);
  system (cmd);
  sprintf (cmd, "mv 1.in.sum_pairs %s.in.sum_pairs", fname);
  system (cmd);

  free (cmd);
  cmd = NULL;
  return;
}


/*********************************************************************/
/*                      TCP/IP not for Windows                       */
/*********************************************************************/
#ifndef _MSC_VER

/**********************************************************************/
/*                                                                    */
/* The client part of a client-server pair. This simply takes two     */
/* numbers and adds them together, returning the result to the client */
/*                                                                    */
/* User types:                                                        */
/*                   3 + 5                                            */
/*                   a + b                                            */
/*                   halt + server                                    */
/**********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <unistd.h>

/* #define PORT 0x1234          define a port number for the service */
#define PORT "3490"  /* define a port number for the service */
/* #define HOST "macbook" */
#define HOST "desktop4"
/* #define HOST "localhost" */
#define bufsize 20


int sendall(int s, char *buf, int *len)
{
  int total = 0;        // how many bytes we've sent
  int bytesleft = *len; // how many we have left to send
  int n;

  while(total < *len) 
    {
      n = send (s, buf+total, (size_t) bytesleft, 0);
      if (n == -1) { break; }
      total += n;
      bytesleft -= n;
    }

  *len = total; // return number actually sent here
  return n==-1?-1:0; // return -1 on failure, 0 on success
} 

int recvall(int s, char *buf, int *len)
{
  int total = 0;        // how many bytes we've receive
  int bytesleft = *len; // how many we have left to receive
  int n;

  while(total < *len) 
    {
      n = recv (s, buf+total, (size_t) bytesleft, 0);
      if (n == -1) { break; }
      total += n;
      bytesleft -= n;
    }

  *len = total; // return number actually sent here
  return n==-1?-1:0; // return -1 on failure, 0 on success
} 

int
IMA_tcp_send (int sd, char *cmd)
{
  int len;
  int lenbuf;
  char buf[bufsize];
  FILE *fp;
  int i;

  fp = fopen (cmd, "r");
  lenbuf = bufsize;
  while (lenbuf == bufsize)
    {
      lenbuf = 0;
      for (i = 0; i < bufsize; i++)
        {
          if (feof (fp))
            {
              break;
            }
          buf[i] = fgetc (fp);
          lenbuf++;
        }
      /* Send the length of buffer to send. */
      if (lenbuf != bufsize)
        {
          lenbuf--;
        }
      if (lenbuf > 0)
        {
          /* printf ("Sending: %d\n", lenbuf); */
          len = sizeof (int);
          if (sendall (sd, (char *)&lenbuf, &len) == -1)
            {
              perror ("send");
              exit (1);
            }

          if (sendall (sd, buf, &lenbuf) == -1)
            {
              perror ("send");
              exit (1);
            }
        }
    }
  len = sizeof (int);
  lenbuf = 999;
  if (sendall (sd, (char *)&lenbuf, &len) == -1)
    {
      perror ("send");
      exit (1);
    }


  fclose (fp);
  fp = NULL;
  printf ("Sent: %s\n", cmd);
  return 0;
}


int
IMA_tcp_recv (int sd, char *cmd)
{
  int len;
  int lenbuf;
  char buf[bufsize];
  FILE *fp;
  int i;

  fp = fopen (cmd, "w");
  lenbuf = bufsize;
  while (lenbuf == bufsize)
    {
      /* Receive the length of buffer. */
      len = sizeof (int);
      if (recvall (sd, (char *)&lenbuf, &len) == -1)
        {
          perror ("recv");
          exit (1);
        }
      if (lenbuf == 999)
        {
          break;
        }
      /* printf ("Received: %d\n", lenbuf); */
      if (recvall (sd, buf, &lenbuf) == -1)
        {
          perror ("recv");
          exit (1);
        }
      for (i = 0; i < lenbuf; i++)
        {
          fputc (buf[i], fp);
        }
    }
  if (lenbuf < bufsize)
    {
      len = sizeof (int);
      if (recvall (sd, (char *)&lenbuf, &len) == -1)
        {
          perror ("recv");
          exit (1);
        }
      assert (lenbuf == 999);
    }
  fclose (fp);
  fp = NULL;
  printf ("Received: %s\n", cmd);
  return 0;

}


// get sockaddr, IPv4 or IPv6:
void *
get_in_addr (struct sockaddr *sa)
{
  if (sa->sa_family == AF_INET)
    {
      return &(((struct sockaddr_in *) sa)->sin_addr);
    }

  return &(((struct sockaddr_in6 *) sa)->sin6_addr);
}

void 
IMA_structurama_run (char *fname)
{
  int lencmd;
  char *cmd;
  char *bname;
  int len;
  int lenbname;

  int sockfd;
  struct addrinfo hints, *servinfo, *p;
  int rv;
  char s[INET6_ADDRSTRLEN];

  memset (&hints, 0, sizeof hints);
  hints.ai_family = AF_UNSPEC;
  hints.ai_socktype = SOCK_STREAM;

  /* if ((rv = getaddrinfo ("172.17.34.102", PORT, &hints, &servinfo)) != 0) for Hey Lab ``cluster'' */
  /* if ((rv = getaddrinfo ("10.0.0.1", PORT, &hints, &servinfo)) != 0) for Genetics cluster. */
  if ((rv = getaddrinfo ("172.17.34.101", PORT, &hints, &servinfo)) != 0)
    {
      fprintf (stderr, "getaddrinfo: %s\n", gai_strerror (rv));
      return;
    }

  // loop through all the results and connect to the first we can
  for (p = servinfo; p != NULL; p = p->ai_next)
    {
      if ((sockfd = socket (p->ai_family, p->ai_socktype,
                            p->ai_protocol)) == -1)
        {
          perror ("client: socket");
          continue;
        }

      if (connect (sockfd, p->ai_addr, p->ai_addrlen) == -1)
        {
          close (sockfd);
          perror ("client: connect");
          continue;
        }

      break;
    }

  if (p == NULL)
    {
      fprintf (stderr, "client: failed to connect\n");
      return;
    }

  inet_ntop (p->ai_family, get_in_addr ((struct sockaddr *) p->ai_addr),
             s, sizeof s);
  printf ("client: connecting to %s\n", s);

  freeaddrinfo (servinfo);      // all done with this structure

  lencmd = 65 + strlen (fname);

  bname = strrchr (fname, '/');
  if (bname == NULL)
    {
      bname = fname;
    }
  else
    {
      bname = bname + 1;
    }
  lenbname = strlen (bname);
  len = sizeof (int);
  if (sendall (sockfd, (char *)&lenbname, &len) == -1)
    {
      perror ("send");
      exit (1);
    }

  if (sendall (sockfd, bname, &lenbname) == -1)
    {
      perror ("send");
      exit (1);
    }


  cmd = malloc (lencmd * sizeof (char));
  sprintf (cmd, "%s.in", fname);
  IMA_tcp_send (sockfd, cmd); 
  sprintf (cmd, "%s.in.p", fname);
  IMA_tcp_send (sockfd, cmd); 
/*
  sprintf (cmd, "%s.bat", fname);
  IMA_tcp_send (sockfd, cmd); 
*/

  sprintf (cmd, "%s.in.sum_assignments", fname);
  IMA_tcp_recv (sockfd, cmd); 
  sprintf (cmd, "%s.in.sum_dist.nex", fname);
  IMA_tcp_recv (sockfd, cmd); 
  sprintf (cmd, "%s.in.sum_pairs", fname);
  IMA_tcp_recv (sockfd, cmd); 

  printf ("Finished!\n");

  free (cmd);
  cmd = NULL;

  close (sockfd);

}
#undef PORT
#undef HOST 
#undef bufsize

#else
void 
IMA_structurama_run (char *fname)
{
  IMA_structurama_localrun (fname);
}
#endif /* _MSC_VER */

/*********************************************************************/
/*                      TCP/IP not for Windows                       */
/*********************************************************************/





