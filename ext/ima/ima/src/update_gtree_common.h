/* IMa  2007-2009  Jody Hey, Rasmus Nielsen and Sang Chul Choi*/
#ifndef _UPDATE_GTREE_COMMON_H_
#define _UPDATE_GTREE_COMMON_H_
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

__BEGIN_DECLS
/* updating ancestral alleles at nodes,  under stepwise model */
#define updateAfrac 0.05        // inverse of updateAi, just a useful proportion,  don't need to do all of the nodes every step
#define updateAi  20
void init_update_assignment (void);
void free_update_assignment (void);

/* prototypes */



double likelihoodDG (int ci, int li);
int findperiod_Assignment (int ci, double t);
void treeweight_Assignment (int ci, int li);
void integrate_tree_prob_Assignment (int ci,
                                     struct genealogy_weights *gweight,
                                     struct probcalc *pcalc);
void copyfraclike_Assignment (int ci, int li);
void storescalefactors_Assignment (int ci, int li);
void restorescalefactors_Assignment (int ci, int li);
double finishSWupdateA_Assignment (int ci, int li, int ai, int edge,
                                   int downedge, int sisedge, int newsisedge,
                                   double u, double *Aterm);
double updateA_Assignment (int ci, int li, int ai, double u, int *count);
void init_gtreecommon (void);
void free_gtreecommon (void);
double calcmrate (int mc, double mt);
void joinsisdown (int ci, int li, int sis, int *tmrcachange);
void splitsisdown (int ci, int li, int slidingedge, int down, int newsis);
void getm (int ci, double mrate, struct edgemiginfo *edgem,
           struct edgemiginfo *sisem);
void slider (int ci, int li, int slidingedge, int *sis, double *timepoint,
             double *slidedist);
void slider_nomigration (int ci, int li, int slidingedge, int *sis,
                         double *timepoint, double *slidedist);
void IMA_reset_edgemiginfo (struct edgemiginfo *em);
void storeoldedges (int ci, int li, int edge, int sisedge, int downedge);
void restoreedges (int ci, int li, int edge, int sisedge, int downedge,
                   int newsisedge);
double getmprob (int ci, double mrate, struct edgemiginfo *edgem,
                 struct edgemiginfo *sisem);
void storeAinfo (int ci, int li, struct edge *gtree, int edge, int sisedge,
                 int downedge);
void fillmiginfoperiods (int ci, struct edgemiginfo *em);
void fillmiginfo (int ci, int li, struct edge *gtree, int edge, int sisedge);
void copynewmig_to_gtree (int ci, int li);
void storetreestats (int ci, int li, int mode);
int picktopop (int nowpop, int plist[], int numpops);
int picktopop2 (int nowpop, int plist[], int numpops, int notother);
int mwork (int ci, struct edgemiginfo *edgem, int lastmigperiod,
           double mrate);
double pathcondition (int issame, int moves, int popsless1);

double getnewt (int timeperiod, double t_u_prior, double t_d_prior,
                double oldt, int whichupdate);


__END_DECLS
#endif /* _UPDATE_GTREE_COMMON_H_ */
