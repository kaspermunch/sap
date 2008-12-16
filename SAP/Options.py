
import re, time, os, sys
from optparse import OptionParser

class Options:

    def __init__(self):

        usage="""%prog [options] <input files>

The program does statistical assignment of DNA sequences to taxonomic
groups represented in an annnotated sequence database.

Type 'sap --onlinehelp' to open the online manual in your default browser."""

        self.parser = OptionParser(usage=usage, version="%prog 1.0")

        # General options:
        self.parser.add_option("--onlinehelp",
                          action="store_true",
                          default=False,
                          help="Open online manual in default browser.")
        self.parser.add_option("-d", "--project",
                          type="string",
                          default = os.path.abspath(time.strftime('project-%H.%M-%m.%d.%y')),
                          #default=None,
                          help="Output directory for this project. Results will be written to this directory and cached results will be read from and written to this directory.")
        self.parser.add_option("-D", "--database",
                          type="string",
                          default = "GenBank",
                          help="Name of database plugin to use. Default is online GenBank. You can also specify the path to a file in FASTA format to serve as database.")
        self.parser.add_option("-S", "--sampler", "--sampling",
                          type="string",
                          default = "Barcoder",
                          help="Name of sampling plugin to use. Default is Barcoder. The alternative is ConstrainedNJ")
        self.parser.add_option("-A", "--alignment",
                          type="string",
                          default = "Clustalw2",
                          help="Name of alignment plugin to use. Default is ClustalW2.")
        self.parser.add_option("--inputformat",
                          type="string",
                          default="fasta",
                          help="Use this option if input file are in Nexus format. Sequences will be read from the Nexus alignment.")
        self.parser.add_option("-v", "--viewresults",
                          type="string",
                          default=False,
                          help="View results of in specified project folder in default browser.")
        self.parser.add_option("--installdependencies",
                          action="store_true",
                          default=False,
                          help="Used for installing dependencies if you want to do that before you run the program for the first time.")

        # Parallel options:
        self.parser.add_option("--hostfile",
                          type="string",
                          default=False,
                          help="Use this option to specify a file with a newline seperated list of host nanmes. The program will then use this pool of machines to distribute the assignment of individual seuences between the machines using ssh.")
        self.parser.add_option("--sge",
                          type="int",
                          default=0,
                          help="Use this option to run parallel assignments on a Sun Grid Engine cluster with a shared file system. Use the qrun.sh script to do this. The program uses the SGE queue to mannage its sub-tasks. This means that sub-tasks are submitted to the queue individually. The sge option is given two arguments. The argument is the maximal number of queued sub-jobs allowed.")

        # Blast options:
        self.parser.add_option("-t", "--minsignificance",
                          type="string",
                          default='10',
                          help="Lower sinificance bound (E-value if using GenBank) for accepting homologue sequences.")
        self.parser.add_option("-s", "--significance",
                          type="float",
                          default=0.1,
                          help="Value above which hits are considered significant")
        self.parser.add_option("-n", "--nrsignificant",
                          type="int",
                          default=5,
                          help="Minimum number of accepted significant homologues required.")
#         self.parser.add_option("-e", "--evaluecutoff",
#                           type="string",
#                           default='10',
#                           help="Lower E-value bound for accepting homologue sequences.")
#         self.parser.add_option("-t", "--evaluesignificance",
#                           type="float",
#                           default=0.1,
#                           help="E-value cutoff for significant hits")
#         self.parser.add_option("-s", "--minsignificant",
#                           type="int",
#                           default=5,
#                           help="Minimum number of accepted significant homologues required.")
        self.parser.add_option("-l", "--limitquery",
                          type="string",
                          default="all[FILT]",
                          help="Entrez query to limit blast database. Eg. 'COI[GENE] AND Hominidae[ORGANISM]' to compare only to Human COI sequences, or 'keyword[barcode]' to compare to only sequences labeled barcode entries.")
        self.parser.add_option("--minidentity",
                          type="float",
                          default=0.95,
                          help="Minimum global alignment similarity of best blast hit.")
        self.parser.add_option("-r", "--relbitscore",
                          type="float",
                          default=0.5,
                          help="Ratio of best to worst length normalized bit score accepted")
        self.parser.add_option("-m", "--maxblasthits",
                          type="int",
                          default=200,
                          help="Number of blast hits to retrieve.")
        self.parser.add_option("--nolowcomplexfilter",
                          action="store_true",
                          default=False,
                          help="Disable low-complexity filtering of query sequences.")

        # Homologue set compilation options:
        self.parser.add_option("-q", "--quickcompile",
                          action="store_true",
                          default=False,
                          help="Do not make pairwise alignments to sample-sequence to further ensure that the homologue maches up.")
        self.parser.add_option("-b", "--besthits", "--softalignmentlimit",
                          type="int",
                          default=30,
                          help="Maximal numbber of homologues to add to  alignment purely based on significant E-value")
        self.parser.add_option("-a", "--alignmentlimit",
                          dest="alignmentlimit",
                          type="int",
                          default=50,
                          help="Maximal number of sequences allowed in the homologue set.")    
        self.parser.add_option("-i", "--individuals",
                          type="int",
                          default=1,
                          help="Number of best matching individuals from a species to include in homologue set.")
        self.parser.add_option("--subspecieslevel",
                          action="store_true",
                          default=False,
                          help="Include one of each relevant subspecies (and not only species) if possible.")
        self.parser.add_option("-p", "--phyla",
                          type="int",
                          default=2,
                          help="Number of phyla to try to get into homology set.")
        self.parser.add_option("-c", "--classes",
                          type="int",
                          default=3,
                          help="")
        self.parser.add_option("-o", "--orders",
                          type="int",
                          default=5,
                          help="Number of orders to try to get into homology set.")
        self.parser.add_option("-f", "--families",
                          type="int",
                          default=6,
                          help="Number of families to try to get into homology set.")
        self.parser.add_option("-g", "--genera",
                          type="int",
                          default=10,
                          help="Number of genera to try to get into homology set.")
        self.parser.add_option("--minimaltaxonomy",
                          type="int",
                          default=5,
                          help="Minimal number of annotated taxonomic levels in accepted taxonomy annotation.")    
        self.parser.add_option("--harddiversity",
                          action="store_true",
                          default=False,
                          help="Discard homologue set if diversity specifications are not met.")    
        self.parser.add_option("--forceexcludegilist",
                          type="string",
                          #default=False,
                          default=None,
                          help="Space seperated list of gi numbers to force exclude in each homology set. Takes priority over forceincludegilist")
        self.parser.add_option("--forceincludegilist",
                          type="string",
                          #default=False,
                          default=None,
                          help="Space seperated list of gi numbers to force include in each homology set.")
        self.parser.add_option("--forceidentity",
                          type="float",
                          default=0.9,
                          help="Minimal accepted match identity of forced included GIs to query sequence.")
        self.parser.add_option("--forceincludefile",
                          type="string",
                          #default=False,
                          default=None,
                          help="File name of entries to force include in each homology set.")
        self.parser.add_option("--nofillin",
                          action="store_true",
                          default=False,
                          help="Do not Fill in more individuals up to alignment limit even if nr. of individuals can be kept even.")
        self.parser.add_option("--fillinall",
                          action="store_true",
                          default=False,
                          help="Fill in more individuals up to alignment limit disrespecting nr. of individuals.")
        self.parser.add_option("--fillintomatch",
                          type="string",
                          default=False,
                          help="Used with the force include options to fill in individuals up to the number of of homologues for the specied species or subspecies.")
        self.parser.add_option("--flanks",
                          type="int",
                          default=0,
                          help="Extra flanks to add to the homologue before aligning the set.")
        self.parser.add_option("--notruncate",
                          action="store_true",
                          default=False,
                          help="Do not truncate homologue sequence to match of sample-sequence. Used when the sample-sequence and database are prealigned.")

        # Alignment options:
        self.parser.add_option("--alignmentoption",
                          type="string",
                          action="append",
                          default=['-gapopen=50'],
                          help="Options passed to ClustalW2. Secify once for every option. Defaults to '-gapopen=50'")

        # Tree statistics options:
        self.parser.add_option("--prunelevel",
                          dest="prunelevel",
                          type="float",
                          default=0,
                          help="Prune level... float")

        # Output options:
        self.parser.add_option("-x", "--ppcutoff",
                          action="append",
                          type="int",
                          default=[95],
                          help="""Posterior probability cut-off for assignments reported on the main summary page. Default is 95% Give the options more than once to get reports for more than one cut-off""")
        self.parser.add_option("--diffs",
                          action="store_true",
                          default=False,
                          help="Calculate and show pairwise differences statistics for query sequences.")
        self.parser.add_option("--cleanlook",
                          action="store_true",
                          default=False,
                          help="Do not produce red text in summary and clone pages.")
        self.parser.add_option("--nocopycache",
                          action="store_true",
                          default=False,
                          help="Do not copy cache for duplicate query sequences.")
        self.parser.add_option("--warnongaps",
                          action="store_true",
                          default=False,
                          help="Issue a warning if the querySequence has internal gaps in the alignment.")
        self.parser.add_option("--bayesfactors",
                          action="store_true",
                          default=False,
                          help="Print bayes factors as well.")
        self.parser.add_option("--svg",
                          action="store_true",
                          default=False,
                          help="Generate publication grade SVG of the taxonomy summaries instead of ASCI trees. The SVG pictures are easily converted to EPS, PS or PSF format using the free program Incscape.")


        self.parser.add_option("--blastcache",
                          type="string",
                          #default=None,
                          default='blastcache',
                          help="Where to put cached blast files. FOR INTERNAL USE ONLY.")
        self.parser.add_option("--dbcache",
                          type="string",
                          #default=None,
                          default='dbcache',
                          help="Where to put cached GenBank files.")
        self.parser.add_option("--homologcache",
                          type="string",
                          #default=None,
                           default='homologcache',
                          help="Where to put cached homologue files. FOR INTERNAL USE ONLY.")
        self.parser.add_option("--alignmentcache",
                          type="string",
                          #default=None,
                          default='alignmentcache',
                          help="Where to put cached alignment files. FOR INTERNAL USE ONLY.")
        self.parser.add_option("--treestatscache",
                          type="string",
                          #default=None,
                          default='treestatscache',
                          help="Where to put cached tree files. FOR INTERNAL USE ONLY.")
        self.parser.add_option("--treescache",
                          type="string",
                          #default=None,
                          default='treescache',
                          help="Where to put cached tree statistics files. FOR INTERNAL USE ONLY.")
        self.parser.add_option("--statsdir",
                          type="string",
                          #default=None,
                          default='stats',
                          help="Where to put cached sequence statistics files. FOR INTERNAL USE ONLY.")
        self.parser.add_option("--resultdir",
                          type="string",
                          #default=None,
                          default='html',
                          help="Where to put cached result HTML files.")
        self.parser.add_option("--datadir",
                          type="string",
                          #default=None,
                          default='data',
                          help="Where to put the fixed input files.")

        self.parser.add_option("--_align",
                          action="store_true",
                          default=False,
                          help="Runs alignment. For internal use ony.")
        self.parser.add_option("--_sample",
                          action="store_true",
                          default=False,
                          help="Runs sampling. For internal use ony.")
        self.parser.add_option("--_stats",
                          action="store_true",
                          default=False,
                          help="Runs tree statistics. For internal use ony.")
        self.parser.add_option("--_clustalw",
                          action="store_true",
                          default=False,
                          help="Runs clustalw. For internal use ony.")

        self.parser.add_option("--mb",
                          action="store_true",
                          default=False,
                          help="Runs mrbayes.")


        self.parser.add_option("--TESTING",
                          action="store_true",
                          default=False)

        self.options, self.args = self.parser.parse_args()

    def showMessageAndExit(self, msg, guiParent=None):

        if guiParent is None:
            print msg
            sys.exit()            
        else:
            import wx
            dlg = wx.MessageDialog(guiParent, msg, 'Invalid option specification', wx.OK | wx.ICON_INFORMATION)
            dlg.ShowModal()
            dlg.Destroy()
            sys.exit()

    def postProcess(self, guiParent=None):
        """
        Post process options
        """
        
        if os.name == 'nt' and (self.options.sge or self.options.hostfile):
            self.showMessageAndExit("Parallel computing is not available on Windows.", guiParent=guiParent)

        if re.search(r' ', self.options.project):
            self.showMessageAndExit("Name of project directory can not contain space characters.", guiParent=guiParent)
        self.options.project = os.path.abspath(self.options.project)

        self.options.blastcache = os.path.join(self.options.project, self.options.blastcache)
        self.options.dbcache = os.path.join(self.options.project, self.options.dbcache)
        self.options.homologcache = os.path.join(self.options.project, self.options.homologcache)
        self.options.alignmentcache = os.path.join(self.options.project, self.options.alignmentcache)
        self.options.treestatscache = os.path.join(self.options.project, self.options.treestatscache)
        self.options.treescache = os.path.join(self.options.project, self.options.treescache)
        self.options.resultdir = os.path.join(self.options.project, self.options.resultdir)
        self.options.statsdir = os.path.join(self.options.project, self.options.statsdir)
        self.options.datadir = os.path.join(self.options.project, self.options.datadir)

        if self.options.forceincludegilist:
            self.options.forceincludegilist = re.split(r'\s+', self.options.forceincludegilist)
        else:
            self.options.forceincludegilist = []

        if self.options.forceexcludegilist:
            self.options.forceexcludegilist = re.split(r'\s+', self.options.forceexcludegilist)
        else:
            self.options.forceexcludegilist = []

        if self.options.forceexcludegilist and self.options.forceincludegilist:
            self.options.forceincludegilist = [x for x in self.options.forceincludegilist if x not in self.options.forceexcludegilist]

        if self.options.nofillin and self.options.fillinall:
            self.showMessageAndExit("Don't use the options fillinall and nofillin at the same time.", guiParent=guiParent)

        # Make sure alignment options are unique and nonoverlapping:
        alignmentoption = []
        nrGapOpenOpt = 0
        for opt in self.options.alignmentoption:
            if re.search('-gapopen', opt):
                nrGapOpenOpt += 1
        for opt in self.options.alignmentoption:
            if nrGapOpenOpt > 1 and opt == '-gapopen=50':
                nrGapOpenOpt -= 1
                continue
            else:
                alignmentoption.append(opt)
        self.options.alignmentoption = alignmentoption


        # Make sure ppcutoffs are unique and sorted:
        u = {}
        for x in self.options.ppcutoff:
            u[x] = 1
        self.options.ppcutoff = u.keys()
        self.options.ppcutoff.sort()

        return (self.options, self.args)


