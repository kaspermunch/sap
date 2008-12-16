
import sys, os, tempfile
from SAP.Bio.Nexus import Nexus
from SAP.Bio.Nexus import Trees as NexusTrees

from SAP.UtilityFunctions import *
from ConstrainedNJ import ConstrainedNJ

# The sampler writes the following files to the cache dir (existence of the nex file should indication proper completion):
# <filename>_<query>.<samplername>.nex
# <filename>_<query>.<samplername>.out
# <filename>_<query>.<samplername>.err

class Assignment:

    def __init__(self, options):
        self.options = options

        self.name = 'ConstrainedNJ'

    def run(self, alignmentFileName):
        """
        Run neighbour join bootstrap on alignment
        Returns name of tree file
        """

        baseName = os.path.splitext(os.path.split(alignmentFileName)[-1])[0]
        treeFileName = "%s/%s.%s" % (self.options.treescache, baseName, self.name + ".nex")

        print "%s: Sampling trees: " % baseName,
        sys.stdout.flush()

        # TODO: The way we check if there are cached results is not safe. The file may not be complete.
        if os.path.exists(treeFileName):
            print "Using cached results."
            sys.stdout.flush()
        else:
            print "Computing...",
            sys.stdout.flush()

            constraintTreeFileName = os.path.join(self.options.homologcache, baseName + ".nex")
            constraintTreeFile = open(constraintTreeFileName, 'r')
            constraintTree = Nexus.Nexus(constraintTreeFile).trees[0]

            bootstraps = 1000

            queryName = baseName
            alignment = Nexus.Nexus(alignmentFileName)
            cnj = ConstrainedNJ(alignment, constraintTree=constraintTree)

            cnj.dumpBootstrapTreesNexus(bootstraps, treeFileName)

            # Make a consensus tree and write it to a seperate file:
#             bootstrapTrees = Nexus.Nexus(treeFileName)
#             consensusTree = NexusTrees.consensus(bootstrapTrees.trees)
#             consensusTreeFileName = os.path.join(self.options.treescache, baseName + ".nj.con")
# 
#             s = "#NEXUS\n\nbegin trees;\n   " + consensusTree.to_string() + "\nend;\n"
#             writeFile(consensusTreeFileName, s)

            print "done."
            sys.stdout.flush()
