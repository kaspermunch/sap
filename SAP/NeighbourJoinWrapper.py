
import sys, os
from SAP.Bio.Nexus import Nexus
from SAP.Bio.Nexus import Trees as NexusTrees

from UtilityFunctions import *
import NeighbourJoin

class NeighbourJoinWrapper:

    def __init__(self, options):
        self.options = options

    def run(self, alignmentFileName):
        """
        Run neighbour join bootstrap on alignment
        Returns name of tree file
        """

        baseName = os.path.splitext(os.path.split(alignmentFileName)[-1])[0]
        treeFileName = os.path.join(self.options.treescache, baseName + ".nex.nj.tree")

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

            queryName = baseName
            alignment = Nexus.Nexus(alignmentFileName)
            nj = NeighbourJoin.NeighbourJoin(alignment, queryName=queryName, constraintTree=constraintTree)
            nj.dumpBootstrapTrees(self.options.bootstraps, treeFileName)

            bootstrapTrees = Nexus.Nexus(treeFileName)
            consensusTree = NexusTrees.consensus(bootstrapTrees.trees)
            consensusTreeFileName = os.path.join(self.options.treescache, baseName + ".nj.con")

            s = "#NEXUS\n\nbegin trees;\n   " + consensusTree.to_string() + "\nend;\n"
            writeFile(consensusTreeFileName, s)

            print "done."
            sys.stdout.flush()
