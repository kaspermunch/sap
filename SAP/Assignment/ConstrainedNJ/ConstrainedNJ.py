
import os, random, copy, sys, re, math
from SAP.Bio.Nexus import Nexus

try:
    # Use the fast C++ implementation:
    from cConstrainedNJlib import computeTree, initialize
except:
    # Fall back on the python implementation:
    print "WARNING: C++ version of ConstrainedNJlib not found. Using the slower python implementation."
    from ConstrainedNJlib import computeTree, initialize

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass

class BootstrapError(Error):
    """
    Exception raised when bootstrap failes.

    Attributes:
        expression -- input expression in which
                      the error occurred
        message -- explanation of the error
    """
    def __init__(self, message):
        self.message = message
        print self.message

class ConstrainedNJ(object):

    def __init__(self, alignment, constraintTree=None):
        """
        Takes an alignment object as input and s
        """        
        self.nodeList = None
        
        # How many times the class re-tries each bootstrap iteration
        # of a sequence pair in the resampled alignment only has gaps:
        self.triesOnTooFewAlnBases = 10

        # Minimum number of aligned bases in re-sampled alignment that
        # each sequence pair must have for the bootstrap iteration not
        # to be aborted and restarted.
        self.minPairAlnBases = 50

        self.keyList = alignment.matrix.keys()
        self.keyList.sort()

        # Make a more efficient representation:
        self.alignment = []
        self.translation = {}
        for i, key in enumerate(self.keyList):
            #self.translation[key] = str(i+1)
            self.translation[key] = str(i+1)
            self.alignment.append(alignment.matrix[key].tostring())
        self.alignmentLength = len(self.alignment[0])

        self.constraintTree = constraintTree

        self.constraintList = []
                    
        if constraintTree is not None:            
            # Have a list self.partitionList of partitions as sets sorted by size
            rootNode=constraintTree.root
            self.constraintList = self.constraintPartitions(constraintTree, rootNode)
            self.constraintList = filter(lambda x: len(x) > 1, self.constraintList)
            self.constraintList.sort(lambda a, b: cmp(len(a),len(b)))
            self.constraintList = [" ".join(x) for x in self.constraintList]
            self.origConstraintList = copy.deepcopy(self.constraintList)  # extra copy for safe keeping

        initialize(int(self.alignmentLength/2))

    def constraintPartitions(self, tree, nodeID, partList=[]):
        """
        Creates a list of all partitions in the constraint tree.
        """
        node = tree.chain[nodeID]
        if node.succ != []:
            otus = tree.get_taxa(nodeID)
            otus = map(lambda x: self.translation[x], otus)
            partList.append(otus)
            for childID in node.succ:
                self.constraintPartitions(tree, childID, partList=partList)
        return partList

    def resetConstraintList(self):
        self.constraintList = copy.deepcopy(self.origConstraintList)

    def getTreeString(self, resample=False):
        if resample:
            resample = 1
        else:
            resample = 0
        return computeTree(self.alignment, self.constraintList, resample)

    def getNexusString(self):

        treeString = self.getTreeString()
        s = "#NEXUS\n\nbegin trees;\n"
        s += "   translate\n"
        for i, key in enumerate(self.keyList):
            if i == len(self.keyList)-1:
                char = ';'
            else:
                char = ','
            s += "       %d %s%s\n" % (i+1, key, char)
        s += "   tree CNJ_tree = %s;\nend;\n" % treeString
        return s

    def dumpNexusString(self, fileName):
        s = self.getNexusString()
        fp = open(fileName, 'w')        
        fp.write(s)
        fp.close()

    def getBootstrapTreesNexus(self, iterations):

        # Generate the header:
        s = "#NEXUS\n\nbegin trees;\n"
        s += "   translate\n"
        for i, key in enumerate(self.keyList):
            if i == len(self.keyList)-1:
                char = ';'
            else:
                char = ','
            s += "       %d %s%s\n" % (i+1, key, char)
        # Loop over the bootstraps:
        for i in range(iterations):
            tries = 0
            while tries < self.triesOnTooFewAlnBases:
                try:
                    treeStr = self.getTreeString(resample=True)
                    break
                except RuntimeError, e:
                    if str(e) == "TooFewAlnBasesError":
                        tries += 1
                    else:
                        raise RuntimeError(str(e))
            else:
                raise BootstrapError("Not enough aligned bases for reliable bootstrapping.")
            s += "   tree rep_%d = %s;\n" % (i+1, treeStr)
            if self.constraintTree:
                self.resetConstraintList()
        s += "end;\n"
        return s

    def dumpBootstrapTreesNexus(self, iterations, fileName):

        s = self.getBootstrapTreesNexus(iterations)
        fp = open(fileName, 'w')
        fp.write(s)
        fp.close()    
        
    def consensusTree(self, iterations):

        import tempfile
        from SAP.Bio.Nexus import Trees as NexusTrees

        tmpFile, tmpFileName = tempfile.mkstemp()
        self.dumpBootstrapTreesNexus(iterations, tmpFileName)
        bootstrapTrees = Nexus.Nexus(tmpFileName)
        os.remove(tmpFileName)        
        consensusTree = NexusTrees.consensus(bootstrapTrees.trees)

        return "#NEXUS\n\nbegin trees;\n   " + consensusTree.to_string() + "\nend;\n"

    ## def bootstrapIterator(self, iterations):        
    ##     for i in range(iterations):
    ##         self.resampleIndexList()
    ##         tries = 0
    ##         while tries < self.triesOnTooFewAlnBases:
    ##             try:
    ##                 treeString = self.getTreeString()
    ##                 break
    ##             except TooFewAlnBasesError:
    ##                 tries += 1                
    ##         else:
    ##             raise BootstrapError("Too many gaps. Gave up on the alignment.")
    ##         #yield self.nexusObj
    ##         yield treeString
    ##     self.resetIndexList()

if __name__ == "__main__":

    # tests:
    # alignment with all identical sequences should return equal bootstrap values.

#     # REMEMBER TO REMOVE THIS AGAIN.....
#     random.seed(7)
#     # REMEMBER TO REMOVE THIS AGAIN.....

    alignmentFileName = sys.argv[1]
    if not os.path.exists(alignmentFileName):
        raise Exception
    alignment = Nexus.Nexus(alignmentFileName)

    constraintTree = None
    if len(sys.argv) > 2:
        constraintFileName = sys.argv[2]
        if not os.path.exists(constraintFileName):
            raise Exception
        nexus = Nexus.Nexus(constraintFileName)
        constraintTree = nexus.trees[0]
 
    cnj = ConstrainedNJ(alignment, constraintTree=constraintTree)

    print cnj.getTreeString()

#     cnj.dumpNexusString('/Users/kasper/Desktop/cnj.nex')
#     cnj.dumpBootstrapTreesNexus(1000, '/Users/kasper/Desktop/tmp.nex')
#     cnj.dumpBootstrapTreesNexus(1000, 'c.nex')
#     print cnj.consensusTree(1000)
#     print cnj.getNexusString()
#     print cnj.getBootstrapTrees(1000, '/Users/kasper/Desktop/tmp.nex')




