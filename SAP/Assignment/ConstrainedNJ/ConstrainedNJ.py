
import os, random, copy, sys, re, math
from Bio.Nexus import Nexus

#from cConstrainedNJlib import computeTree, initialize
import ctypes

# try:
#     # Use the fast C++ implementation:
#     from cConstrainedNJlib import computeTree, initialize
# except:
#     # Fall back on the python implementation:
#     print "WARNING: C++ version of ConstrainedNJlib not found. Using the slower python implementation."
#     from ConstrainedNJlib import computeTree, initialize


class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    def __init__(self, message):
        self.msg = message
    pass

class BootstrapError(Error):
    """
    Exception raised when bootstrap failes.
    """
    def __init__(self, message):
        self.msg = message

class ConstrainedNJ(object):

    def __init__(self, alignment, constraintTree=None):
        """
        Takes an Nexus alignment object as input (or a diatance matrix).
        """        

        self.cnjlib = ctypes.CDLL(os.path.join(os.path.dirname(__file__), "cConstrainedNJlib.so"))
#         self.cnjlib.initialize.restype = None
#         self.cnjlib.initialize.argTypes = [ctypes.c_int]
#         self.cnjlib.cleanup.restype = None
#         self.cnjlib.cleanup.argTypes = [ctypes.c_int]
        self.cnjlib.computeTree.restype = ctypes.c_char_p
        self.cnjlib.computeTree.argTypes = [ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]

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
            self.alignment.append(str(alignment.matrix[key]))
        self.alignmentLength = int(len(self.alignment[0]))

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

# #         initialize(int(self.alignmentLength))
#         self.cnjlib.initialize(int(self.alignmentLength))
# #         cache = ((ctypes.c_double * self.alignmentLength) * self.alignmentLength)()
# #         print cache
# #         self.cnjlib.initialize(len(self.keyList), ctypes.byref(cache))
# 
#     def __del__(self):
#         self.cnjlib.cleanup(self.alignmentLength)

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

    def getTreeString(self, resample=False, branchlengths=False):
        # return computeTree(self.alignment, self.constraintList, int(resample), int(branchlengths))
        def f(L):
            a = (ctypes.c_char_p * len(L))()
            a[:] = L
            return len(L), a
#        print [a for a in self.alignment]
        return self.cnjlib.computeTree(*f(self.alignment) + f(self.constraintList) + (int(resample), int(branchlengths)))

    def getNexusString(self, resample=False, branchlengths=False):

        treeString = self.getTreeString(resample=resample, branchlengths=branchlengths)
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

    def dumpNexusString(self, fileName, branchlengths=False):
        s = self.getNexusString(branchlengths=branchlengths)
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
        tmpFile.close()
        os.remove(tmpFileName)        
        consensusTree = NexusTrees.consensus(bootstrapTrees.trees)

        return "#NEXUS\n\nbegin trees;\n   " + consensusTree.to_string() + "\nend;\n"


    def supportedTree(self, iterations, branchlengths=True, outgroup=None):

        import tempfile
        from SAP.Bio.Nexus import Trees as NexusTrees

        tree = Nexus.Nexus(self.getNexusString(branchlengths=branchlengths)).trees[0]        
        bootstrapTrees = Nexus.Nexus(self.getBootstrapTreesNexus(iterations)).trees

        total=len(bootstrapTrees)
        if total==0:
            return None
        # shouldn't we make sure that it's NodeData or subclass??
        dataclass=bootstrapTrees[0].dataclass
        max_support=bootstrapTrees[0].max_support
        clades={}
        #countclades={}
        alltaxa=set(bootstrapTrees[0].get_taxa())
        # calculate calde frequencies
        c=0
        for t in bootstrapTrees:
            c+=1
            #if c%50==0:
            #    print c
            if alltaxa!=set(t.get_taxa()):
                raise TreeError, 'Trees for consensus must contain the same taxa'
            t.root_with_outgroup(outgroup=outgroup)
            for st_node in t._walk(t.root):
                subclade_taxa=t.get_taxa(st_node)
                subclade_taxa.sort()
                subclade_taxa=str(subclade_taxa) # lists are not hashable
                if subclade_taxa in clades:
                    clades[subclade_taxa]+=float(t.weight)/total
                else:
                    clades[subclade_taxa]=float(t.weight)/total

        if alltaxa!=set(tree.get_taxa()):
            raise TreeError, 'Tree must contain the same taxa as trees for consensus'
        tree.root_with_outgroup(outgroup=outgroup)
        for st_node in tree._walk(tree.root):
            subclade_taxa=tree.get_taxa(st_node)
            subclade_taxa.sort()
            subclade_taxa=str(subclade_taxa) # lists are not hashable
            n = tree.node(st_node)
            if subclade_taxa in clades:
                n.data.support = clades[subclade_taxa]
            else:
                n.data.support = 0.0

        return "#NEXUS\n\nbegin trees;\n   " + tree.to_string(plain=False) + "\nend;\n"
        

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

#    print cnj.getNexusString(branchlengths=True)
        
#     for i in range(1):
#         print cnj.getTreeString()
        
    #print cnj.getBootstrapTreesNexus(1000)



#     cnj.dumpNexusString('/Users/kasper/Desktop/cnj.nex')
#     cnj.dumpBootstrapTreesNexus(100, '/Users/kasper/Desktop/tmp.nex')
#     cnj.dumpBootstrapTreesNexus(1000, 'c.nex')
#     print cnj.consensusTree(10)
    print cnj.supportedTree(100, outgroup="ANOCA")
#     print cnj.getNexusString()
#     print cnj.getBootstrapTrees(1000, '/Users/kasper/Desktop/tmp.nex')




