
# Neighbour or trees, built by the neighbour joining algorithm

import os
from SAP.Bio.Nexus import Nexus
import random, copy, sys, re, math

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass
    
class TooFewAlnBasesError(Error):
    """
    Exception raised when there are only gaps between two gaps in the
    resampled alignment.
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

class Node:
    """
    Container class for a cluster of nodes.
    """
    def __init__(self, distList, name=None, left=None, distLeft=None, right=None, distRight=None, leafSet=None):
        self.name = name
        self.distList = distList # Distances to lower-numbered nodes, or None
        self.left = left 	# Left and right children, or None
        self.right = right
        self.distLeft  = distLeft # Length of edges to the children, if any
        self.distRight = distRight
        self.leafSet = leafSet # Store list of subnodes as set in each node

    def live(self):
        """
        Returns true if distances to lower ordered nodes are defined
        which means that the node has yet to be joined.
        """
        return self.distList is not None

    def kill(self):
        """
        Remove the node from the set we are collapsing by setting the
        distance to lower-numbered nodes to None
        """
        self.distList = None
# 
#     def __str__(self):
#         sList = []
#         if self.left and self.right:
#             ## sList.append("%s:%f" % (self.left, self.distLeft))
#             ## sList.append("%s:%f" % (self.right, self.distRight))
#             sList.append(str(self.left))
#             sList.append(str(self.right))
#         else:
#             return self.name
#         return '(' + ','.join(sList) + ')'

    def __str__(self):
        sList = []
        if self.left and self.right:
            if len(self.left.leafSet) > len(self.right.leafSet):
                return "(%s,%s)" % (str(self.left), str(self.right))
            else:
                return "(%s,%s)" % (str(self.right), str(self.left))
        else:
            return str(self.name)

class NeighbourJoin:
    """
    The neighbour-joining algorithm.  Make a rooted tree by arbitrarily
    adding a root node with edges to the last two leaves.
    """    
    def __init__(self, alignment, queryName=None, constraintTree=None):
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
        self.minPairAlnBases = 20

        self.keyList = alignment.matrix.keys()
        self.keyList.sort()

        # Make a more efficient representation:
        self.alignment = []
        self.translation = {}
        for i, key in enumerate(self.keyList):
            self.translation[key] = str(i+1)
            self.alignment.append({'name':str(i+1), 'seq':str(alignment.matrix[key])})
        self.alignmentLength = len(self.alignment[0]['seq'])
        self.indexList = range(self.alignmentLength)

        self.constraintTree = constraintTree
                    
        if constraintTree is not None and queryName is not None:
            # Find bounds of query seq and make corresponding index lists:
            self.queryName = queryName
            querySeq = str(alignment.matrix[self.queryName])
            leftBound = len(re.search("^(-*)", querySeq).groups()[0])
            rightBound = len(querySeq) - len(re.search("(-*)$", querySeq).groups()[0])
            self.queryIndexList = range(leftBound, rightBound)
            self.nonQueryIndexList = range(leftBound) + range(rightBound, self.alignmentLength)

            self.queryName = self.translation[self.queryName]
            
            # Have a list self.partitionList of partitions as sets sorted by size
            rootNode=constraintTree.root
            self.constraintList = self.constraintPartitions(constraintTree, rootNode)
            self.constraintList = filter(lambda x: len(x) > 1, self.constraintList)
            self.constraintList.sort(lambda a, b: cmp(len(a),len(b)))
            self.origConstraintList = copy.deepcopy(self.constraintList)  # extra copy for safe keeping
        else:
            self.queryIndexList = range(self.alignmentLength)
            self.nonQueryIndexList = []


            
    def constraintPartitions(self, tree, nodeID, partList=[]):
        """
        Creates a list of all partitions in the constraint tree.
        """
        node = tree.chain[nodeID]
        if node.succ != []:
            otus = tree.get_taxa(nodeID)
            otus = map(lambda x: self.translation[x], otus)
            partList.append(set(otus))
            for childID in node.succ:
                self.constraintPartitions(tree, childID, partList=partList)
        return partList

    def __str__(self):
        if not self.nodeList:
            raise Exception
        root = self.nodeList[-1]
        return str(root)

    def calculateDistanceVector(self, i):
        distances = []
        for j in range(len(self.alignment)):
            distances.append(self.distance(i, j))
        return distances

    def distance(self, i, j):
        """
        Calculates distance using Kimuras two-parameter model.
        """
        extraGaps = 0
        indels = 0
        distance = 0.0
        transitions = 0
        transversions = 0
        uninformative = 0

        # A = adenine (purine), C = cytosine (pyrimidine), G = guanine (purine), T = thymine (pyrimidine), U = uracil
        # R = G A (purine), Y = T C (pyrimidine), K = G T (keto), M = A C (amino)
        # S = G C (strong bonds), W = A T (weak bonds)
        # B = G T C (all but A), D = G A T (all but C), H = A C T (all but G), V = G C A (all but T),
        # N = A G C T (any)

        ## score = 0
        ## for idx in self.indexList:
        ##     x = self.alignment[i]['seq'][idx]
        ##     y = self.alignment[j]['seq'][idx]

        ##     if x == "-" and y == "-":
        ##         extraGaps += 1
        ##         continue
        ##     elif x == "-" or y == "-":
        ##         indels += 1
        ##         continue
        ##     elif x == y:
        ##         score += 1

        ## ident = float(score)/float(self.alignmentLength - extraGaps - indels)
        ## print ident
        ## return 1 - ident

        for idx in self.indexList:
            x = self.alignment[i]['seq'][idx]
            y = self.alignment[j]['seq'][idx]

            if x == "-" and y == "-":
                extraGaps += 1
                continue
            elif x == "-" or y == "-":
                indels += 1
                continue
            elif x in 'MKWSBDHVN' or y in 'MKWSBDHVN':
                uninformative += 1
                continue
            elif x == y:
                continue

            if x == 'A':
                if y == 'G' or y == 'R':
                    transitions += 1
                elif y == 'T' or y == 'C' or y == 'Y':
                    transversions += 1
            elif x == 'C' or y == 'Y':
                if y == 'T':
                    transitions += 1
                elif y == 'A' or y == 'G' or y == 'R':
                    transversions += 1
            elif x == 'T':
                if y == 'C' or y == 'Y':
                    transitions += 1
                elif y == 'A' or y == 'G' or y == 'R':
                    transversions += 1
            elif x == 'G':
                if y == 'A' or y == 'R':
                    transitions += 1
                elif y == 'T' or y == 'C' or y == 'Y':
                    transversions += 1
            elif x == 'R':
                if y == 'Y':
                    transversions += 1
            elif x == 'Y':
                if y == 'R':
                    transversions += 1


        if self.alignmentLength - extraGaps - indels - uninformative < self.minPairAlnBases:
            # If resampled alignment for the two seqeuences noly
            # contains gaps we throw and exception so the bootstrap
            # iteration is abborted and started over with a new
            # resampling:
            raise TooFewAlnBasesError()

        P = float(transitions) / (self.alignmentLength - extraGaps - indels - uninformative)
        Q = float(transversions) / (self.alignmentLength - extraGaps - indels - uninformative)        


        try:
            distance = 0.5 * math.log(1 / (1 - 2 * P - Q)) + 0.25 * math.log(1 / (1 - 2 * Q))
        except ValueError:
            # The distance is more than expected between randomly sampled sequences.
            # So we jus set it to a very large value:
            distance = 100
        except ZeroDivisionError:
            # The distance is more than expected between randomly sampled sequences.
            # So we jus set it to a very large value:
            distance = 100

        # Check that distance is not 'nan':
        if not (distance < 0 or distance >= 0):
            # The distance is more than expected between randomly sampled sequences.
            # So we jus set it to a very large value:
            distance = 100

        return distance

    def makeTree(self, demandAdditiveness=False):

        #self.liveNodes = []
        self.N = len(self.alignment)   # The initial number of leaves
        self.K = self.N # The number of clusters created so far
        self.r = []  # The average distance to other leaves
        self.nodeList = []
        for i in range(2*self.N-1):
            self.nodeList.append(None)
            self.r.append(None)

        self.ds = []
        for i, d in enumerate(self.alignment):
            self.ds.append(self.calculateDistanceVector(i))
            self.nodeList[i] = Node(self.ds[i], d['name'], leafSet=set([d['name']]))

        if demandAdditiveness and not self.isAdditive():
            return None

        ## ######################################
        ## width = 8
        ## print
        ## for i, key in enumerate(self.keyList):
        ##     print key[:width],
        ##     for j in range(i):
        ##         print '%-*f' % (width, self.d(i, j)),
        ##     print
        ## print ' ' * width,
        ## for i in range(len(self.keyList)-1):
        ##     print self.keyList[i][:width],
        ## print
        ## ######################################

        while self.K < 2 * self.N - 2:
            self.findAndJoin()
        # Two leaves remain; self.nodeList[self.K-1] is one of them, go find the other
        # Arbitrarily add a root node at this point
        K2 = self.K - 2
        while not self.nodeList[K2].live():
            K2 -= 1
        dij = self.d(K2, self.K-1) / 2
        self.nodeList[self.K] = Node(None, None, self.nodeList[K2], dij, self.nodeList[self.K-1], dij,
                                     self.nodeList[K2].leafSet.union(self.nodeList[self.K-1].leafSet))
        self.K += 1        

    def d(self, i, j):
        return self.nodeList[max(i, j)].distList[min(i, j)]

    def computeR(self):
        for i in range(self.K):
            if self.nodeList[i].live():
                sum = 0
                for k in range(self.K):
                    if self.nodeList[k].live()and k != i:
                        sum += self.d(i, k)
                L = 2 * self.N - self.K	    # The current number of leaves
                self.r[i] = sum / (L - 2);  # Strange, but the book says so (p 171)

    def constraintPartition(self, i, j):
        """
        Check if the taxonomic constraints allow i and j to be joined at this point.
        """
        if self.constraintTree is None:
            # Always accept:
            return True, False
        else:
            # Make sure there is no overlap in the node sets:
            assert not self.nodeList[i].leafSet.intersection(self.nodeList[j].leafSet)

            leafSet = self.nodeList[i].leafSet.union(self.nodeList[j].leafSet)
            if self.queryName in leafSet:
                leafSet.remove(self.queryName) 
            partitionIndex = None
            deletePartition = False

            for k, partition in enumerate(self.constraintList):
                if len(leafSet.intersection(partition)):
                    if leafSet.issubset(partition):
                        partitionIndex = k
                        if leafSet == partition:
                            deletePartition = True
                    break

            ## for k, partition in enumerate(self.constraintList):
            ##     if len(leafSet.intersection(partition)) and leafSet.issubset(partition):
            ##         partitionIndex = k
            ##         if leafSet == partition:
            ##             deletePartition = True
            ##     break
                
            return partitionIndex, deletePartition
    
    def findAndJoin(self):
        """
        Find closest two live clusters and join them.
        """
        self.computeR()
        mini = -1
        minj = -1
        mind = 1000000000
        pairList = []
        for i in range(self.K): 
            if self.nodeList[i].live():
                for j in range(i):
                    if self.nodeList[j].live():
                        #print i, j
                        d = self.d(i, j) - (self.r[i] + self.r[j])
                        partitionIndex, deletePartitionBool = self.constraintPartition(i, j)
                        if partitionIndex is not None:
                            if d < mind:
                                mind = d
                                pairList = [[i, j, partitionIndex, deletePartitionBool]]
                            elif d == mind:
                                pairList.append([i, j, partitionIndex, deletePartitionBool])

        assert pairList
        #print mind
        pair = random.choice(pairList)        

#         if self.K == 34:
#             pair[0] = 29
#             pair[1] = 11
# #         if self.K == 33:
# #             pair[0] = 13
# #             pair[1] = 3
#             print "#", pair[0]+1, pair[1]+1, self.d(pair[0], pair[1]), self.r[pair[0]], self.r[pair[1]], self.d(pair[0], pair[1]) - self.r[pair[0]] - self.r[pair[1]]
#             s1 = 0.0
#             s2 = 0.0
#             for f in range(self.K): 
#                 if self.nodeList[f].live():
#                     print f, self.d(pair[1], f)
#                     s1 += self.d(pair[1], f)
#                     s2 += self.d(pair[0], f)
#             print "##", s1 / ((2 * self.N - self.K) - 2), s2 / ((2 * self.N - self.K) - 2)
#             print "==", self.r[pair[0]], self.r[pair[1]]
#             sys.exit()

#         print "#", pair[0]+1, pair[1]+1, mind, self.d(pair[0], pair[1]), self.r[pair[0]], self.r[pair[1]], self.d(pair[0], pair[1]) - self.r[pair[0]] - self.r[pair[1]]

        self.join(pair[0], pair[1])
        if self.constraintTree and pair[3]:
            del self.constraintList[pair[2]]


    def join(self, i, j):
        """
        Join i and j to form node self.K
        """
        distList = []
        for x in range(self.K):
            distList.append(None)
        dij = self.d(i, j)
        for m in range(self.K):
            if self.nodeList[m].live() and m != i and m != j:
                distList[m] = (self.d(i, m) + self.d(j, m) - dij) / 2
        dik = (dij + self.r[i] - self.r[j]) / 2
        djk = dij - dik
        #print self.nodeList[i].leafSet, self.nodeList[j].leafSet,
        self.nodeList[self.K] = Node(distList, None, self.nodeList[i], dik, self.nodeList[j], djk,
                                     self.nodeList[i].leafSet.union(self.nodeList[j].leafSet))
        self.nodeList[i].kill() 
        self.nodeList[j].kill()
        self.K += 1

        ## ######################################
        ## for i, n in enumerate(self.nodeList):
        ##     if n is None:
        ##         break
        ## print self.nodeList[i-1]
        ## ######################################

    def isAdditive(self):
        for i in range(self.N):
            for j in range(i+1,self.N):
                for k in range(j+1, self.N):
                    for m in range(k+1, self.N):
                        dijdkm = self.d(i, j) + self.d(k, m)
                        dikdjm = self.d(i, k) + self.d(j, m)
                        dimdjk = self.d(i, m) + self.d(j, k)
                        if not (dijdkm == dikdjm and dijdkm >= dimdjk \
                             or dijdkm == dimdjk and dijdkm >= dikdjm \
                             or dikdjm == dimdjk and dikdjm >= dijdkm):
                            #print "(i, j, k, m) = (%d, %d, %d, %d)" % (i, j, k, m)
                            return False
        return True

    def resampleIndexList(self):
        indexList = []
        queryIndexList = self.queryIndexList
        nonQueryIndexList = self.nonQueryIndexList
        for i in self.queryIndexList:
            indexList.append(random.choice(queryIndexList))
        if not self.constraintTree and self.nonQueryIndexList:
            for i in self.nonQueryIndexList:
                indexList.append(random.choice(nonQueryIndexList))
        self.indexList = indexList

    ## def resampleIndexList(self):
    ##     indexList = []
    ##     for i in range(self.alignmentLength):
    ##         indexList.append(random.randint(0,self.alignmentLength-1))
    ##     self.indexList = indexList

    def resetIndexList(self):
        self.indexList = range(self.alignmentLength)

    def resetConstraintList(self):
        self.constraintList = copy.deepcopy(self.origConstraintList)

    def bootstrapIterator(self, iterations):        
        for i in range(iterations):
            self.resampleIndexList()
            tries = 0
            while tries < self.triesOnTooFewAlnBases:
                try:
                    self.makeTree()
                    break
                except TooFewAlnBasesError:
                    tries += 1                
            else:
                raise BootstrapError("Too many gaps. Gave up on the alignment.")
            #yield self.nexusObj
            yield str(self)
        self.resetIndexList()

    def dumpBootstrapTrees(self, iterations, fileName):

        s = "#NEXUS\n\nbegin trees;\n"
        s += "   translate\n"
        for i, key in enumerate(self.keyList):
            if i == len(self.keyList)-1:
                char = ';'
            else:
                char = ','
            s += "       %d %s%s\n" % (i+1, key, char)
        for i in range(iterations):
            ## print i, 
            self.resampleIndexList()            
            tries = 0
            while tries < self.triesOnTooFewAlnBases:
                try:
                    self.makeTree()
                    break
                except TooFewAlnBasesError:
                    tries += 1                
            else:
                raise BootstrapError("Too many gaps. Gave up on the alignment.")
            ## print "##############"
            s += "   tree rep_%d = %s;\n" % (i+1, str(self))
            if self.constraintTree:
                self.resetConstraintList()
        s += "end;\n"
        fp = open(fileName, 'w')
        fp.write(s)
        fp.close()
        self.resetIndexList()

    ## def dumpBootstrapTrees(self, iterations, fileName):
    ##     fp = open(fileName, 'w')
    ##     fp.write("#NEXUS\n\nbegin trees;\n")
    ##     fp.write("   translate\n")
    ##     for i, key in enumerate(self.keyList):
    ##         if i == len(self.keyList)-1:
    ##             char = ';'
    ##         else:
    ##             char = ','
    ##         fp.write("       %d %s%s\n" % (i+1, key, char))
    ##     for i in range(iterations):
    ##         self.resampleIndexList()
    ##         self.makeTree()
    ##         fp.write("   tree rep_%d = %s;\n" % (i+1, str(self)))
    ##         if self.constraintTree:
    ##             self.resetConstraintList()
    ##     fp.write("end;\n")
    ##     fp.close()
    ##     self.resetIndexList()

    def nexus(self):
        s = "#NEXUS\n\nbegin trees;\n   tree NJ_tree = %s;\nend;\n" % str(self)
        return s

    def dumpNexus(self, fileName):
        fp = open(fileName, 'w')
        s = "#NEXUS\n\nbegin trees;\n"
        s += "   translate\n"
        for i, key in enumerate(self.keyList):
            if i == len(self.keyList)-1:
                char = ';'
            else:
                char = ','
            s += "       %d %s%s\n" % (i+1, key, char)
        s += "   tree NJ_tree = %s;\nend;\n" % str(self)
        fp.write(s)
        fp.close()

    def nexusString(self):
        s = "#NEXUS\n\nbegin trees;\n"
        s += "   translate\n"
        for i, key in enumerate(self.keyList):
            if i == len(self.keyList)-1:
                char = ';'
            else:
                char = ','
            s += "       %d %s%s\n" % (i+1, key, char)
        s += "   tree NJ_tree = %s;\nend;\n" % str(self)
        print s

if __name__ == "__main__":


    # tests:
    # alignment with all identical sequences should return equal bootstrap values.
    # check that generation of indexLists for query and rest work as I think



    # python newversion.py ~/projects/sap/analyses/benchmark/taxonophy/benchmarkAb1BC/alignments/InsectaAb1_15077724.nex ~/projects/sap/analyses/benchmark/taxonophy/benchmarkAb1BC/homologues/InsectaAb1_15077724.nex


    random.seed(7)

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
 
#    cnj = ConstrainedNJ(alignment, queryName='release_JH46t_text_so', constraintTree=constraintTree)
    cnj = NeighbourJoin(alignment, queryName='InsectaAb1_15077724', constraintTree=constraintTree)    


    cnj.makeTree()
    print "#NEXUS\n\nbegin trees;\n   tree NJ_tree = %s;\nend;" % str(cnj)
#    print cnj.nexusString()

# 
#     alignmentFileName = sys.argv.pop()
#     if not os.path.exists(alignmentFileName):
#         raise Exception
#     alignment = Nexus.Nexus(alignmentFileName)
# 
#     nj = NeighbourJoin(alignment)    
# #     nj.makeTree()
# #     nj.dumpNexus('~/Deskopt/nei.nex')
#     nj.dumpBootstrapTrees(10, '~/Deskopt/tmp.nex')


    ## from SAP.Bio.Nexus import Trees as NexusTrees
    ## nex = Nexus.Nexus('tmp.nex')
    ## cons = NexusTrees.consensus(nex.trees)
    ## print cons.to_string()

    ## print nj

    ## THE ITERAATOR SHOULD RETURN NEXUS OBJECTS. WE SHOULD HAVE
    ## ANOTHER MEMTHOD - dumpBootstrapTrees - THAT DUMPS NEXUS TO A
    ## FILE.

    ## fp = open('boot.nex', 'w')
    ## for tree in nj.bootstrapIterator(10):
    ##     pass
    ##     fp.write(tree + "\n")
    ## fp.close()

    
    ## nj.checkAdditivity()


    ## The resampling should sample proportionally from the query part of the
    ## alignmen and the rest of the alignment.

    ## Implement Kimulra 2 parameter model.
