
import os
from SAP.Bio.Nexus import Nexus
import random, copy, sys, re, math
from decimal import Decimal

# Array of array alignment[<sequence_nr>][base_nr]
alignment = None          

# Array of arrays each specifying a constraint.  (name integers)
backboneSetsList = None

# Array of tuples specifying allowed pairs.     (index integers)
allowedPairsList = None

# Array of otus not part of back bone constraint.
unconstrainedOTUs = None

N = None
K = None
L = None

# Sequence alignment length:
alignmentLength = None

# list of length alignmentLength for indexing into the alignment
alignmentIndexList = None

# Minimum number of aligned bases allowed. 
minPairAlnBases = 20

# List of nodes specifying the growing tree:
nodeList = []

# index list to mapping 0 -> L integers to indexes in the matrices.
idxL = []

# index list to mapping 0 -> L integers to indexes in the node list.
idxN = []

# matrix of d values
Qmatrix = []

# matrix of distances  - replaces getDist
distMatrix =[]

# Vector for holding a list of r values:
rSums = []

# Boolean indicating if resamplingMatrix has been built
resamplingMatrixComputed = False

# A matrix of [transition, transversion, length] lists we then can sample distances from
resamplingMatrix = []

# For cached distances:
k2pCacheMatrix = []

# class Error(Exception):
#     """
#     Base class for exceptions in this module.
#     """
#     pass
#     
# class TooFewAlnBasesError(Error):
#     """
#     Exception raised when there are only gaps between two gaps in the
#     resampled alignment.
#     """
#     pass


class Node:
    """
    Container class for a cluster of nodes.
    """
    def __init__(self, name=None, left=None, distLeft=None, right=None, distRight=None, leafSet=None):
        self.name = name
        self.left = left 	# Left and right children, or None
        self.right = right
        self.distLeft  = distLeft # Length of edges to the children, if any
        self.distRight = distRight
        self.leafSet = leafSet # Store list of subnodes as set in each node

#     def __str__(self):
#         if self.left and self.right:
#             if len(self.left.leafSet) > len(self.right.leafSet):
#                 return "(%s,%s)" % (str(self.left), str(self.right))
#             else:
#                 return "(%s,%s)" % (str(self.right), str(self.left))
#         else:
#             return str(self.name)

    def __str__(self):
        if self.left and self.right:
            if len(self.left.leafSet) > len(self.right.leafSet):
                return "(%s,%s):%s" % (str(self.left), str(self.right), str(int(self.name)-1)) # <- substracting one is a hack
            else:
                return "(%s,%s):%s" % (str(self.right), str(self.left), str(int(self.name)-1)) # <- substracting one is a hack
        else:
            return str(int(self.name)-1) # <- substracting one is a hack

#     def __str__(self):
#         if self.left and self.right:
#             if len(self.left.leafSet) > len(self.right.leafSet):
#                 return "(%s:%.3f,%s:%.3f)" % (str(self.left), self.distLeft, str(self.right), self.distRight)
#             else:
#                 return "(%s:%.3f,%s:%.3f)" % (str(self.right), self.distRight, str(self.left), self.distLeft)
#         else:
#             return str(self.name)

def xuniqueCombinations(items, n):
    """
    Retruns the unique combinations of the items in the list.
    """
    if n==0: yield []
    else:
        for i in xrange(len(items)):
            for cc in xuniqueCombinations(items[i+1:],n-1):
                yield [items[i]]+cc

def printDistMatrix():
    """
    Prints the distance matrix.
    """
    print
    for i in range(N):
        print "%-3d" % i,
        #for j in range(i):
        for j in range(N):
            print '%.3f' % distMatrix[i][j],
        print " %d" % i
    print ' ' * 3,
    for i in range(N):
        print "%-5d" % i,
    print

def printQmatrix():
    """
    Prints the distance matrix.
    """
    print
    for i in range(N):
        print "%-3d" % i,
        #for j in range(i):
        for j in range(N):
            print '%.3f' % Qmatrix[i][j],
        print " %d" % i
    print ' ' * 4,
    for i in range(N):
        print "%-6d" % i,
    print

## Calculating K2P distances or getting them from a cache: #######################################

def k2pStats(i, j):
    """
    Calculates distance using Kimuras two-parameter model.
    """

    global alignment, alignmentLength, alignmentIndexList
    
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

    #for idx in range(alignmentLength):
    for idx in alignmentIndexList:
        x = alignment[i][idx]
        y = alignment[j][idx]

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
        elif x == 'C':
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
            if y == 'T' or y == 'C' or y == 'Y':
                transversions += 1
        elif x == 'Y':
            if y == 'A' or y == 'G' or y == 'R':
                transversions += 1

    alignedBases = alignmentLength - extraGaps - indels - uninformative

    if alignedBases < minPairAlnBases:
        # If resampled alignment for the two seqeuences noly
        # contains gaps we throw and exception so the bootstrap
        # iteration is abborted and started over with a new
        # resampling:
        #raise TooFewAlnBasesError()
        raise RuntimeError("TooFewAlnBasesError")

    return (transitions, transversions, alignedBases)


def k2pDist(transitions, transversions, alignedBases):

    global k2pCacheMatrix

    if alignedBases == alignmentLength and k2pCacheMatrix[transitions][transversions] != -1:
        distance = k2pCacheMatrix[transitions][transversions]
    else:
        P = float(transitions) / alignedBases
        Q = float(transversions) / alignedBases

        try:
            distance = 0.5 * math.log(1 / (1 - 2 * P - Q)) + 0.25 * math.log(1 / (1 - 2 * Q))
        except ValueError, ZeroDivisionError:
            # The distance is more than expected between randomly sampled sequences.
            # So we jus set it to a very large value:
            distance = Decimal("Infinity")
            #distance = 100

        if (distance == distance) == False: # This is a NaN (not a number)
            # The distance is more than expected between randomly sampled sequences.
            # So we jus set it to a very large value:
            distance = Decimal("Infinity")
            #distance = 100

        k2pCacheMatrix[transitions][transversions] = distance

    return distance

def computeDistance(i, j):

    transitions, transversions, alignedBases = k2pStats(i, j)

    return k2pDist(transitions, transversions, alignedBases)

def resampleDistance(i, j):

    global resamplingMatrix

    transitions, transversions, alignedBases = resamplingMatrix[i][j]

    transitions = int(round(random.randint(1, alignedBases) * (float(transitions) / alignedBases)))
    transversions = int(round(random.randint(1, alignedBases) * (float(transversions) / alignedBases)))

    return k2pDist(transitions, transversions, alignedBases)

## Maintaining the Q matrix: ######################################################################

def calcRsum(i):

    global L, distMatrix
  
    r = 0.0
    for j in range(L):
        r += distMatrix[i][idxL[j]]
    return r / (L - 2)

def updateMatrices(old_i, old_j):

    global L, idxL, Qmatrix, distMatrix, rSums

    # Remove the node with the largest integer name from the list of nodes to be joined:
    for k in range(L):
        if idxL[k] == old_j:
            del idxL[k]
            break
    # Decrement the number of nodes to be joined correspondingly:
    L -= 1

    # Update the column old_i in the distance matrix that now represents the joined pair:
    for k in range(L):
        distMatrix[idxL[k]][old_i] = ( -distMatrix[old_i][old_j] + distMatrix[idxL[k]][old_j] + distMatrix[idxL[k]][old_i] ) / 2.0
    # Updata the corresponding row:
    for k in range(L):
        distMatrix[old_i][idxL[k]] = distMatrix[idxL[k]][old_i]

    # If there are only two nodes left we don't need to update the Qmatrix:
    if L == 2:
        return

    # Calculate r sums:
    for k in range(N):
        rSums[k] = 0
    # Recalculate relevant entries in Qmatrix:
    for k, l in allowedPairsList:
        if not rSums[k]:
            rSums[k] = calcRsum(k)
        if not rSums[l]:
            rSums[l] = calcRsum(l)
        Qmatrix[k][l] = distMatrix[k][l] - (rSums[k] + rSums[l])
        Qmatrix[l][k] = Qmatrix[k][l]

#     # Calculate r sums:
#     for k in range(L):
#         rSums[k] = calcRsum(idxL[k]);
# 
#     # Recalculate Qmatrix:
#     for k in range(L):
#         for l in range(L):
#             Qmatrix[idxL[k]][idxL[l]] = distMatrix[idxL[k]][idxL[l]] - (rSums[k] + rSums[l])
#             Qmatrix[idxL[l]][idxL[k]] = Qmatrix[idxL[k]][idxL[l]]

## Calculating and maintaining the list of pairs allowed to join: ################################

def updateConstraints(i, j):

    global backboneSetsList, K
        
    removedConstraint = False
    
    for k in xrange(len(backboneSetsList)-1, -1, -1):
        if j+1 in backboneSetsList[k]:
            backboneSetsList[k].remove(j+1)
            if len(backboneSetsList[k]) == 1:
                removedConstraint = True
                del backboneSetsList[k]

    return removedConstraint

def computeAllowedPairs():

    global N, allowedPairsList, backboneSetsList

    uDict = {}

    # Initialize allowedPairsList
    allowedPairsList = []
    for i in range(N):
        inConstraint = False
        constrainedSet = set()
        for s in backboneSetsList:
            if i+1 in s:
                inConstraint = True
                nodeIndeces = [x-1 for x in s if x not in constrainedSet] # we subctract one to convert to node indeces
                for p in xuniqueCombinations(nodeIndeces, 2):
                    p.sort()
                    uDict[str(p)] = p
                    ## if i in p and p not in allowedPairsList: # THIS IS KIND OF EXPENSIVE.....
                    ##     allowedPairsList.append(p)
                break
            else:
                constrainedSet = constrainedSet.union(s)

        if not inConstraint:
            for j in range(N):
                if i != j:
                    p = [i, j]
                    p.sort()
                    uDict[str(p)] = p
                    ## if i in p and p not in allowedPairsList: # THIS IS KIND OF EXPENSIVE.....
                    ##     allowedPairsList.append(t)

    allowedPairsList = uDict.values()

def updateAllowedPairs(i, j):

    global N, allowedPairsList, backboneSetsList, unconstrainedOTUs

#     if i in unconstrainedOTUs:
#         # If i is not in the backbone sets we don't want i to represent the child
#         # nodes. If they are both unconstrained swopping is not needed but does not matter.
#         i, j = j, i 

    removedConstraint = updateConstraints(i, j)

    uDict = {}

    for k in xrange(len(allowedPairsList)-1, -1, -1):
        
        if allowedPairsList[k][0] == j or allowedPairsList[k][1] == j:

            if j in unconstrainedOTUs:
                continue

            if allowedPairsList[k] == [i, j]  or allowedPairsList[k] == [j, i]: # we may have swopped them.
                continue
            elif allowedPairsList[k][0] == j:
                p = [i, allowedPairsList[k][1]]
                p.sort()
                uDict[str(p)] = p
            elif allowedPairsList[k][1] == j:
                p = [allowedPairsList[k][0], i]
                p.sort()
                uDict[str(p)] = p
            else:
                raise Exception
        else:
            uDict[str(allowedPairsList[k])] = allowedPairsList[k]
                        
    if removedConstraint:
        # if we removed a set it means that i now represents all otus
        # from that set. So we find the first set in the updated
        # constraints list that includes i and then add all
        # combinations of i to items in this set:
        constrainedSet = set()
        for s in backboneSetsList:
            if i+1 in s:
                nodeIndeces = [x-1 for x in s if x not in constrainedSet] # we subctract one to convert to node indeces
                for p in xuniqueCombinations(nodeIndeces, 2):
                    p.sort()
                    uDict[str(p)] = p
                break
            else:
                constrainedSet = constrainedSet.union(s)

    allowedPairsList = uDict.values()

# def check(i, j):
#     """
#     Check if the taxonomic constraints allow i and j to be joined at this point.
#     """
#     # Make sure there is no overlap in the node sets:
#     assert not nodeList[i].leafSet.intersection(nodeList[j].leafSet)
# 
#     leafSet = nodeList[i].leafSet.union(nodeList[j].leafSet)
#     partitionIndex = None
#     deletePartition = False
# 
#     for k, partition in enumerate(backboneSetsList):
#         if len(leafSet.intersection(partition)):
#             if leafSet.issubset(partition):
#                 partitionIndex = k
#                 if leafSet == partition:
#                     deletePartition = True
#             break
# 
#     return partitionIndex

## Finding a pair to join and create a new node: ##################################################

def findPair():
    """
    Find closest two live clusters and join them.
    """
    global allowedPairsList, backboneSetsList, K, N, r, nodeList

    mini = -1
    minj = -1
    mind = 1000000000
    pairList = []

    for i, j in allowedPairsList:
        if Qmatrix[i][j] < mind:
            mind = Qmatrix[i][j]
            pairList = [[i, j]]
            pairList[0].sort()
        elif Qmatrix[i][j] == mind:
            pairList.append([i, j])
            pairList[-1].sort()
    assert pairList, allowedPairsList

    return random.choice(pairList)

def createParentNode(i, j):
    """
    Join i and j to form node K
    """
    global K, nodeList, idxN, distMatrix

    dik = (distMatrix[i][j] + calcRsum(i) - calcRsum(j)) / 2
    djk = distMatrix[i][j] - dik

    #print idxN, i, j, idxN[i], idxN[j], nodeList[idxN[i]].leafSet, nodeList[idxN[j]].leafSet,nodeList[idxN[i]].leafSet.union(nodeList[idxN[j]].leafSet)
    #print i, j, idxN

    nodeList[K] = Node(K+1, nodeList[idxN[i]], dik, nodeList[idxN[j]], djk,
                       nodeList[idxN[i]].leafSet.union(nodeList[idxN[j]].leafSet))
    idxN[i] = K

## The main functions: ############################################################################
     
def initialize(dim):
    """
    Initialize matrix for cached distances.
    """
    global k2pCacheMatrix

    for i in range(dim):
        k2pCacheMatrix.append([])
        for j in range(dim):
            k2pCacheMatrix[i].append(-1)

def computeTree(_alignment, _backboneSetsList, resample):
    """
    Main function to compute a neighbour-joining tree under a back bone constraint
    specified by a set of constraints. Constraints are specified as a list of lists of
    nodes allowed to join. Whether a pair is allowed to join is established by finding the
    first list that any of the two nodes is part of. If the other node is not part of that
    list too the join is not allowed.
    """

#     # REMEMBER TO REMOVE THIS AGAIN.....
#     random.seed(7)
#     # REMEMBER TO REMOVE THIS AGAIN.....

    global alignment, backboneSetsList, N, K, L, idxL, idxN, nodeList, \
           alignmentLength, resamplingMatrixComputed, resamplingMatrix, unconstrainedOTUs, alignmentIndexList

    # Populate the global variables:
    alignment =  _alignment

    # Sequence length of the alignment:
    alignmentLength = len(alignment[0])

    # The initial number of leaves
    N = len(alignment)
    # The number of clusters created so far (including the first N "leaf culsters" we start with)
    K = N
    # The number of nodes left to be joined:
    L = N

    # Allocate for the list of nodes:
    nodeList = []
    for i in range(2*N-1):
        nodeList.append(None)

    # Turn lists specifying unconstrained sets in to set objects:
    backboneSetsList = []
    if len(_backboneSetsList) > 0:
        for s in _backboneSetsList:
            l = s.split(' ')
            l = [int(x) for x in l]
            backboneSetsList.append(set(l))
    
    # Initialize the vector that keeps track of indeces in the distMatrix:
    idxL = range(N)
    # Initialize the vector that keeps track of indeces in the nodeList:
    idxN = range(N)

    # Find the otus that are not part of the backbone constraint:
    unionSet = set()
    for s in backboneSetsList:
        unionSet = unionSet.union(s)
    unconstrainedOTUs = []
    for i in range(N):
        if i+1 not in unionSet:
            unconstrainedOTUs.append(i)

    # Initialize rSums:
    for i in range(N):
        rSums.append(None)

        
    # Create distMatrix - like alocating memory in C
    for i in range(N):
        distMatrix.append([])
        for j in range(N):
            distMatrix[i].append(None)

    # Initialize the index list:
    alignmentIndexList = range(alignmentLength)

    if resample:

        ####################################################
        # Resample the index list:
        tmpList = []
        for i in range(alignmentLength):
            tmpList.append(random.choice(alignmentIndexList))
        alignmentIndexList = tmpList

        # Fill in distMatrix
        for i in range(N):
            for j in range(N):
                if i < j:
                    distMatrix[i][j] = computeDistance(i, j)
                elif i > j:
                    distMatrix[i][j] = distMatrix[j][i]
                else:
                    distMatrix[i][j] = 0.0
#         ####################################################
#         if not resamplingMatrixComputed:
#             # Create resamplingMatrix:
#             for i in range(N):
#                 resamplingMatrix.append([])
#                 for j in range(N):
#                     resamplingMatrix[i].append(None)
#             # Fill in the resampling matrix if this has not been done already:
#             for i in range(N):
#                 for j in range(N):
#                     if i < j:
#                         resamplingMatrix[i][j] = k2pStats(i, j)
#                     elif i > j:
#                         resamplingMatrix[i][j] = resamplingMatrix[j][i]
#                     else:
#                         resamplingMatrix[i][j] = [0, 0, 1]
#             resamplingMatrixComputed = True
#
#         # Fill in distMatrix
#         for i in range(N):
#             for j in range(N):
#                 if i < j:
#                     distMatrix[i][j] = resampleDistance(i, j)
#                 elif i > j:
#                     distMatrix[i][j] = distMatrix[j][i]
#                 else:
#                     distMatrix[i][j] = 0.0
#         ####################################################

    else:
        # Fill in distMatrix
        for i in range(N):
            for j in range(N):
                if i < j:
                    distMatrix[i][j] = computeDistance(i, j)
                elif i > j:
                    distMatrix[i][j] = distMatrix[j][i]
                else:
                    distMatrix[i][j] = 0.0
    
    # Create Qmatrix - like alocating memory in C
    for i in range(N):
        Qmatrix.append([])
        for j in range(N):
            Qmatrix[i].append(None)
    # Fill in Qmatrix:
    for i in range(N):
        for j in range(N):
            if i <= j:
                Qmatrix[i][j] = distMatrix[i][j] - (calcRsum(i) + calcRsum(j))
            elif i > j:
                Qmatrix[i][j] = Qmatrix[j][i]

    #printQmatrix()
    #printDistMatrix()

    # Create initial N leaf nodes:
    for i in range(N):
        nodeList[i] = Node(i+1, leafSet=set([i+1]))

    # Compute the pairs that are allowed to join at first:
    computeAllowedPairs()

    # Main loop:
    while K < 2 * N - 2:
        i, j = findPair()

        if i in unconstrainedOTUs:
            # If i is not in the backbone sets we don't want i to represent the child
            # nodes. If they are both unconstrained swopping is not needed but does not matter.
            i, j = j, i 

        createParentNode(i, j)
        updateAllowedPairs(i, j)
        updateMatrices(i, j)
        
        K += 1
    assert len(allowedPairsList) == 1, allowedPairsList     # only one join remains and this should be reflected in allowedPairsList..
    n1, n2 = allowedPairsList[0]
    dij = distMatrix[n1][n2] / 2
    nodeList[K] = Node(K+1, nodeList[idxN[n2]], dij, nodeList[idxN[n1]], dij, nodeList[idxN[n2]].leafSet.union(nodeList[idxN[n1]].leafSet))
    K += 1

    # Make a string representation and return it:
    return str(nodeList[-1])
