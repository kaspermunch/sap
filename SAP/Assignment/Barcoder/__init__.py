
try:
   import cPickle as pickle
except:
   import pickle
import re, os, sys, tempfile, shutil

from SAP.Bio.Nexus import Nexus

import Barcoder
from SAP.UtilityFunctions import *

class Assignment:

    def __init__(self, options):
        self.options = options

        self.name = "Barcoder"

    def run(self, alignmentFileName):
        """
        Call Maxlike on alignment
        Returns name of tree file
        """

        # Make a temp dir for output files:
        tmpDirName = tempfile.mkdtemp()

        baseName = os.path.splitext(os.path.split(alignmentFileName)[-1])[0]        
        outputPrefix = os.path.join(tmpDirName, "%s.%s" % (baseName, self.name))
        constraintsFileName = os.path.join(tmpDirName, baseName + ".constr")
        phylipFileName = os.path.join(tmpDirName, baseName + ".phylip")
        treesFileName = os.path.join(self.options.treescache, "%s.%s.nex" % (baseName, self.name))

        print "%s: Sampling trees: " % baseName,
        sys.stdout.flush()

        if os.path.exists(treesFileName) and os.path.getsize(treesFileName) > 0:
            print "Using cached results."
            sys.stdout.flush()
        else:
            print "Computing...",
            sys.stdout.flush()

            constraintTreeFileName = os.path.join(self.options.homologcache, baseName + ".nex")
            constraintTreeFile = open(constraintTreeFileName, 'r')
            constraintTree = Nexus.Nexus(constraintTreeFile).trees[0]

            if not os.path.exists(alignmentFileName):
                raise IOError, "Alignment file %s not found" % alignmentFileName

            # Produce a list of groups that needs to be enforced in each sampled tree:
            rootNode=constraintTree.root
            #constraintList = self.constraintPartitions(constraintTree, rootNode)
            constraintList = self.constraintPartitions(constraintTree, rootNode, partList=[])
            #constraintList = filter(lambda x: len(x) > 1, constraintList)
            constraintList.sort(lambda a, b: cmp(len(a),len(b)))            

            pruneIndexes = []
            for i in range(len(constraintList) - 1):   # (consstraintList[-1] contains all nodes)
                # Remove partitions that are defined by another partition as partition.difference(allnodes):
                for j in range(i):
                    if constraintList[i] == constraintList[-1].difference(constraintList[j]):
                        pruneIndexes.append(i)
                # remove non-sensible constraints of len 1 or all len otus - 1:
                if len(constraintList[i]) == 1 or len(constraintList[i]) == len(constraintList[-1]) - 1:
                    if i not in pruneIndexes:
                        pruneIndexes.append(i)                    
            for idx in pruneIndexes:
                del constraintList[idx]
            
            # Make a list of partition constraint strings of the form: a b c | d e f g
            constraintPartitionList = []
            for i in range(len(constraintList) - 1):
                assert constraintList[-1].issuperset(constraintList[i])                   
                partitionList = []
                for otus in constraintList[i]:
                    partitionList.append(otus)
                complementList = []
                for otus in constraintList[-1].difference(constraintList[i]):
                    complementList.append(otus)
                constraintPartitionList.append("%s | %s" % (" ".join(partitionList), " ".join(complementList)))
            # Write the constraints:
            writeFile(constraintsFileName, '\n'.join(constraintPartitionList) + '\n')

            # Write a phylip formatted alignment:
            alignment = Nexus.Nexus(alignmentFileName)
            seqList = alignment.matrix.items()
            
            writePhylipFile(phylipFileName, seqList)

            # Run Barcoder:
            chainLength = 1100000
            printFreq = 1000
            sampleFreq = 1000
            burninLength = 100000
               
            cmd = "bc1 -i %s -c %s -o %s -l %d -p %d -s %d" % (phylipFileName, constraintsFileName, outputPrefix,
                                                               chainLength, printFreq, sampleFreq)
            arguments = cmd.split(' ')
            
            retval = Barcoder.runprogram(arguments, outputPrefix)

            ## outFile = os.path.join(tmpDirName, baseName + ".out")
            ## retval = os.system("bc1 -i %s -c %s -o %s -l %d -p %d -s %d &> %s" % (phylipFileName, constraintsFileName, outputPrefix,
            ##                                                     self.options.mcmcgenerations, printFreq, self.options.mcmcsamplefreq, outFile))
            ## if retval == 0:
            ##     # Write a flag file to indicate the run completed properly:
            ##     writeFile(flagFile, '')
            ## elif os.path.exists(flagFile):
            ##     # Remove a previous flag file if it should exist:
            ##     os.remove(flagFile)            

            if retval == 0:
               # If successfull, write the trees to the cache without the burnin:
               treesFile = open(treesFileName, 'w')
               tmpTreesFile = open(outputPrefix + '.tree', 'r')
               burninTrees = burninLength / sampleFreq + 1 # we add one because the very first tree is always sampled.
               for line in tmpTreesFile.xreadlines():
                  if burninTrees and line.startswith('   tree'):
                      burninTrees -= 1
                      continue
                  treesFile.write(line)
            else:
               # Move all the output files to the cache dir for inspection:
               shutil.move(outputPrefix + '.tree', self.options.treescache)
               print "Sampler: %s failed. Output files moved to trees cache for inspection" % self.name

            # Remove the tempfiles:
            shutil.rmtree(tmpDirName)

            print "done."
            sys.stdout.flush()

    #def constraintPartitions(self, tree, nodeID, partList=[]):
    def constraintPartitions(self, tree, nodeID, partList):
        """
        Creates a list of all partitions in the constraint tree.
        """
        node = tree.chain[nodeID]
        if node.succ != []:
            otus = tree.get_taxa(nodeID)
            partList.append(set(otus))
            for childID in node.succ:
                self.constraintPartitions(tree, childID, partList=partList)
        return partList

