try:
   import cPickle as pickle
except:
   import pickle
import os, copy
from math import floor

import Fasta
from UtilityFunctions import *

class PairWiseDiffs:

    def __init__(self, options):
        self.options = options

    def updateAlignmentMatrix(self, matrices, sets):
        """
        Update alignment matrix read from file with potentially new sequences
        """

        for fastaFileBaseName in sets.keys():

            set = sets[fastaFileBaseName]

            # Sort the sequence names:
            keyList = set.keys()
            keyList.sort()

            # Number of keys in list:
            n = len(keyList)

            # Number of sequence pairs in set:
            pairs = n*(n-1)/2.0

            counter = 0

            # Initialise dict entry for each sequence name:
            for key in set.keys():
                if not matrices[fastaFileBaseName].has_key(key):
                    matrices[fastaFileBaseName][key] = {}

            # Keeping track of progress:
            percentage = 0.0
            oldPercentage = 0.0

            # Loop over all pairs key1 and key2 (including self pairs):
            for i in range(n):
                key1 = keyList[i]
                for j in range(i, n):
                    key2 = keyList[j]

                    counter += 1

                    # Initialise pair entry if not present:
                    if not matrices[fastaFileBaseName][key1].has_key(key2):
                        matrices[fastaFileBaseName][key1][key2] = {}
                    if not matrices[fastaFileBaseName][key2].has_key(key1):
                        matrices[fastaFileBaseName][key2][key1] = {}

                    # Trun pair entry into dict if it just contains the identity float:
                    if isinstance(matrices[fastaFileBaseName][key1][key2], float):
                        value = matrices[fastaFileBaseName][key1][key2]
                        matrices[fastaFileBaseName][key1][key2] = {}
                        matrices[fastaFileBaseName][key1][key2]['alignment'] = value
                    if isinstance(matrices[fastaFileBaseName][key2][key1], float):
                        value = matrices[fastaFileBaseName][key2][key1]
                        matrices[fastaFileBaseName][key2][key1] = {}
                        matrices[fastaFileBaseName][key2][key1]['alignment'] = value

                    # Add alignment info to the entry if this does not exist:
                    if (not matrices[fastaFileBaseName][key1][key2].has_key('score') or
                        not matrices[fastaFileBaseName][key1][key2].has_key('alignedSequences')) :

                        seq1 = set[key1]
                        seq2 = set[key2]

                        if len(seq1) >= len(seq2):
                            maximumLength = len(seq1)
                            minimumLength = len(seq2)
                        else:
                            maximumLength = len(seq2)
                            minimumLength = len(seq1)

                        #mem0 = MemoryUsage.memory()

                        tmpBaseName = randomString(8) + '.tmp'
                        fastaTmpFileName = tmpBaseName + '.fasta'
                        alignmentTmpFileName = tmpBaseName + '.nex'
                        dndTmpFileName = tmpBaseName + '.dnd'

                        writeFile(fastaTmpFileName, ">%s\n%s\n>%s\n%s\n" % ('tmp1', seq1, 'tmp2', seq2))

                        cline = Clustalw.MultipleAlignCL(fastaTmpFileName)
                        cline.set_output(alignmentTmpFileName)
                        align = Clustalw.do_alignment(cline)

                        os.remove(fastaTmpFileName)
                        os.remove(alignmentTmpFileName)
                        os.remove(dndTmpFileName)

                        aligned1 = align.get_seq_by_num(0)
                        aligned2 = align.get_seq_by_num(1)
                        assert len(aligned1) == len(aligned2)

#                         alignmentMatches = 0
#                         for i in range(len(aligned1)):
#                             if aligned1.data[i] == aligned2.data[i]:
#                                 alignmentMatches += 1 
# 
#                         normalizedScore = float(alignmentMatches)/float(maximumLength)
#                         normalizedScore_norm_shortest = float(alignmentMatches)/float(minimumLength)
# 
#                         matrices[fastaFileBaseName][key1][key2]['alignment'] = normalizedScore
#                         matrices[fastaFileBaseName][key2][key1]['alignment'] = normalizedScore
# 
#                         # alignment_norm_shortest score is normalized based on shortest sequence.
#                         # This means it gives 100% if sequences are not of the same length but one is contained in the other
#                         matrices[fastaFileBaseName][key1][key2]['alignment_norm_shortest'] = normalizedScore_norm_shortest
#                         matrices[fastaFileBaseName][key2][key1]['alignment_norm_shortest'] = normalizedScore_norm_shortest

                        simScore = similarityScore(aligned1.data, aligned2.data)
                        matrices[fastaFileBaseName][key1][key2]['score'] = simScore
                        matrices[fastaFileBaseName][key2][key1]['score'] = simScore


                        matrices[fastaFileBaseName][key1][key2]['alignedSequences'] = (aligned1.data, aligned2.data)
                        matrices[fastaFileBaseName][key2][key1]['alignedSequences'] = (aligned2.data, aligned1.data)

                        #mem1 = MemoryUsage.memory(mem0)

                    # Add taxonomy info to the entry if this does not exist:
                    if (not matrices[fastaFileBaseName][key1][key2].has_key('class') or
                        not matrices[fastaFileBaseName][key2][key1].has_key('class') or
                        not matrices[fastaFileBaseName][key1][key2].has_key('order') or
                        not matrices[fastaFileBaseName][key2][key1].has_key('order') or
                        not matrices[fastaFileBaseName][key1][key2].has_key('family') or
                        not matrices[fastaFileBaseName][key2][key1].has_key('family') or
                        not matrices[fastaFileBaseName][key1][key2].has_key('genus') or
                        not matrices[fastaFileBaseName][key2][key1].has_key('genus')):

                        # Read the most probable taxonomic level from the taxonomy statistics pickle files:
                        pickleFileName1 = "treeStatistics/%s_%s.pickle" % (fastaFileBaseName, key1)
                        if os.path.exists(pickleFileName1):
                            pickleFile1 = open(pickleFileName1, 'r')
                            try:
                                treeStat1 = pickle.load(pickleFile1)
                            except:
                                pickleFile1.close()
                                continue
                            mostProbableTaxonomyNames1 = treeStat1.mostProbableTaxonomyNames()
                            pickleFile1.close()
                        else:
                            mostProbableTaxonomyNames1 = {'class':{}, 'order':{}, 'family':{}, 'genus':{}}

                        pickleFileName2 = "treeStatistics/%s_%s.pickle" % (fastaFileBaseName, key2)
                        if os.path.exists(pickleFileName2):
                            pickleFile2 = open(pickleFileName2, 'r')
                            try:
                                treeStat2 = pickle.load(pickleFile2)
                            except:
                                pickleFile2.close()
                                continue
                            mostProbableTaxonomyNames2 = treeStat2.mostProbableTaxonomyNames()
                            pickleFile2.close()
                        else:
                            mostProbableTaxonomyNames2 = {'class':{}, 'order':{}, 'family':{}, 'genus':{}}

                        # Identify the names of the most probable class, order, family, genus:
                        if ((not matrices[fastaFileBaseName][key1][key2].has_key('class')) or
                            (not matrices[fastaFileBaseName][key2][key1].has_key('class'))):
                            mostProbableClass1 = mostProbableTaxonomyNames1['class'].get('name')
                            mostProbableClass2 = mostProbableTaxonomyNames2['class'].get('name')
                            if (mostProbableClass1 != None) and (mostProbableClass2 != None):
                                if mostProbableClass1 == mostProbableClass2:
                                    matrices[fastaFileBaseName][key1][key2]['class'] = 1
                                    matrices[fastaFileBaseName][key2][key1]['class'] = 1
                                else:
                                    matrices[fastaFileBaseName][key1][key2]['class'] = 0
                                    matrices[fastaFileBaseName][key2][key1]['class'] = 0

                        if ((not matrices[fastaFileBaseName][key1][key2].has_key('order')) or
                            (not matrices[fastaFileBaseName][key2][key1].has_key('order'))):
                            mostProbableOrder1 = mostProbableTaxonomyNames1['order'].get('name')
                            mostProbableOrder2 = mostProbableTaxonomyNames2['order'].get('name')
                            if (mostProbableOrder1 != None) and (mostProbableOrder2 != None):
                                if mostProbableOrder1 == mostProbableOrder2:
                                    matrices[fastaFileBaseName][key1][key2]['order'] = 1
                                    matrices[fastaFileBaseName][key2][key1]['order'] = 1
                                else:
                                    matrices[fastaFileBaseName][key1][key2]['order'] = 0
                                    matrices[fastaFileBaseName][key2][key1]['order'] = 0

                        if ((not matrices[fastaFileBaseName][key1][key2].has_key('family')) or
                            (not matrices[fastaFileBaseName][key2][key1].has_key('family'))):
                            mostProbableFamily1 = mostProbableTaxonomyNames1['family'].get('name')
                            mostProbableFamily2 = mostProbableTaxonomyNames2['family'].get('name')
                            if (mostProbableFamily1 != None) and (mostProbableFamily2 != None):
                                if mostProbableFamily1 == mostProbableFamily2:
                                    matrices[fastaFileBaseName][key1][key2]['family'] = 1
                                    matrices[fastaFileBaseName][key2][key1]['family'] = 1
                                else:
                                    matrices[fastaFileBaseName][key1][key2]['family'] = 0
                                    matrices[fastaFileBaseName][key2][key1]['family'] = 0

                        if ((not matrices[fastaFileBaseName][key1][key2].has_key('genus')) or
                            (not matrices[fastaFileBaseName][key2][key1].has_key('genus'))):
                            mostProbableGenus1 = mostProbableTaxonomyNames1['genus'].get('name')
                            mostProbableGenus2 = mostProbableTaxonomyNames2['genus'].get('name')
                            if (mostProbableGenus1 != None) and (mostProbableGenus2 != None):
                                if mostProbableGenus1 == mostProbableGenus2:
                                    matrices[fastaFileBaseName][key1][key2]['genus'] = 1
                                    matrices[fastaFileBaseName][key2][key1]['genus'] = 1
                                else:
                                    matrices[fastaFileBaseName][key1][key2]['genus'] = 0
                                    matrices[fastaFileBaseName][key2][key1]['genus'] = 0


                        if matrices[fastaFileBaseName][key1][key2].has_key('class'):
                            classScore = matrices[fastaFileBaseName][key1][key2]['class']
                        else:
                            classScore = None

                        if matrices[fastaFileBaseName][key1][key2].has_key('order'):
                            orderScore = matrices[fastaFileBaseName][key1][key2]['order']
                        else:
                            orderScore = None

                        if matrices[fastaFileBaseName][key1][key2].has_key('family'):
                            familyScore = matrices[fastaFileBaseName][key1][key2]['family']
                        else:
                            familyScore = None

                        if matrices[fastaFileBaseName][key1][key2].has_key('genus'):
                            genusScore = matrices[fastaFileBaseName][key1][key2]['genus']
                        else:
                            genusScore = None

                    if pairs:
                        percentage = floor(counter/float(pairs) * 100)
                        if percentage and not percentage % 10 and percentage != oldPercentage:
                            print "\t%s: %.0f%%" % (fastaFileBaseName, percentage)
                            oldPercentage = percentage

                # Dump the updated version of the matrix:
                alignmentMatrixFileName = os.path.join(self.options.statsdir, fastaFileBaseName + "_matrix.pickle")
                alignmentMatrixFile = open(alignmentMatrixFileName, 'w')
                pickle.dump(matrices[fastaFileBaseName], alignmentMatrixFile)
                alignmentMatrixFile.close()

    #         # Dump the updated version of the matrix:
    #         alignmentMatrixFile = open(alignmentMatrixFileName, 'w')
    #         pickle.dump(matrices[fastaFileBaseName], alignmentMatrixFile)
    #         alignmentMatrixFile.close()


    def runPairWiseDiffs(self, fastaFileNames):

        print 'Calculating pairwise diffenrences:'

    #     # Read in fasta sequences into a dictionary:
    #     completeSets = {}
    #     for fastaFileName in fastaFileNames:
    #         baseName = os.path.splitext(os.path.basename(fastaFileName))[0]
    #         #baseName = os.path.basename(fastaFileName).split(".")[0]
    #         completeSets[baseName] = {}
    # 
    #         fastaFile = open(fastaFileName, 'r')
    #         fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        
    #         for fastaRecord in fastaIterator:
    #             newName = safeName(copy.copy(fastaRecord.title))
    #             #completeSets[baseName][fastaRecord.title.strip()] = fastaRecord.sequence
    #             completeSets[baseName][newName] = fastaRecord.sequence
    #         fastaFile.close()
    # 
    #     # Load existing alignment matrix
    #     alignmentMatrices = {}
    #     for fastaFileBaseName in completeSets.keys():
    #         if not alignmentMatrices.has_key(fastaFileBaseName):
    #             alignmentMatrices[fastaFileBaseName] = {}
    #         
    #         alignmentMatrixFileName = os.path.join(options.statsdir, fastaFileBaseName + "_matrix.pickle")
    #         if os.path.exists(alignmentMatrixFileName) and os.path.getsize(alignmentMatrixFileName)>0:
    #             alignmentMatrixFile = open(alignmentMatrixFileName, 'r')
    #             alignmentMatrices[fastaFileBaseName] = pickle.load(alignmentMatrixFile)
    #             alignmentMatrixFile.close()
    # 
    #     # Add any new alignments to alignment matrix (and save to them to file)
    #     self.updateAlignmentMatrix(alignmentMatrices, completeSets, alignmentMatrixFileName)

        # Read in fasta sequences into a dictionary:
        completeSets = {}
        for fastaFileName in fastaFileNames:
            baseName = os.path.splitext(os.path.basename(fastaFileName))[0]
            #baseName = os.path.basename(fastaFileName).split(".")[0]
            completeSets[baseName] = {}

            fastaFile = open(fastaFileName, 'r')
            fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        
            for fastaRecord in fastaIterator:
                newName = safeName(copy.copy(fastaRecord.title))
                #completeSets[baseName][fastaRecord.title.strip()] = fastaRecord.sequence
                completeSets[baseName][newName] = fastaRecord.sequence
            fastaFile.close()

        # Load existing alignment matrix
        alignmentMatrices = {}
        for fastaFileBaseName in completeSets.keys():
            if not alignmentMatrices.has_key(fastaFileBaseName):
                alignmentMatrices[fastaFileBaseName] = {}

            alignmentMatrixFileName = os.path.join(self.options.statsdir, fastaFileBaseName + "_matrix.pickle")
            if os.path.exists(alignmentMatrixFileName) and os.path.getsize(alignmentMatrixFileName)>0:
                alignmentMatrixFile = open(alignmentMatrixFileName, 'r')
                alignmentMatrices[fastaFileBaseName] = pickle.load(alignmentMatrixFile)
                alignmentMatrixFile.close()

        # Add any new alignments to alignment matrix (and save to them to file)
        self.updateAlignmentMatrix(alignmentMatrices, completeSets)
