

try:
    import cPickle as pickle
except:
    import pickle
import re, os, sys, tempfile, shutil, copy

####################################################
from SAP.Bio import Entrez
####################################################

from SAP import Fasta, Table
from SAP.Bio.Nexus import Nexus

from SAP.Databases import Native

from SAP.Assignment.ConstrainedNJ import ConstrainedNJ

import IMa
from SAP.UtilityFunctions import *

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass

class AssignmentError(Error):
    """
    Exception raised when assignment fails.
    """
    def __init__(self, message):
        self.msg = message

class Assignment:

    def __init__(self, options):
        self.options = options
        
        self.name = "IMa"

    def _getCandiateSpecies(self, fileBaseName):
        """
        Extract the two species with highest posterior probability from the assignment result.
        """
        pickleFileName = os.path.join(self.options.treestatscache, fileBaseName + ".pickle")

        pickleFile = open(pickleFileName, 'r')
        taxonomySummary = pickle.load(pickleFile)
        pickleFile.close()

        speciesTupleList = taxonomySummary.getLevelProbs('species', tupleList=[])

        if not speciesTupleList:
            return None
        
#             # load the taxonomy summary of homologues:
#             homologsPickleFileName = os.path.join(self.options.homologcache, fileBaseName + ".pickle")
#             pickleFile = open(homologsPickleFileName, 'r')
#             homologyResult = pickle.load(pickleFile)
#             pickleFile.close()
# 
#             # Get highest scoring genus:
#             genusTupleList = homologyResult.priorSummary.getLevelProbs('genus', tupleList=[])
#             if not genusTupleList:
#                 print "WARNING: no candidate genus"    
#             genusTupleList.sort(lambda x, y: cmp(x[1]), y[1])
#             genusSummary, prob = genusTupleList.pop()
# 
#             # Get all species homologs from that genus:
#             speciesTupleList = genusSummary.getLevelProbs('species', tupleList=[])
#             if speciesTupleList > 5:
#                 print "WARNING: too many species"    

        speciesTupleList.sort(cmp=lambda x,y: cmp(y[1], x[1]))

        speciesList = []
        for name, prob in speciesTupleList:
#             if len(speciesList) == 1 or prob < min(self.options.ppcutoff)/float(100):
#                 break            
            speciesList.append(name)

        if not speciesList:
            return None

        return speciesList

    def _retrieveSequences(self, speciesList):

        sequenceLists = {}
        if os.path.exists(self.options.database):
            # local database:
            db = Native.DB(self.options.database, self.options)
            for species in speciesList:
                for seqID in db.index[species]:
                    homologue, retrievalStat = db.get(seqID)
                    sequenceLists.setdefault(species, []).append(Fasta.Record(homologue.gi, homologue.sequence))
        else:
            # genbank
            for species in speciesList:

                Entrez.email = self.options.email

                handle = Entrez.esearch(db="nucleotide", retmax=10, term="%s[ORGN]" % species)
                # handle = Entrez.esearch(db="nucleotide", retmax=10, term="%s[ORGN] AND barcode[keyword]" % species)
                record = Entrez.read(handle)

                if record["Count"] < 5:
                    print "WARNING: only %d sequences representing %s" % (record["Count"], species)

                success = False

                for tries in range(10):
                    try:
                        handle = Entrez.efetch(db="nucleotide", id=','.join(record["IdList"]), rettype="fasta", retmax=10)
                        fastaIterator = Fasta.Iterator(handle, Fasta.RecordParser())
                        for entry in fastaIterator:
                            entry.sequence = re.sub('[^ATGC-]', 'N', entry.sequence)
                            entry.title = entry.title.split('|')[1]
                            sequenceLists.setdefault(species, []).append(entry)
                        success = True
                    except:
                        time.sleep(tries * 5)
                        continue
                    break

                if not success:
                    return None

        return sequenceLists


    def _alignSequences(self, sequenceLists, queryFastaRecord):

        # Make a temp dir for output files:
        tmpDirName = tempfile.mkdtemp()

        # generate file names:
        outputTmpFileName = os.path.abspath(randomString(6) + ".nex")
        fastaFileName = os.path.abspath(randomString(6) + ".fasta")

        # contatenate sequence output in the right order:
        fastaStr = ""
        speciesList = sequenceLists.keys()
        speciesList.sort()
        for species in speciesList:
            for entry in sequenceLists[species]:
                fastaStr += str(entry)

        # add the query fasta to the end and write it:
        #fastaStr += str(queryFastaRecord)
        title = queryFastaRecord.title
        fastaStr += ">%s\n%s\n" % (title[:9], queryFastaRecord.sequence)
        writeFile(fastaFileName, fastaStr)

        # call clustalw2 and read back nexus alignmnet:
        alignmentoptions = "-gapopen=50"
        commandLine = "clustalw2 -infile=%s -output=NEXUS -outfile=%s %s" % (os.path.basename(fastaFileName),
                                                                             os.path.basename(outputTmpFileName),
                                                                             alignmentoptions)
        systemCall(commandLine, stdout='IGNORE', stderr='IGNORE')
        #systemCall(commandLine)        

        alignment = Nexus.Nexus(outputTmpFileName)

        os.unlink(outputTmpFileName)
        os.unlink(fastaFileName)

        return alignment

#     def _getLargestTheta(self, sequenceLists):
#         """
#         get the largest within population pi from the diagonal of the matrix
#         """
#         largest_theta = 0
#         speciesList = sequenceLists.keys()
#         speciesList.sort()
# 
#         for species in speciesList:
#             pairs = 0
#             distSum = 0
#             for i in range(len(sequenceLists[species])):
#                 for j in range(len(sequenceLists[species][i+1:])):
#                     distSum += self._getPairwiseDiff(sequenceLists[species][i], sequenceLists[species][j])            
#                     pairs += 1
#                     
#             if distSum == 0:
#                 theta = 0
#             else:
#                 theta = distSum / float(pairs)
# 
#             if theta > largest_theta:
#                 largest_theta = theta
# 
#         return largest_theta


    def _getTheta(self, sequences):
        total = 0.0
        count = 0.0
        for i, s1 in enumerate(sequences):
            for j, s2 in enumerate(sequences[i:]):
                total += self._getPairwiseDiff(s1, s2)                                
                count += 1
        print total, count
        return total / count

    def _getAverageDiv(self, sequenceLists):
        """
        get the largest divergence between two populations.
        """
        assert len(sequenceLists) == 2
        total = 0.0
        count = 0.0
        for i, s1 in enumerate(sequenceLists[0]):
            for j, s2 in enumerate(sequenceLists[1][i:]):
                total += self._getPairwiseDiff(s1, s2)
                count += 1
        return total / count

    def _getPairwiseDiff(self, s1, s2):

        # produce random sequence names:
        n1 = randomString(6)
        n2 = randomString(6)

        # align using clustalw2:
        alignment = pairwiseClustalw2(n1, s1.sequence, n2, s2.sequence)

        # extract aligned sequences and calculate distance:
        s1Aligned = alignment.matrix[n1].tostring()
        s2Aligned = alignment.matrix[n2].tostring()
        pi = 1 - similarityScore(s1Aligned, s2Aligned)

        return pi

#     def _getDistMatrix(self, sequenceLists, queryFastaRecord):
# 
#         matrix = [[0 for x in range(len(sequenceLists) + 1)] for x in range(len(sequenceLists) + 1)] # plus one for query.
# 
#         speciesList = sequenceLists.keys()
#         speciesList.sort()
#         #speciesList.append('unknown') # the unknown is the last in the list.
#         speciesList = ['unknown'] + speciesList # the unknown is the first in the list.
#         
#         sequenceListsWithQuery = copy.copy(sequenceLists)
#         sequenceListsWithQuery['unknown'] = [queryFastaRecord]
#         
#         for i in range(len(speciesList)):
#             for j in range(i, len(speciesList)):
# 
#                 if i == j:
#                     matrix[i][j] = 0.0
#                     matrix[j][i] = 0.0
#                     continue
#                 
#                 pairs = 0
# 
#                 if len(sequenceListsWithQuery[speciesList[i]]) < len(sequenceListsWithQuery[speciesList[j]]):
#                     shortest, longest = sequenceListsWithQuery[speciesList[i]], sequenceListsWithQuery[speciesList[j]]
#                 else:
#                     longest, shortest = sequenceListsWithQuery[speciesList[i]], sequenceListsWithQuery[speciesList[j]]
# 
#                 for a in range(len(shortest)):
#                     for b in range(len(longest[a:])):
#                         matrix[i][j] += self._getDist(shortest[a], longest[b])
#                         matrix[j][i]
#                         pairs += 1
# 
#                 matrix[i][j] /= float(pairs)
#                 matrix[j][i] = matrix[i][j]
# 
#         return matrix

    def _getTree(self, sequenceLists, queryFastaRecord):

        raise "Not implemented"
    
        # what we could do here is do an NJ of all seqs using constraints to cluster
        # seqs from same species. After that we could then collapse all nodes from
        # same species to get the species topology.
        
        cnj = ConstrainedNJ()
        tree = cnj.getTreeString()
        
        return tree

    def _formatInput(self, queryFastaRecord, sequenceLists, tree, alignment):
        
        alignmentLength = None
        alignmentString = ""

        speciesList = sequenceLists.keys()
        speciesList.sort()

        # sort the sequences according to the species list:
        for species in speciesList:
            for gi in [x.title for x in sequenceLists[species]]:              
                alignmentString += "%-9s %s\n" % (gi[:9], alignment.matrix[gi].tostring())
                alignmentLength = len(alignment.matrix[gi].tostring())

        # add the query sequence to the end:
        title = queryFastaRecord.title
        alignmentString += "%-9s %s\n" % (title[:9], alignment.matrix[queryFastaRecord.title[:9]].tostring())

        generaList = [x.split(' ')[0] for x in speciesList]
        genera = "_".join(generaList)

        nameList = [x.replace(' ', '_') for x in speciesList]

        # format the input:
        inputStr = genera + "\n# For finding appropriate priors\n%d\n" % (len(speciesList) + 1) # plus one is the unknown pop.
        inputStr += "unknown " + " ".join(nameList) + "\n"
        inputStr += str(tree) + "\n"
        inputStr += "1\n%s 0" % genera
        for species in speciesList:
            inputStr += " %d" % len(sequenceLists[species])
        inputStr += " %d H 1.000000 A1\n%s" % (alignmentLength, alignmentString)
        
        return inputStr

#    def run(self, queryFastaRecord, fastaFileBaseName):
    def run(self, args):

        assert len(args) == 1

        homologPickleFileName = args[0]
        homologPickleFile = open(args[0])
        homologySet = pickle.load(homologPickleFile)
        homologPickleFile.close()

        queryFastaRecord = homologySet.queryFasta

        # Make a temp dir for output files:
        tmpDirName = tempfile.mkdtemp()

        outputFileName = os.path.join(self.options.statsdir, self.name, "%s.%s.txt" % (queryFastaRecord.title, self.name))

        print "%s: Running IMa: " % queryFastaRecord.title,

        if os.path.exists(outputFileName) and os.path.getsize(outputFileName) > 0:
            print "Using cached results."
            sys.stdout.flush()

        else:
            print "Computing...",
            sys.stdout.flush()

            # calculate required input components:

            # get the most likely species (the returned list is length one)
            candidateSpecies = self._getCandiateSpecies(queryFastaRecord.title)

            if candidateSpecies is None:
                print "no candidates...",
                return None

            # retrieve sequences from GenBank:
            sequenceLists = self._retrieveSequences(candidateSpecies)
            if sequenceLists is None:
                print "no sequenes retrieved..."
                return None

            for sList in sequenceLists.values():
                if len(sList) < 3:
                    print "not enough seqs...",
                    return None

            # compute alignment:
            alignment = self._alignSequences(sequenceLists, queryFastaRecord)

            # compute the tree topology:

            tree = "(0,1):2"
#             if len(sequenceLists) == 1:
#                 tree = "(0,1):2"
#             else:
#                 tree = self._getTree(sequenceLists, queryFastaRecord)


            if len(sequenceLists) >= 2:
                # divergence between two most supported candidate species:
                u = 2e-8
                maxSplittime = 2.0 * self._getAverageDiv(sequenceLists[candidateSpecies[0]], sequenceLists[candidateSpecies[1]]) / u # mult by two to be safe...
            else:
                maxSplittime = 20.0

            # Theta for candidate species
            maxTheta = 2.0 * self._getTheta(sequenceLists[candidateSpecies[0]]) # mult by two to be safe...
            print maxTheta
            if maxTheta == 0:
                print "maxTheta changed to 1"
                maxTheta = 1
            print "remove this"

            maxMigration = 1.0

#             # Max 4Nu:
#             #maxTheta = 5.0 * maxPi
#             maxTheta = 10.0
# 
#             # Max scaled migration rate:
#             #maxMigration = 0.0001
#             maxMigration = 1.0
#             #maxMigration = 5.0 / float(maxPi)
# 
#             # Max scaled split time:
#             maxSplittime = 60.0
#             #maxSplittime = maxTheta

            sequenceLists = { candidateSpecies[0]: sequenceLists[candidateSpecies[0]] }

            outputPrefix = os.path.join(tmpDirName, "%s.%s" % (queryFastaRecord.title, self.name))
            datFileName = os.path.join(self.options.statsdir, self.name, queryFastaRecord.title + ".dat")

            # format and write input:
            inputFileContent = self._formatInput(queryFastaRecord, sequenceLists, tree, alignment)
            writeFile(datFileName, inputFileContent)
            
            # Find best two species, get the constrainttree, all sequences for the relevant species and put them in datFileName
            #cmd = "ima2 -h"

            cmd = "ima2 -i%s -o%s -a124 -q%f -m%f -t%f -b1000000 -l10000 -d100 -z500000" % (datFileName, outputFileName, maxTheta, maxMigration, maxSplittime)

            #cmd = "ima2 -iamr1.dat -otest.out -a124 -q10.0 -m1.0 -t60.0 -s5900 -b1000000 -l10000 -d100 -z500000"
            #cmd = "ima2 -iamr1.dat -otest.out -a124 -q10.0 -m0.00001 -t100.0 -s5900 -b1000000 -l10000 -d100 -z500000"

            #cmd = r'ima2 -iC:\Users\Administrator\Desktop\amr1.dat -oC:\Users\Administrator\Desktop\test.out -a124 -q10.0 -m1.0 -t60.0 -s5900 -b10000 -l1000 -d100 -z5000'

            print cmd
            print outputPrefix
            arguments = cmd.split(' ')            
            retval = IMa.runprogram(arguments, outputPrefix)
            print retval
            
            f = open(outputFileName)
            content = f.read()
            f.close()
            candidateAssignmentProb, ghostAssignmentProb = re.findall(r'\[\d+\] (\S+)',  re.search(r'^(gene\[\d+\].*)', content, re.M).group(0))

            t = Table.Table().add_row(IMa2candidatePop=candidateAssignmentProb,
                                      IMa2ghostPop=ghostAssignmentProb,
                                      IMa2maxTheta=maxTheta,
                                      IMa2maxMigration=maxMigration,
                                      IMa2maxSplitTime=maxSplittime)
            
#            t = Table.Table().load_kv(inputTableFileName).join(t)

            tableFileName = os.path.join(self.options.statsdir, self.name, "%s.%s.tbl" % (queryFastaRecord.title, self.name))
            t.write(tableFileName)


        # Remove the tempfiles:
        shutil.rmtree(tmpDirName)

        print "done."
        sys.stdout.flush()


    def _benchmark(seq, **kwarg):


        homologPickleFileName = args[0]
        homologPickleFile = open(args[0])
        homologySet = pickle.load(homologPickleFile)
        homologPickleFile.close()

        queryFastaRecord = homologySet.queryFasta



        # Make a temp dir for output files:
        tmpDirName = tempfile.mkdtemp()

        outputFileName = os.path.join(self.options.statsdir, self.name, "%s.%s.txt" % (queryFastaRecord.title, self.name))


        # calculate required input components:

        # get the most likely species (the returned list is length one)
        candidateSpecies = self._getCandiateSpecies(queryFastaRecord.title)

        if candidateSpecies is None:
            print "no candidates...",
            return None

        # retrieve sequences from GenBank:
        sequenceLists = self._retrieveSequences(candidateSpecies)
        if sequenceLists is None:
            return None

        for sList in sequenceLists.values():
            if len(sList) < 3:
                print "not enough seqs...",
                return None

        # compute alignment:
        alignment = self._alignSequences(sequenceLists, queryFastaRecord)

        # compute the tree topology:
        tree = "(0,1):2"

        if len(sequenceLists) >= 2:
            # divergence between two most supported candidate species:
            u = 2e-8
            maxSplittime = 2.0 * self._getAverageDiv(sequenceLists[candidateSpecies[0]], sequenceLists[candidateSpecies[1]]) / u # mult by two to be safe...
        else:
            maxSplittime = 20.0

        # Theta for candidate species
        maxTheta = 2.0 * self._getTheta(sequenceLists[candidateSpecies[0]]) # mult by two to be safe...

        maxMigration = 1.0

        sequenceLists = { candidateSpecies[0]: sequenceLists[candidateSpecies[0]] }

        # format and write input:
        inputFileContent = self._formatInput(queryFastaRecord, sequenceLists, tree, alignment)
        writeFile(datFileName, inputFileContent)

        cmd = "ima2 -i%s -o%s -a124 -q%f -m%f -t%f -b1000000 -l10000 -d100 -z500000" % (datFileName, outputFileName, maxTheta, maxMigration, maxSplittime)


        arguments = cmd.split(' ')            
        retval = IMa.runprogram(arguments, outputPrefix)

        f = open(outputFileName)
        content = f.read()
        f.close()
        candidateAssignmentProb, ghostAssignmentProb = re.findall(r'\[\d+\] (\S+)',  re.search(r'^(gene\[\d+\].*)', f, re.M).group(0))


        print "Candidate species assignment", candidateAssignmentProb
        print "Ghost species assignment    ", ghostAssignmentProb


        # Remove the tempfiles:
        shutil.rmtree(tmpDirName)

        print "done."
        sys.stdout.flush()
        


# ima:
# 
# 
# How do we get a phylogeny to start with?
# 
# What is going to happen to the results if we sample rooted genalogies and not unrooted?
# 
# How do we handle the ancestral population sizes? (each is correlated with a splitting time)
# 
# Can we do the LTR like this: Compute MLE of theta given data using all genoalogies (with
# the unknown labeled as 0, 1, or 2) vs. using only using genealogies where it is labeled 1?
# 
# 
# 
# run ./configure in ima/ima to get a config.h that works on linx. The hack in imagsl.h should be enough to make it work on windows
# 
# 
# Seems we do not even need the crossplatform.h file when using EPD Python with -c mingw32




#         lowest_pi = 0
#         closest_species_idx = None
# 
#         # compare query average pi to other populations:
#         candidateSpecies = sequenceLists.keys()
#         candidateSpecies.sort()
#         for i in range(len(candidateSpecies)):
#             pairs = 0.0
#             for entry in sequenceLists[candidateSpecies[i]]:
#                species_pi +=  self._getPi(queryFastaRecord, entry)
#                pairs += 1               
#             species_pi /= pairs
# 
#             if species_pi < lowest_pi:
#                 lowest_pi = species_pi
#                 closest_species_idx = candidateSpecies[i]
# 
#         assert closest_species_idx is not None
# 
#         # the lowest betweeen pop pi to a db population should be lower than that  bewteen the two db populations.
#         assert lowest_pi < matrix[0][1]
# 
#         # db pops 0 and 1
#         # query 2
# 
#         if closest_species_idx == 0:
#             tree = "((0,2):3,1):4"
#         elif closest_species_idx == 0:
#             tree = "((1,2):3,0):4"
#         else:
#             raise Exception


#         prev_min = 10000000
#         for i in range(len(matrix)):
#             for j in range(i+1, len(matrix)):
#                 if matrix[i][j] < prev_min:
#                     prev_min = matrix[i][j]
#                     min_tuple = (str(i), str(j))
# 
#         # Find the one not in that pair:
#         outgroup = [x for x in range(len(matrix[0])) if str(x) not in min_tuple][0]
# 
#         tree = "((%s):3,%s):4" % (",".join(min_tuple), outgroup)




# Turdus                                                                                              
# # For finding appropriate priors                                                                    
# 3                                                                                                   
# amaurochalinus migratorius rufiventris                                                              
# ((0,1):3,2):4                                                                                       
# 1                                                                                                   
# Turdus 0 7 11 901 H 1.000000 A1                                                                     
# BOTW184m  ------------------------------------------------CTAATCTTCGGCGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KBNA026m  ------------------------------------------------------------------------------------------
# KBNA025m  ------------------------------------------------------------------------------------------
# TZBNA187m ---------------------------------------------TACCTAATCTTCGGCGCATGAGCCGGAATAGTGGGTACTGCCCTA
# TZBNA107m ------------------------------------------CTCTACCTAATCTTCGGCGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KBNA221m  ------------------------------------------CTCTACCTANTCTTCGGCGCATGAGCCGGAATAGTGGGTACTGCCCTA
# BOTW010m  ------------------------------------------------------------------GCCGGAATAGTGGGTACTGCCCTA
# KBARG314r ------------------------------------------CTCTATCTAATCTTCGGCGCATGGGCCGGAATAGTGGGTACTGCCCTA
# KBARG330r ------------------------------------------CTCTATCTAATCTTCGGTGCATGAGCCGGAATAGTGGGTACCGCCCTA
# KBAR621r  ------------------------------------------CTCTATCTAATCTTCGGTGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KAARG071r ---------------------------------------------------------GGTGCATGAGCCGGNATAGTAGGTACTGCCCTA
# KBAR202r  ---------------------------------------------TATCTAATCTTCGGTGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KBAR331r  ------------------------------------------CTCTATCTAATCTTCGGTGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KAARG610r ------------------------------------------CTCTATATAATCTTCGGTGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KBAR351r  ---------------------------------------------------ATCTTCGGTGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KAARG079r ------------------------------------------CTCTATATAATCTTTGGTGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KBAR780r  ------------------------------------------CTCTATCTAATCTTCGGTGCATGAGCCGGAATAGTGGGTACTGCCCTA
# KAARG550r ---------------------------------------------------------------------GGAATAGTGGGTACTGCCCTA
# KAARG203a ------------------------------------------CTCTACCTAATNTTCGGCGCATGAGCCGGAATAGTGGGTACTGCCCTA
