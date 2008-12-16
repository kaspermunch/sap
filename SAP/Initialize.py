#!/usr/bin/python

try:
   import cPickle as pickle
except:
   import pickle
import sys, os, re, string, glob
from optparse import OptionParser

from SAP.Bio.Alphabet import IUPAC
from SAP.Bio.Nexus import Nexus
from UtilityFunctions import *

import Fasta

class Initialize:

    def __init__(self, options):
        self.options = options

    def createDirs(self):

        # Create cache and data directories if not already there
        if not os.path.exists(self.options.blastcache):
            os.makedirs(self.options.blastcache)
        if not os.path.exists(self.options.dbcache):
            os.makedirs(self.options.dbcache)
        if not os.path.exists(self.options.alignmentcache):
            os.makedirs(self.options.alignmentcache)
        if not os.path.exists(self.options.homologcache):
            os.makedirs(self.options.homologcache)
        if not os.path.exists(self.options.treescache):
            os.makedirs(self.options.treescache)
        if not os.path.exists(self.options.treestatscache):
            os.makedirs(self.options.treestatscache)
        if not os.path.exists(self.options.treestatscache + '/errorLogs'):
            os.makedirs(self.options.treestatscache + '/errorLogs')
        if not os.path.exists(self.options.resultdir):
            os.makedirs(self.options.resultdir)
        if not os.path.exists(self.options.resultdir + '/clones'):
            os.makedirs(self.options.resultdir + '/clones')
        if not os.path.exists(self.options.resultdir + '/stats'):
            os.makedirs(self.options.resultdir + '/stats')
        if not os.path.exists(self.options.statsdir):
            os.makedirs(self.options.statsdir)
        if not os.path.exists(self.options.datadir):
            os.makedirs(self.options.datadir)

    def fixAndMoveInput(self, files):

        allowedLetters = IUPAC.IUPACAmbiguousDNA().letters

        sequenceCount = 0
        inputFileList = []
        sequenceNameMap = {}
        for inputFileName in files:

            if not os.path.exists(inputFileName):
                print " Input filename %s does not exist!" % inputFileName
                sys.exit(1)

            if self.options.inputformat == 'nexus':
               nex = Nexus.Nexus(inputFileName)
               fileContent = ''
               for name, seq in nex.matrix.items():
                  fileContent += ">%s\n%s\n" % (name, seq.tostring())
            elif self.options.inputformat == 'fasta':
               fileContent = readFile(inputFileName)               
            else:
               print "The supported input formats are 'fasta' and 'nexus'"
               sys.exit()

            fileContent = re.sub(r'\r+', '\n', fileContent)
            tmpOutputFileName = os.path.join(self.options.datadir, "%s.tmp" % os.path.split(inputFileName)[-1])
            writeFile(tmpOutputFileName, fileContent)

            usedIDs = {}
            inputFile = open(tmpOutputFileName, 'r')
            fastaIterator = Fasta.Iterator(inputFile, parser=Fasta.RecordParser())        
            outputFileBaseName = os.path.split(inputFileName)[-1]
            newOutputFileBaseName = re.sub(r'[^0-9a-zA-Z.]+', '', outputFileBaseName)
            if re.match(r'\d', newOutputFileBaseName):
               newOutputFileBaseName = 'n' + newOutputFileBaseName

            outputFileName = os.path.join(self.options.datadir, newOutputFileBaseName)

            baseName = os.path.splitext(os.path.basename(outputFileName))[0]
            sequenceNameMap[baseName] = {}

            inputFileList.append(outputFileName)
            outputFile = open(outputFileName, 'w')
            for fastaRecord in fastaIterator:            
                sequenceCount += 1
                origName = fastaRecord.title
                fastaRecord.title = safeName(fastaRecord.title)

                # Make sure ids are unique:
                if usedIDs.has_key(fastaRecord.title.lower()): # we use lower to make sure they don't just differ in case.
                    i = 1
                    while usedIDs.has_key("%d_%s" % (fastaRecord.title, i)):
                        i += 1
                    fastaRecord.title = "%d_%s" % (fastaRecord.title, i)
                usedIDs[fastaRecord.title.lower()] = True

                sequenceNameMap[baseName][fastaRecord.title] = origName

                # Strip sequence of gap chars:
                #fastaRecord.sequence = fastaRecord.sequence.replace('-', '')
                fastaRecord.sequence = fastaRecord.sequence.replace('~', '')
                fastaRecord.sequence = fastaRecord.sequence.replace('.', '')

                # fastaRecord.sequence = re.sub('[^%s]' % allowedLetters, 'N', fastaRecord.sequence)
                fastaRecord.sequence = re.sub('[^%s-]' % allowedLetters, 'N', fastaRecord.sequence)

                
                # Print only if there is some sequence left:
                if len(fastaRecord.sequence) > 0:
                    outputFile.write(str(fastaRecord) + "\n")
            inputFile.close()
            outputFile.close()
            os.remove(tmpOutputFileName)

        return inputFileList, sequenceCount, sequenceNameMap

    def checkCacheConsistency(self, files):

        idList = []
        for inputFileName in files:
            inputFile = open(inputFileName, 'r')            
            fastaIterator = Fasta.Iterator(inputFile, parser=Fasta.RecordParser())        
            baseName = os.path.splitext(os.path.split(inputFileName)[-1])[0]
            for fastaRecord in fastaIterator:            
                idList.append("%s_%s" %  (baseName, fastaRecord.title))

        pickleFileName = os.path.join(self.options.project, os.path.split(self.options.project)[1] + '.sap')

        if os.path.exists(pickleFileName):
            # See if the options have changed - and if they have:
            # remove the selected parts of the cache.
            pickleFile = open(pickleFileName, 'r')
            prevOptions = pickle.load(pickleFile)
            pickleFile.close()

            # Lists of options that deprecates cache entries:
            deleteBlastCacheList = ["database", "maxblasthits", "limitquery", "minsignificance", "nolowcomplexfilter"]

            deleteHomologueCacheList = [ "notruncate", "quickcompile", "minidentity", "forceidentity", "subspecieslevel", "fillinall", "fillineven", "fillintomatch", "individuals", "significance", "nrsignificant",
                                         "relbitscore", "phyla", "classes", "orders", "families", "genera",
                                         "besthits", "alignmentlimit", "minimaltaxonomy", "harddiversity", "forceincludefile", "forceincludegilist", "forceexcludegilist"]
            deleteHomologueCacheList.extend(deleteBlastCacheList)

            deleteAlignmentCacheList = ["alignment", "alignmentoption"]
            deleteAlignmentCacheList.extend(deleteHomologueCacheList)

            deleteTreeStatsCacheList = ["sampler", "prunelevel"]

            print "Checking cache for deprecated entries"

            for option in self.options.__dict__.keys():

                # This serves to map between new and older option names:
                prevVersionOption = option;
                if not prevOptions.__dict__.has_key(option):
                    if option == 'dbcache':
                       prevVersionOption = 'genbankcache'
                    elif option == 'minsignificance':
                       prevVersionOption = 'evaluecutoff'
                    elif option == 'significance':
                       prevVersionOption = 'evaluesignificance'
                    elif option == 'nrsignificant':
                       prevVersionOption = 'minsignificant'
                    else:
                       print "New option for SAP \"%s\" will take default value \"%s\"" % (option, self.options.__dict__[option])

                if option in deleteBlastCacheList and self.options.__dict__[option] != prevOptions.__dict__[prevVersionOption]:
                    print '\tBlast cache'
                    for queryID in idList:
                        for entry in glob.glob(os.path.join(self.options.blastcache, queryID) + '*'):
                            print "\t\t" + os.path.split(entry)[-1]
                            os.remove(entry)
                    deleteBlastCacheList = []

                if option in deleteHomologueCacheList and self.options.__dict__[option] != prevOptions.__dict__[prevVersionOption]:
                    # Delete the homologcache for the entries in the input files:
                    print '\tHomologue and alignment cache'
                    for queryID in idList:
                        for entry in glob.glob(os.path.join(self.options.homologcache, queryID) + '*'):
                            print "\t\t" + os.path.split(entry)[-1]
                            os.remove(entry)
                        for entry in glob.glob(os.path.join(self.options.alignmentcache, queryID) + '*'):
                            print "\t\t" + os.path.split(entry)[-1]
                            os.remove(entry)
                    deleteHomologueCacheList = []

                if option in deleteTreeStatsCacheList and self.options.__dict__[option] != prevOptions.__dict__[prevVersionOption]:
                    # Delete the tree statistics cache for the entries in the input files:
                    print '\tTree statistics cache'
                    for queryID in idList:
                        for entry in glob.glob(os.path.join(self.options.treestatscache, queryID) + '*'):
                            print "\t\t" + os.path.split(entry)[-1]
                            os.remove(entry)
                    deleteTreeStatsCacheList = []

        # Dump the options specified:
        pickleFile = open(pickleFileName, 'w')
        pickle.dump(self.options, pickleFile)        
        pickleFile.close()

        print

if __name__ == "__main__":
    main()


