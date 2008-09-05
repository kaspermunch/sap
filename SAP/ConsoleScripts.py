
# Standard libs:
import sys, re, os, pickle, copy, time
from optparse import OptionParser
from math import floor
from types import ListType

# BioPython modules:
from Bio.Nexus import Nexus, Trees, Nodes

import Fasta

# Custom modules:
import SAP
# import NeighbourJoin
import MachinePool
import SGE
import Options
from XML2Obj import XML2Obj
from Homology import HomolCompiler, HomologySet, Homologue
from TreeStatistics import TreeStatistics
# from ClustalWrapper import ClustalWrapper
# from MrBayesWrapper import MrBayesWrapper
# from BarcoderWrapper import BarcoderWrapper
# from NeighbourJoinWrapper import NeighbourJoinWrapper
from PairWiseDiffs import PairWiseDiffs
from ResultHTML import ResultHTML
from Initialize import Initialize        
from UtilityFunctions import *
from InstallDependencies import assertNetblastInstalled

def sap():

    optionsParser = Options.Options()
    options, args = optionsParser.postProcess()

    if options.onlinehelp:
        import webbrowser
        webbrowser.open('http://people.binf.ku.dk/kasper/wiki/SAP.html', new=2, autoraise=1)
        sys.exit()

    # Make a string of all options except of the ones with a '_'
    # prefix which are for internal use only:
    optionStr = ''
    for k, v in options.__dict__.items():
        if not re.match('_', k):
            if type(v) is ListType:
                if len(v) > 0:
                    joinString = ' --%s ' % k
                    v = [str(x) for x in v]
                    optionStr += joinString + joinString.join(v)
            elif v is True:
                optionStr += ' --%s' % k
            elif v is not False:
                optionStr += ' --%s %s' % (k, v)

    if options._align:
        try:
            plugin = findPlugin(options.alignment, 'sap.alignment')
        except PluginNotFoundError, X:
            print "The plugin or file %s was not found." % X.plugin
        aligner = plugin.Aligner(options)
        for fastaFileName in args:
            aligner.align(fastaFileName)
    elif options._sample:
        try:
            plugin = findPlugin(options.sampler, 'sap.sampler')
        except PluginNotFoundError, X:
            print "The plugin or file %s was not found." % X.plugin
        sampler = plugin.Sampler(options)
        for alignmentFileName in args:
            sampler.run(alignmentFileName)        
    elif options._stats:
        treeStatistics = TreeStatistics(options)
        treeStatistics.runTreeStatistics(args, generateSummary=False)
    else:

        # Check that netblast and clustalw2 are installed:
        print "Locating dependencies"
        assertNetblastInstalled()
        if not findInSystemPath('clustalw2'):
            print 'clustalw2 is not installed'
            sys.exit()

        # Make directories and write fixed inputfiles:
        init = Initialize(options)
        init.createDirs()

        # Fix the format of input files and copy them to the project directory:
        args, seqCount, sequenceNameMap = init.fixAndMoveInput(args)
        if not args or not seqCount:
            print "You need to specify file name of at least one non-empty sequence file in Fasta format."
            sys.exit()
                    
        # Make sure cache is consistent with specified options:
        init.checkCacheConsistency(args)

        if options.hostfile:
            pool = MachinePool.MachinePool(options.hostfile)
        elif options.sge:
            pool = SGE.SGE(nodes=options.sge)

        fastaFileBaseNames = []

        uniqueDict = {}
        copyLaterDict = {}


#         def getAllFromQueue(self, Q):
#             """Generator to yield one after the others all item currently in
#             the queue Q, without any waiting"""
#             try:
#                 while True:
#                     yield Q.get_nowait()
#             except Queue.Empty:
#                 raise StopIteration
# 
#         import Queue
#         outputQueue = Queue.Queue()
#         maxNrOfThreads = 5
#         homologyThreadPool = HomologyThreadPool(options, outputQueue, maxNrOfThreads)
# 
#         # For each fasta file execute pipeline
#         for fastaFileName in args:
#             records += 1
#             homologyThreadPool.put([fastaRecord, fastaFileName])
#             homologyPoolStatus(homologyThreadPool)
# 
#             # Homology.py would need the same functionality that
#             # SGE.py has so that the queue can be monitored and it can
#             # be made sure that at most some number of threads are
#             # running. The run method should then call compileHomologueset()
# 
#         
#         while records:
#             for record in self.getAllFromQueue(outputQueue):
#                 records -= 1
#                 # submit or run sub commands:
#             
#                 homologyPoolStatus(homologyThreadPool)
#                 
#                 time.sleep(1)
# 
#         # Close the thread pool:
#         homologyThreadPool.close()
#         homologyPoolStatus(homologyThreadPool)


        # For each fasta file execute pipeline
        for fastaFileName in args:
    
            fastaFile = open(fastaFileName, 'r')
            fastaFileBaseName, suffix = os.path.splitext(os.path.basename(fastaFileName))
            fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())
            fastaFileBaseNames.append(fastaFileBaseName)

            for fastaRecord in fastaIterator:

                homolcompiler = HomolCompiler(options)

                # Discard the header except for the first id word:
                fastaRecord.title = re.search(r'^(\S+)', fastaRecord.title).group(1)

                print "%s -> %s: " % (fastaFileBaseName, fastaRecord.title)
                
                # See if the sequence is been encountered before and if so skip it for now:
                if uniqueDict.has_key(fastaRecord.sequence):
                    copyLaterDict.setdefault(uniqueDict[fastaRecord.sequence], []).append('%s_%s' % (fastaFileBaseName, fastaRecord.title))
                    print '\tsequence double - skipping...\n'
                    continue
                else:
                    uniqueDict[fastaRecord.sequence] = '%s_%s' % (fastaFileBaseName, fastaRecord.title)

                # Find homologues: Fasta files and pickled homologyResult objects are written to homologcache
                homologyResult = homolcompiler.compileHomologueSet(fastaRecord, fastaFileBaseName)

                cmd = ''
                if homologyResult != None:
                    # The homologyResult object serves as a job carrying the relevant information.

                    print '\tIssuing sub-tasks:'
                    # Alignment using ClustalW. (Reads the fasta files
                    # in homologcache and puts alignments in
                    # options.alignmentcache)
                    print "\t\tClustalW2 alignment"
                    cmd += "sap %s --_align %s ; " \
                          % (optionStr, os.path.join(options.homologcache, homologyResult.homologuesFileName))


                    print "\t\tTree sampling using", options.sampler
                    cmd += "sap %s --_sample %s ; " % (optionStr, os.path.join(options.alignmentcache, homologyResult.alignmentFileName))
        
                    # Calculation of tree statistics. (Reads pickled
                    # blastresult objects from homologcache and writes
                    # to options.treestatscache)
                    print "\t\tTree statstics computation"
                    cmd += "sap %s --_stats %s" % (optionStr, os.path.join(options.homologcache, homologyResult.homologuesPickleFileName))

                    if options.hostfile or options.sge:
                        try:
                            pool.enqueue(cmd)
                        except SGE.QsubError:
                            print "Error in submission of %s" % cmd
                            pass
                    else:
                        os.system(cmd)

                # Output current status of parallel jobs
                if options.hostfile or options.sge:
                    poolStatus(pool)
    
                print ""

        if options.hostfile or options.sge:
            # Wait for all jobs to finish:
            pool.close()
            # Output current status of parallel jobs
            poolStatus(pool)
            print "\tPool closed"

        # Make dictionary to map doubles the ones analyzed:
        doubleToAnalyzedDict = {}
        for k, l in copyLaterDict.items():
            doubleToAnalyzedDict.update(dict([[v,k] for v in l]))

        if not options.nocopycache and len(doubleToAnalyzedDict):
            # Copy cache files for sequences that occoured more than once:
            print "Copying cached results for %d doubles" % len(doubleToAnalyzedDict)
            copyCacheForSequenceDoubles(copyLaterDict, options)
            
        # Calculate the pairwise differences between sequences in each file:
        if options.diffs:
            pairwisediffs = PairWiseDiffs(options)
            pairwisediffs.runPairWiseDiffs(args)
            #runPairWiseDiffs(args)

        # Summary tree stats:
        print 'Computing tree statistics summary...'
        treeStatistics = TreeStatistics(options)
        treeStatistics.runTreeStatistics(args, generateSummary=True, doubleToAnalyzedDict=doubleToAnalyzedDict)
        print "done"

        # Make HTML output:
        print '\tGenerating HTML output...'

        resultHTML = ResultHTML(options)
        resultHTML.webify([options.treestatscache + '/summary.pickle'], fastaFileBaseNames, doubleToAnalyzedDict, sequenceNameMap)
        print 'done'
