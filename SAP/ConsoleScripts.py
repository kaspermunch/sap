
# Standard libs:
import sys, re, os, pickle, copy, time, traceback, webbrowser
from optparse import OptionParser
from math import floor
from types import ListType

# BioPython modules:
from SAP.Bio.Nexus import Nexus, Trees, Nodes

from SAP import Fasta

# Custom modules:
from SAP import MachinePool, SGE, Options
from SAP.XML2Obj import XML2Obj
from SAP.Homology import HomolCompiler, HomologySet, Homologue
from SAP.TreeStatistics import TreeStatistics
from SAP.PairWiseDiffs import PairWiseDiffs
from SAP.ResultHTML import ResultHTML
from SAP.Initialize import Initialize        
from SAP.UtilityFunctions import *
from SAP.FindPlugins import *
from SAP.InstallDependencies import assertNetblastInstalled, assertClustalw2Installed

def sap():

    try:

        optionsParser = Options.Options()
        options, args = optionsParser.postProcess()

        if options.installdependencies:
            assertNetblastInstalled()
            assertClustalw2Installed()

        if options.onlinehelp:
            webbrowser.open('http://ib.berkeley.edu/labs/slatkin/munch/StatisticalAssignmentPackage.html', new=2, autoraise=1)
            sys.exit()

        if options.viewresults:
            try:
                webbrowser.open('file://' + os.path.abspath(os.path.join(options.viewresults, 'html', 'index.html')), new=2, autoraise=1)
            except:
                if os.path.exists(options.viewresults):
                    print "The anlysis has not completed and no results are available."
                else:
                    print "The spcified project folder does not exist."
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
            assertClustalw2Installed()
    
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
                        cmd += "%s %s --_align %s ; " \
                              % (sys.argv[0], optionStr, os.path.join(options.homologcache, homologyResult.homologuesFileName))
    
    
                        print "\t\tTree sampling using", options.sampler
                        cmd += "%s %s --_sample %s ; " % (sys.argv[0], optionStr, os.path.join(options.alignmentcache, homologyResult.alignmentFileName))
            
                        # Calculation of tree statistics. (Reads pickled
                        # blastresult objects from homologcache and writes
                        # to options.treestatscache)
                        print "\t\tTree statistics computation"
                        cmd += "%s %s --_stats %s" % (sys.argv[0], optionStr, os.path.join(options.homologcache, homologyResult.homologuesPickleFileName))
    
                        if options.hostfile or options.sge:
                            try:
                                pool.enqueue(cmd)
                            except SGE.QsubError:
                                print "Error in submission of %s" % cmd
                                pass
                        else:
                            if sys.platform == 'win32':
                                # Windows CMD won't take long command lines:
                                cmds = cmd.split(';')
                                for cmd in cmds:
                                    os.system(cmd)
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

    except SystemExit, exitVal:
        sys.exit(exitVal)
    except Exception, exe: 
        print """
## SAP crached, sorry ###############################################
Help creating a more stable program by sending the debugging informaion
listed below to kaspermunch@gmail.com along with *.sap file in the project
folder and the sequence input file used.
"""
        print "".join(traceback.format_tb(sys.exc_info()[2]))
        print exe




def sap_boli_backend():

    try:

        optionsParser = Options.Options()
        options, args = optionsParser.postProcess()

        if options.installdependencies:
            assertNetblastInstalled()
            assertClustalw2Installed()

        if options.onlinehelp:
            webbrowser.open('http://ib.berkeley.edu/labs/slatkin/munch/StatisticalAssignmentPackage.html', new=2, autoraise=1)
            sys.exit()

        if options.viewresults:
            try:
                webbrowser.open('file://' + os.path.abspath(os.path.join(options.viewresults, 'html', 'index.html')), new=2, autoraise=1)
            except:
                if os.path.exists(options.viewresults):
                    print "The anlysis has not completed and no results are available."
                else:
                    print "The spcified project folder does not exist."
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

                ############################################
                # Write to step3 output file:
                baseName = os.path.splitext(os.path.split(fastaFileName)[-1])[0]
                alignmentFileName = os.path.join(options.alignmentcache, baseName + ".nex")
                alignment = Nexus.Nexus(alignmentFileName)
                alignment.export_fasta(os.path.join(options.project, "BDP3.txt"))            
                ############################################

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
            assertClustalw2Installed()
    
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

            ############################################
            step3CommandFileName = os.path.join(options.project, "step3commands.txt")
            step4CommandFileName = os.path.join(options.project, "step4commands.txt")                                

            # Flag indicating if this is the first normal call to sap:
            firstDefaultCall = None
            if os.path.exists(step3CommandFileName) or os.path.exists(step4CommandFileName):
                firstDefaultCall = False
            else:
                firstDefaultCall = True

            # Open files for writing commends or subsequent steps:
            step3CommandFile = open(step3CommandFileName, 'w')
            step4CommandFile = open(step4CommandFileName, 'w')

            step1OutputFile = open(os.path.join(options.project, "BDP1.txt"), 'w')
            ############################################

            # For each fasta file execute pipeline
            for fastaFileName in args:

                ############################################
                if firstDefaultCall is False:
                    break
                ############################################

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
                    
                    ############################################
                    if homologyResult != None:
                        step1OutputFile.write("%s\t%s\t-1\n" % (homologyResult.queryName, homologyResult.homologues.keys()))
                    ############################################

                    cmd = ''
                    if homologyResult != None:
                        # The homologyResult object serves as a job carrying the relevant information.

                        print '\tIssuing sub-tasks:'
                        # Alignment using ClustalW. (Reads the fasta files
                        # in homologcache and puts alignments in
                        # options.alignmentcache)
                        print "\t\tClustalW2 alignment"
                        ############################################
                        # Instead of executing we write the relevant commands to files:                        
                        #cmd += "%s %s --_align %s ; " \
                        cmd = "%s %s --_align %s" \
                              % (sys.argv[0], optionStr, os.path.join(options.homologcache, homologyResult.homologuesFileName))
                        step3CommandFile.write(cmd)
                        ############################################
    
                        print "\t\tTree sampling using", options.sampler
                        ############################################
                        # Instead of executing we write the relevant commands to files:                        
                        #cmd += "%s %s --_sample %s ; " % (sys.argv[0], optionStr, os.path.join(options.alignmentcache, homologyResult.alignmentFileName))
                        cmd = "%s %s --_sample %s ; " % (sys.argv[0], optionStr, os.path.join(options.alignmentcache, homologyResult.alignmentFileName))
                        step4CommandFile.write(cmd)
                        ############################################
                        
                        # Calculation of tree statistics. (Reads pickled
                        # blastresult objects from homologcache and writes
                        # to options.treestatscache)
                        print "\t\tTree statistics computation"
                        ############################################
                        # Instead of executing we write the relevant commands to files:                        
                        #cmd += "%s %s --_stats %s" % (sys.argv[0], optionStr, os.path.join(options.homologcache, homologyResult.homologuesPickleFileName))
                        cmd = "%s %s --_stats %s" % (sys.argv[0], optionStr, os.path.join(options.homologcache, homologyResult.homologuesPickleFileName))
                        step4CommandFile.write(cmd)
                        ############################################

                        ############################################
                        ## if options.hostfile or options.sge:
                        ##     try:
                        ##         pool.enqueue(cmd)
                        ##     except SGE.QsubError:
                        ##         print "Error in submission of %s" % cmd
                        ##         pass
                        ## else:
                        ##     if sys.platform == 'win32':
                        ##         # Windows CMD won't take long command lines:
                        ##         cmds = cmd.split(';')
                        ##         for cmd in cmds:
                        ##             os.system(cmd)
                        ##     else:
                        ##         os.system(cmd)
                        ############################################
    
                    # Output current status of parallel jobs
                    if options.hostfile or options.sge:
                        poolStatus(pool)
        
                    print ""

            ############################################
            # we exit here unless we are running the last analysis step: tree
            # statistics in which case we go on to wrap the whole thing up:            
            if firstDefaultCall is True:
                sys.exit()
            ############################################
            
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

            ############################################
            print "make BDP4.html"

            writeFile(os.path.join(options.project, 'BDP4.html'), '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN"><HTML><HEAD><TITLE></TITLE><META http-equiv="REFRESH" content="0;url=./html/index.html"></HEAD></HTML>')

            ############################################
            
            print 'done'


    except SystemExit, exitVal:
        sys.exit(exitVal)
    except Exception, exe: 
        print """
## SAP crached, sorry ###############################################
Help creating a more stable program by sending the debugging informaion
listed below to kaspermunch@gmail.com along with *.sap file in the project
folder and the sequence input file used.
"""
        print "".join(traceback.format_tb(sys.exc_info()[2]))
        print exe



def sap_boli():

    from optparse import OptionParser

    usage="""%prog [options] [inputFile [outputFile]]

This program runs SAP on the BOLI server"""


    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("-s", "--step",
                      dest="step",
                      type="int",
                      default=1,
                      help="Step to run")
    parser.add_option("-t", "--tempdir",
                      dest="tempdir",
                      type="string",
                      default=None,
                      help="Temp dir for tmp files")

    (options, args) = parser.parse_args()

    if options.step == 1:
        os.system("sap_boli_backend -m 1000 -q -d %s -D BOLI -S ConstrainedNJ %s " % (options.tempdir, " ".join(args)))
    elif options.step == 2:
        fp = open(os.path.join(options.tempdir, "BDP1.txt"), 'w')
        fp.write("not applicable\n")
        fp.close()
    elif options.step == 3:
        cmdFileName = os.path.join(options.tempdir, 'step3commands.txt')
        cmdFile = open(cmdFileName, 'r')
        cmd = cmdFile.read()
        os.system(cmd)
    elif options.step == 4:
        cmdFileName = os.path.join(options.tempdir, 'step4commands.txt')
        cmdFile = open(cmdFileName, 'r')
        cmd = cmdFile.read()
        os.system(cmd)
        os.system("sap_boli_backend -d %s %s" % (options.tempdir, " ".join(args)))

if __name__ == "__main__":
    main()


