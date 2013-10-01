
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
from SAP.InstallDependencies import assertClustalw2Installed, assertBlastInstalled

from SAP.PostAnalysis import IMa

from SAP.Exceptions import AnalysisTerminated

def sap():

    try:

        optionsParser = Options.Options()
        options, args = optionsParser.postProcess()

        if options.viewresults:
            try:
                webbrowser.open('file://' + os.path.abspath(os.path.join(options.viewresults, 'html', 'index.html')), new=2, autoraise=1)
            except:
                if os.path.exists(options.viewresults):
                    print "The anlysis has not completed and no results are available."
                else:
                    print "The spcified project folder does not exist."
            sys.exit()

        if options.installdependencies:
            print "Checking that dependencies are installed on your system..."
            assertClustalw2Installed()
            assertBlastInstalled()
            sys.exit()
            
        if options.onlinehelp:
            webbrowser.open('http://ib.berkeley.edu/labs/slatkin/munch/StatisticalAssignmentPackage.html', new=2, autoraise=1)
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
                raise AnalysisTerminated(1, "The plugin or file %s was not found." % X.plugin)
            aligner = plugin.Aligner(options)
            for fastaFileName in args:
                aligner.align(fastaFileName)
        elif options._sample:
            try:
                plugin = findPlugin(options.assignment, 'sap.assignment')
            except PluginNotFoundError, X:
                raise AnalysisTerminated(1, "The plugin or file %s was not found." % X.plugin)
            assignment = plugin.Assignment(options)
            for alignmentFileName in args:
                try:
                    assignment.run(alignmentFileName)
                except plugin.AssignmentError, X:
                    print X.msg

        elif options._stats:
            treeStatistics = TreeStatistics(options)
            treeStatistics.runTreeStatistics(args, generateSummary=False)

#             #######################################
#             if options.ghostpopulation:
#                 ima = IMa.Assignment(options)
#                 ima.run(args)
#             #######################################
            
        else:

            # Check that netblast and clustalw2 are installed:
            print "Locating dependencies"
            assertClustalw2Installed()
            assertBlastInstalled()
        
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
    

            homolcompiler = HomolCompiler(options)

            inputQueryNames = {}
    
            # For each fasta file execute pipeline
            for fastaFileName in args:
        
                fastaFile = open(fastaFileName, 'r')
                fastaFileBaseName, suffix = os.path.splitext(os.path.basename(fastaFileName))
                fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())
                fastaFileBaseNames.append(fastaFileBaseName)

                inputQueryNames[fastaFileBaseName] = {}
                    
                for fastaRecord in fastaIterator:
                    #break
                    #homolcompiler = HomolCompiler(options)
    
                    # Discard the header except for the first id word:
                    fastaRecord.title = re.search(r'^(\S+)', fastaRecord.title).group(1)

                    inputQueryNames[fastaFileBaseName][fastaRecord.title] = True
    
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
                              % ('sap', optionStr, os.path.join(options.homologcache, homologyResult.homologuesFileName))
    
    
                        print "\t\tTree sampling using", options.assignment
                        cmd += "%s %s --_sample %s ; " % ('sap', optionStr, os.path.join(options.alignmentcache, homologyResult.alignmentFileName))


                        # Calculation of tree statistics. (Reads pickled
                        # blastresult objects from homologcache and writes
                        # to options.treestatscache)
                        print "\t\tTree statistics computation"
                        cmd += "%s %s --_stats %s" % ('sap', optionStr, os.path.join(options.homologcache, homologyResult.homologuesPickleFileName))

                        cmd = cmd.replace('(', '\(').replace(')', '\)')
                        
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

                fastaFile.close()
    
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
            treeStatistics.runTreeStatistics(args, generateSummary=True, doubleToAnalyzedDict=doubleToAnalyzedDict, inputQueryNames=inputQueryNames)
            print "done"

            # Make HTML output:
            print '\tGenerating HTML output...'
    
            resultHTML = ResultHTML(options)
            resultHTML.webify([options.treestatscache + '/summary.pickle'], fastaFileBaseNames, doubleToAnalyzedDict, sequenceNameMap)
            print 'done'

    except SystemExit, exitVal:
        sys.exit(exitVal)
    except AnalysisTerminated, exe:
        print "\n\n", exe.msg
        sys.exit(exe.exitValue)
#     #########################
#     except IOError, exe:
#         os.system('lsof')
#         print "".join(traceback.format_tb(sys.exc_info()[2]))
#         print exe
#     #########################
    except Exception, exe: 
        print """
## SAP crashed, sorry ###################################################
Help creating a more stable program by sending all the debugging information
between the lines and your SAP version number to kaspermunch@gmail.com along
with *.sap file in the project folder and the sequence input file used.
"""
        print "".join(traceback.format_tb(sys.exc_info()[2]))
        print exe
        print "#########################################################################"



        """
        simulate fasta seqquences using simulation pipeline
        write an estimation procedure for sap that collects IMa results in a table

        from SAP.PostAnalysis import IMa
        ima = IMa.Assignment(options)
        ima._benchmark(args)

        _benchmark should take a sequence file of for each of the two candidate species as extra arguments (instead of downloading them)

        should return a table row with results


        loop over parameter grid
        call my old simulation script
        run IMa on simulated alignment

        simulate two categories:
        two seqs from the same small/large population
        two seqs from different small/large populations with small/large split time
        assess how much it matters that we have to assume the ancestral Ne



Introduction
Introduction about assignmnet, previous work, barcoding, bold progress, nr of species in the world the problem with unsampled species/populations, our idea for addressing this, outline of the method.

Methods
About sap and the framework.
About the IMa2 model (sang chul)
How sap and ima2 interacts
About the benchmark analysis
About the example analysis on Flycatchers 

Results
Benchmark analysis
Dependence on:
diversity/divergence (probability of ILS)
how often do we make a wrong high confidence assignment with sap, and how often is this caught correctly by IMa
how reliable is a low outgroup statistic in assuring us that the true species is in the database.
When the outgroup statistic is low, how often does that mean the true species is missing?
If the outgroup statistic is high it can be so correctly because the sequence is an outgroup to the other in the db sample - which occurs more often with a small number of sequences in the database. When the outgroup statistic is misleadingly high, how often will IMa assign to the candidate species?
How dependent is IMa of the number of candidate sequences in the analysis?
What if we took all bold entries in genbank and assigned them excluding the correct species. What fraction of assignments would be high confidence wrong ones?

Example analysis


Use Coasim with no recombination to simulate homology data sets to pass to alignment, assignment and IMa
Make a benchmark_sap function in Consolescripts sap as sap but with this functionality instead of normal input
Make a way to summarize the result information



        from SimulationPipeline.SimulatePipeHotspotCTMC import simulateCoasim

        assert os.path.exists(os.path.abspath(bppseqgenOptionsFile))

        def spec(**args):
            c1 = args["N1"] / args["NeRef"]
            c2 = args["N2"] / args["NeRef"]
            c3 = args["N3"] / args["NeRef"]
            c12 = args["N12"] / args["NeRef"]
            c123 = args["N123"] / args["NeRef"]
            t1 = args["T1"] / args["g"] / (2*args["NeRef"])
            t12 = args["T12"] / args["g"] / (2*args["NeRef"])
            t123 = args["T123"] / args["g"] / (2*args["NeRef"])
            speciesTree = P(c123, M(t12, [P(c3, S(1)), P(c12, M(t1, [P(c1, S(1)), P(c2, S(1))]))]))
            return speciesTree

        seq = simulateCoasim(spec, length=L, NeRef=NeRef, r=0, g=g, u=u, addOutGroup=addOutGroup, \
                             T1=T1, T12=T12, T123=T123, \
                             N1=N1, N2=N2, N3=N3, N12=N12, N123=N123, \
                             optionsFile=os.path.abspath(bppseqgenOptionsFile))

        print >>tmpOutputFile, seq



            # Make a summary data structure    
            homologyResult = HomologySet(queryName=queryName,
                                            origQueryName=origQueryName,
                                            queryFasta=fastaRecord,
                                            blastRecord=None,
                                            homologues={},
                                            options=self.options)


            taxonomy = Taxonomy()
            taxonomy.populateFromString('family:munch,genus:kasper,species:me')
            #populateFromString(self, str, fieldsep=',', subfieldsep=':'):

            for gi in ...:
                # Make a homology object:
                homologue = Homologue(gi=gi,
                                      sequence=sequence,
                                      significance = 0, 
                                      taxonomy=taxonomy)
                homologyResult.homologues.apppend(homologue)


        """
   



if __name__ == "__main__":
    main()


