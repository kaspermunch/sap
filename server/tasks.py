from app import app, celery

from SAP.Homology import HomolCompiler
from SAP.TreeStatistics import TreeStatistics
from SAP.PairWiseDiffs import PairWiseDiffs
from SAP.ResultHTML import ResultHTML
from SAP.Initialize import Initialize
from SAP.UtilityFunctions import *
from SAP.FindPlugins import *
from SAP.InstallDependencies import *

from email_notification import email_success, email_failure, email_revoked
from celery.exceptions import SoftTimeLimitExceeded

@celery.task(name='app.run_analysis', bind=True)
def run_analysis(self, input_file, options, stdout_file, stderr_file):

    class RedirectStdStreams(object):
        def __init__(self, stdout=None, stderr=None):
            if stdout is not None:
                stdout = open(stdout, 'w')
            if stderr is not None:
                stderr = open(stderr, 'w')
            self.stdout = stdout
            self.stderr = stderr
            self._stdout = stdout or sys.stdout
            self._stderr = stderr or sys.stderr

        def __enter__(self):
            self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
            self.old_stdout.flush()
            self.old_stderr.flush()
            sys.stdout, sys.stderr = self._stdout, self._stderr

        def __exit__(self, exc_type, exc_value, traceback):
            self._stdout.flush(); self._stderr.flush()
            if sys.stdout is self.stdout:
                sys.stdout.close()
            if sys.stderr is self.stderr:
                sys.stderr.close()
            sys.stdout = self.old_stdout
            sys.stderr = self.old_stderr

    with RedirectStdStreams(stdout=stdout_file, stderr=stderr_file):

        # Make directories and write fixed inputfiles:
        init = Initialize(options)
        init.createDirs()

        inputFiles, seqCount, sequenceNameMap = init.fixAndMoveInput([input_file])
        init.checkCacheConsistency(inputFiles)

        progress = 1
        self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

        fastaFileBaseNames = []

        try:
            alignmentPlugin = findPlugin(options.alignment, 'SAP.alignment')
        except PluginNotFoundError:
            from SAP.Alignment import Clustalw2 as alignmentPlugin
            # exec("from SAP.Alignment import %s as alignmentPlugin" % options.alignment)
        aligner = alignmentPlugin.Aligner(options)

        try:
            assignmentPlugin = findPlugin(options.assignment, 'SAP.assignment')
        except PluginNotFoundError:
            if options.assignment == "Barcoder":
                from SAP.Assignment import Barcoder as assignmentPlugin
            elif options.assignment == "ConstrainedNJ":
                from SAP.Assignment import ConstrainedNJ as assignmentPlugin
            else:
                assert 0
           # exec("from SAP.Assignment import %s as assignmentPlugin" % options.assignment)
        assignment = assignmentPlugin.Assignment(options)

        uniqueDict = {}
        copyLaterDict = {}

        homolcompiler = HomolCompiler(options)

        inputQueryNames = {}

        # For each fasta file execute pipeline
        for fastaFileName in inputFiles:

            fastaFile = open(fastaFileName, 'r')
            fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())
            fastaFileBaseName = os.path.splitext(os.path.basename(fastaFileName))[0]
            fastaFileBaseNames.append(fastaFileBaseName)

            inputQueryNames[fastaFileBaseName] = {}

            for fastaRecord in fastaIterator:

                # Discard the header except for the first id word:
                fastaRecord.title = re.search(r'^(\S+)', fastaRecord.title).group(1)

                app.logger.info("file: {}, query: {}".format(fastaFileBaseName, fastaRecord.title))

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

                progress += 1
                self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

                if homologyResult != None:
                    # The homologyResult object serves as a job carrying the relevant information.

                    aligner.align(os.path.join(options.homologcache, homologyResult.homologuesFileName))

                    progress += 1
                    self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

                    try:
                       assignment.run(os.path.join(options.alignmentcache, homologyResult.alignmentFileName))
                    except assignmentPlugin.AssignmentError, X:
                       print X.msg

                    progress += 1
                    self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

                    treeStatistics = TreeStatistics(options)
                    treeStatistics.runTreeStatistics([os.path.join(options.homologcache, homologyResult.homologuesPickleFileName)], generateSummary=False)

                    progress += 1
                    self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})
                else:
                    progress += 3
                    self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

            fastaFile.close()

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
            pairwisediffs.runPairWiseDiffs(inputFiles)

        # Summary tree stats:
        print 'Computing tree statistics summary...'
        treeStatistics = TreeStatistics(options)
        treeStatistics.runTreeStatistics(inputFiles, generateSummary=True, doubleToAnalyzedDict=doubleToAnalyzedDict, inputQueryNames=inputQueryNames)
        print "done"

        progress += 1
        self.update_state(state='PROGRESS', meta={'current': progress, 'total': seqCount*4+2})

        # Make HTML output:
        print '\tGenerating HTML output...'

        resultHTML = ResultHTML(options)
        resultHTML.webify([options.treestatscache + '/summary.pickle'], fastaFileBaseNames, doubleToAnalyzedDict, sequenceNameMap)
        print 'done'

    return options.project


@celery.task(name='app.notify_email', bind=False)
def notify_email(mail, result, email_address):
    if isinstance(result, basestring):
        proj_id = os.path.basename(result)
        email_success(mail, email_address, proj_id=proj_id)
    elif isinstance(result, SoftTimeLimitExceeded):
        email_revoked(mail, email_address, proj_id=None)
    else:
        email_failure(mail, email_address, proj_id=None)
    return result
