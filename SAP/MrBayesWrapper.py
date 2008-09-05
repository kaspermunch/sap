
try:
   import cPickle as pickle
except:
   import pickle
import re, os, sys
from UtilityFunctions import *

class MrBayesWrapper:

    def __init__(self, options):
        self.options = options

    def run(self, alignmentFileName):
        """
        Call MrBayes on alignment
        Returns name of tree file
        """

        baseName = os.path.splitext(os.path.split(alignmentFileName)[-1])[0]
        treeFileName = os.path.join(self.options.treescache, baseName + ".nex")
        treeFileNamePostTree = os.path.join(self.options.treescache, baseName + ".nex.trprobs")

        print "%s: Sampling trees: " % baseName,
        sys.stdout.flush()

        if os.path.exists(treeFileNamePostTree):
            print "Using cached results."
            sys.stdout.flush()
        else:
            print "Computing...",
            sys.stdout.flush()
            commandFileName = treeFileName + ".commands"
            logFile = os.path.join(self.options.treescache, baseName + ".nex.log")

            if not os.path.exists(alignmentFileName):
                raise IOError, "Alignment file %s not found" % alignmentFileName

            # Clustalw guidance tree file:
            guidanceTreeFileName = alignmentFileName.replace(".nex", ".dnd")
            if not os.path.exists(guidanceTreeFileName):
                raise IOError, "Guidance tree file %s not found" % guidanceTreeFileName
            guidanceTree = readFile(guidanceTreeFileName).replace("\n", "")

            guidanceTree = re.sub(":[\-0-9.]+", "", guidanceTree).replace(";","")

            # Get the number of homologues to decide on a temperature for the mcmcmc:
            homologueFileName = os.path.join(self.options.homologcache, baseName + ".pickle")
            homologueFile = open(homologueFileName, 'r')
            homologyResult = pickle.load(homologueFile)
            nrHomologues = len(homologyResult.homologues)
            mcmctemp = max(0.01, 0.01 + (50 - nrHomologues) * 0.004)
            homologueFile.close()

            # print guidanceTree
            mrBayesCommands = "\n".join(("\n",
                                         "begin mrbayes;",
                                         "  set autoclose=yes nowarn=yes; ",
                                         "  execute %s; " % alignmentFileName,
                                         "  lset nst=6 rates=gamma; ",
                                         "  usertree = %s;" % guidanceTree,
                                         #"  mcmcp nruns=2 mcmcdiagn=yes diagnfreq=1000 stoprule=no stopval=0.01 savebrlen=yes filename=%s ngen=7000000 samplefreq=100 relburnin=yes burninfrac=0.04 temp=0.01 nchains=6 nswaps=1;" % treeFileName,
                                         "  mcmcp nruns=2 mcmcdiagn=yes diagnfreq=1000 stoprule=no stopval=0.01 savebrlen=no filename=%s temp=%f ngen=%d samplefreq=%d relburnin=yes burninfrac=%f; " % (treeFileName, mcmctemp, self.options.mcmcgenerations, self.options.mcmcsamplefreq, self.options.mcmcrelburnin),
                                         "  mcmc; " ,
                                         "  sumt; " ,
                                         "  quit; ",
                                         "end;\n"))

            writeFile(commandFileName, mrBayesCommands)
            os.system("mb %s &> %s" % (commandFileName, logFile))

            print "done."
            sys.stdout.flush()
