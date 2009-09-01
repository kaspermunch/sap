
try:
   import cPickle as pickle
except:
   import pickle
import re, os, sys, tempfile, shutil

# from SAP.Bio.Nexus import Nexus

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
        self.message = message

class Assignment:

    def __init__(self, options):
        self.options = options

        self.name = "IMa"

    def run(self, alignmentFileName):

        # Make a temp dir for output files:
        tmpDirName = tempfile.mkdtemp()

        baseName = os.path.splitext(os.path.split(alignmentFileName)[-1])[0]        
        outputPrefix = os.path.join(tmpDirName, "%s.%s" % (baseName, self.name))
        datFileName = os.path.join(tmpDirName, baseName + ".dat")
        outputFileName = os.path.join(tmpDirName, baseName + ".out")

        print "%s: Running IMa: " % baseName,
        sys.stdout.flush()

        if os.path.exists(treesFileName) and os.path.getsize(treesFileName) > 0:
            print "Using cached results."
            sys.stdout.flush()
        else:
            print "Computing...",
            sys.stdout.flush()

            # Find best two species, get the constrainttree, all sequences for the relevant species and put them in datFileName

            #cmd = "ima2 -iamr1.dat -otest.out -a124 -q10.0 -m1.0 -t60.0 -s5900 -b1000000 -l10000 -d100 -z500000"
            cmd = r'ima2 -iC:\Users\Administrator\Desktop\amr1.dat -oC:\Users\Administrator\Desktop\test.out -a124 -q10.0 -m1.0 -t60.0 -s5900 -b10000 -l1000 -d100 -z5000'



            arguments = cmd.split(' ')            
            retval = IMa.runprogram(arguments, outputPrefix)

        # Remove the tempfiles:
        shutil.rmtree(tmpDirName)

        print "done."
        sys.stdout.flush()

