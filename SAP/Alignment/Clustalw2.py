
import os, re, sys, os
from Bio.Nexus import Nexus

from SAP.UtilityFunctions import *

class Aligner:

    def __init__(self, options):
        self.options = options
    
    def align(self, fastaFileName, truncate=True):
        """
        Use clustalw to align fasta files
        """

        baseName = os.path.splitext(os.path.split(fastaFileName)[-1])[0]
        tmpAlignmentFileName = os.path.join(self.options.alignmentcache, baseName + ".nontrunc.nex")
        alignmentFileName = os.path.join(self.options.alignmentcache, baseName + ".nex")
        outputTmpFileName = randomString(6) + '.tmp'

        print "%s: Alignment: " % baseName,
        sys.stdout.flush()

        if os.path.exists(alignmentFileName) and os.path.getsize(alignmentFileName)>0:
            print "Using cached results."
            sys.stdout.flush()
        else:
            print "Computing...", 
            sys.stdout.flush()

            alignmentoptions = " ".join(self.options.alignmentoption)

            if os.name == 'nt':
                commandLine = "clustalw2 %s -output=NEXUS -outfile=%s %s" \
                              % (fastaFileName, outputTmpFileName, alignmentoptions)
            else:
                commandLine = "clustalw2 -infile=%s -output=NEXUS -outfile=%s %s" \
                              % (fastaFileName, outputTmpFileName, alignmentoptions)
            systemCall(commandLine, stdout='IGNORE', stderr='IGNORE')

            writeFile(tmpAlignmentFileName, readFile(outputTmpFileName))
            os.unlink(outputTmpFileName)

            # Move dnd file to alignment directory
            if os.path.exists(os.path.splitext(fastaFileName)[0] + ".dnd"):
                dndFile = os.path.splitext(fastaFileName)[0] + ".dnd"
                writeFile(os.path.join(self.options.alignmentcache, os.path.basename(dndFile)), readFile(dndFile))
                os.unlink(dndFile)
                #os.system("mv %s %s" % (os.path.splitext(fastaFileName)[0] + ".dnd", self.options.alignmentcache))

            # Read alignment back in as string:
            nexusContents = readFile(tmpAlignmentFileName)

            # Comment out symbols header information in nexus file
            nexusContents = re.sub("(symbols=\"[^\"]*\")", r"[\1]", nexusContents)

            # Write the alignment:
            writeFile(tmpAlignmentFileName, nexusContents)

            # Read alignment back in as Nexus obj.
            alignment = Nexus.Nexus(tmpAlignmentFileName)

            if truncate is True:
                # Get borders of query sequence:
                queryName = baseName
                querySeq = str(alignment.matrix[queryName])
                leftBound = len(re.search("^(-*)", querySeq).groups()[0])
                rightBound = len(querySeq) - len(re.search("(-*)$", querySeq).groups()[0])

                # Truncate the alignment to only include sequence relevant to the query:
                for key in alignment.matrix.keys():
                    seq = str(alignment.matrix[key])
                    alignment.matrix[key].data = seq[leftBound:rightBound]

            # we convet to str here to get anound wxpython unicode....
            alignment.write_nexus_data(filename=str(alignmentFileName))
            
            print "done." 
            sys.stdout.flush()


