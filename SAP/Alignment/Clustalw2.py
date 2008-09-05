
import os, re, sys, os
from SAP.Bio.Nexus import Nexus

from SAP.UtilityFunctions import *

class Aligner:

    def __init__(self, options):
        self.options = options
    
    def align(self, fastaFileName):
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
            commandLine = "clustalw2 -infile=%s -output=NEXUS -outfile=%s %s &> /dev/null" % (fastaFileName, outputTmpFileName, alignmentoptions)
            #commandLine = "clustalw2 -infile=%s -output=NEXUS -outfile=%s %s" % (fastaFileName, outputTmpFileName, alignmentoptions)
            os.system(commandLine)
#             stdin, stdout, stderr = os.popen3(commandLine)
#             stdin.close()
            #pipe = os.popen(commandLine, 'r')

            writeFile(tmpAlignmentFileName, readFile(outputTmpFileName))
            os.unlink(outputTmpFileName)

    
#             commandLine = Clustalw.MultipleAlignCL(fastaFileName)
#             # panalise gap ending:
#             #commandLine.is_no_end_pen = 1
#             # increase gap open panalty (default it 15)
#             commandLine.gap_open_pen = 25
#             commandLine.set_output(tmpAlignmentFileName, output_type="NEXUS")
# 
#             # Run ClustalW:
#             Clustalw.do_alignment(commandLine)

            # Move dnd file to alignment directory
            if os.path.exists(os.path.splitext(fastaFileName)[0] + ".dnd"):
                os.system("mv %s %s" % (os.path.splitext(fastaFileName)[0] + ".dnd", self.options.alignmentcache))

            # Read alignment back in as string:
            nexusContents = readFile(tmpAlignmentFileName)

            # Comment out symbols header information in nexus file
            nexusContents = re.sub("(symbols=\"[^\"]*\")", r"[\1]", nexusContents)

            # The following works around a bug in clustalw for Mac:
            nexusContents = re.search(r'.*end;\n', nexusContents, re.S).group(0)

            # Write the alignment:
            writeFile(tmpAlignmentFileName, nexusContents)

            # Read alignment back in as Nexus obj.
            alignment = Nexus.Nexus(tmpAlignmentFileName)

            # Get borders of query sequence:
            queryName = baseName
            querySeq = alignment.matrix[queryName].tostring()
            leftBound = len(re.search("^(-*)", querySeq).groups()[0])
            rightBound = len(querySeq) - len(re.search("(-*)$", querySeq).groups()[0])

            # Truncate the alignment to only include sequence relevant to the query:
            for key in alignment.matrix.keys():
                seq = alignment.matrix[key].tostring()
                alignment.matrix[key].data = seq[leftBound:rightBound]
            alignment.write_nexus_data(filename=alignmentFileName)
#write_nexus_data(self, filename=None, matrix=None, exclude=[], delete=[], blocksize=None, interleave=False, interleave_by_partition=False, comment=None, omit_NEXUS=False, append_sets=True, mrbayes=False)

            print "done", 
            sys.stdout.flush()


#             # Truncate the alignment to only include sequence relevant to the query and that all are unique:
#             uniqueSeqs = []
#             deletedKeys = []
#             for key in alignment.matrix.keys():
#                 seq = alignment.matrix[key].tostring()
#                 if seq in uniqueSeqs:
#                     deletedKeys.append(key)
#                     del alignment.matrix[key]
#                 else:
#                     uniqueSeqs.append(seq)
#                     alignment.matrix[key].data = seq[leftBound:rightBound]
#             alignment.write_nexus_data(filename=alignmentFileName)
# 
#             # Adjust constraint tree accordingly:
#             constraintTreeFileName = "%s/%s" % (self.options.homologcache, baseName + ".nex")
#             constraint = Nexus.Nexus(constraintTreeFileName)
#             rootNode = constraint.trees[0].root
#             for key in deletedKeys:
#                 constraint.trees[0].prune(key)
#             print constraint.trees[0].get_taxa(rootNode)
#             treeStr = "#NEXUS\n\nbegin trees;\n   %s\nend;\n" % constraint.trees[0].to_string()
#             writeFile(constraintTreeFileName, treeStr)

