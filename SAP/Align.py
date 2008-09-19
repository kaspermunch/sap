
import os, re, sys, os
from SAP.Bio.Nexus import Nexus
from SAP.Bio import pairwise2, Seq

import Fasta

from UtilityFunctions import *

class Align:

    def __init__(self, options):
        self.options = options
    
    def align(self, fastaFileName):


        baseName = os.path.splitext(os.path.split(fastaFileName)[-1])[0]
        alignmentFileName = os.path.join(self.options.alignmentcache, baseName + ".nex")

        print "%s: Alignment: " % baseName,
        sys.stdout.flush()

        if os.path.exists(alignmentFileName) and os.path.getsize(alignmentFileName) > 0:
            print "Using cached results."
            sys.stdout.flush()
        else:
            print "Computing...", 
            sys.stdout.flush()

            # Read in all the fasta entries:
            fastaFile = open(fastaFileName, 'r')
            fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        
            fastaDict = {}
            for fastaRecord in fastaIterator:
                fastaDict[fastaRecord.title] = fastaRecord.sequence

            # Get the query sequence:
            queryName = baseName
            querySeq = fastaDict[queryName]
            del fastaDict[queryName]

            alignment = Nexus.Nexus()
    
            seqList = [[queryName, Seq.Seq(querySeq)]]
            for title, sequence in fastaDict.items():
                alignmentList = pairwise2.align.globalms(querySeq, sequence, 1, 0, -10, -.5, one_alignment_only=1)
                queryAlnString = alignmentList[0][0]
                homologueAlnString = alignmentList[0][1]

                # Delete the columns that introcudes gaps in the querySeq:
                deleteList = []
                for i, char in enumerate(queryAlnString):
                    if char == '-':
                        deleteList.append(i)
                prunedHomologueAlnString = ''
                for i, char in enumerate(homologueAlnString):
                    if i not in deleteList:
                        prunedHomologueAlnString += char
                homologueAlnString = prunedHomologueAlnString

                # Make a list of tuple with truncated title and Seq object and add it the list:
                tup = (title, Seq.Seq(homologueAlnString))
                seqList.append(tup)

            # Write the alignment to a file:
            writeNexusFile(alignmentFileName, seqList)
            
            print "done", 
            sys.stdout.flush()


if __name__ == "__main__":

    from optparse import OptionParser

    usage="""%prog [options] [inputFile [outputFile]]

This program ....
It also does ..."""


    parser = OptionParser(usage=usage, version="%prog 1.0")

    parser.add_option("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      #type="string",
                      default=False,
                      help="Print status output to STDERR")

    (options, args) = parser.parse_args()

    inFile = args.pop(0)

    aligner = Align(None)
    aligner.align(inFile)





