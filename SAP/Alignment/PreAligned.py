import os, re, sys, os
from Bio.Nexus import Nexus
from SAP import Fasta
from UtilityFunctions import *

class Aligner:

    def __init__(self, options):
        self.options = options
    
    def align(self, fastaFileName):
        """
        Just reformat from Fasta to Nexus alignment and remove redundant gap columns.
        """

        baseName = os.path.splitext(os.path.split(fastaFileName)[-1])[0]
        alignmentFileName = os.path.join(self.options.alignmentcache, baseName + ".nex")

        print "%s: Alignment: " % baseName,
        sys.stdout.flush()

        if os.path.exists(alignmentFileName) and os.path.getsize(alignmentFileName)>0:
            print "Using cached results."
            sys.stdout.flush()
        else:
            print "Computing...", 
            sys.stdout.flush()

            fastaFile = open(fastaFileName, 'r')
            fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        

            # Queue the query and aligned best homologue for writing:
            seqList = []
            for entry in fastaIterator:
                seqList.append(entry)

            # Write the remaining homologes:
            for fastaRecord in fastaIterator:
                sequence = fastaRecord.sequence
                seqList.append(fastaRecord)

            fastaFile.close()

            # Get rid of columns that are all gaps:
            gapCols = []
            seqlen = len(seqList[0].sequence)
            for i in range(seqlen):
                isGapCol = True
                for j in range(len(seqList)):
                    if seqList[j].sequence[i] != '-':
                        isGapCol = False
                        break
                if isGapCol:
                    gapCols.append(i)            
            for j in range(len(seqList)):
                for i in gapCols[::-1]:
                    seqList[j].sequence = seqList[j].sequence[:i] + seqList[j].sequence[i+1:]


            # Write the nexus alignment:
            writeNexusFile(alignmentFileName, seqList)

            print "done", 
            sys.stdout.flush()


if __name__ == "__main__":

    
    fastaFileName = sys.argv[1]
    if not os.path.exists(fastaFileName):
        raise Exception

    class Bunch(object):
        """
        Generic class for lumping attributes.
        """
        def __init__(self, **keywords):
            self.__dict__.update(keywords)

    options = Bunch(alignmentcache=".")

    aln = Aligner(options)
    aln.align(fastaFileName)

