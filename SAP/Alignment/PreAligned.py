import os, re, sys, os
from SAP.Bio.Nexus import Nexus
from SAP import Fasta
from SAP.Bio import pairwise2
from SAP.UtilityFunctions import *

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
                # Splice in the introduced gaps
                for i in range(len(introducedGaps)):
                    sequence = sequence[:introducedGaps[i]] + '-' + sequence[introducedGaps[i]:]
                # Removing irelevant flanks from sequence:
                fastaRecord.sequence = sequence[leftBoundary:rightBoundary]
                seqList.append(fastaRecord)

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

