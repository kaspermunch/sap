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
        Use clustalw to align fasta files
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

            # We use the fact that all homologues are pre-aligned and that the query has
            # the best match to the first homologue:

            fastaFile = open(fastaFileName, 'r')
            fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        

            # Get the query fasta record:
            queryFasta = fastaIterator.next()

            # Get the best maching homologue fasta record:
            bestMatchFasta = fastaIterator.next()

            if findOnSystem('clustalw2'):
                alignment = pairwiseClustalw2(queryFasta.title, queryFasta.sequence, \
                                              bestMatchFasta.title, bestMatchFasta.sequence)
            else:
                print "(clustalw2 not found - aligning without scoring matrix)...",
                alignmentList = pairwise2.align.globalms(queryFasta.sequence, bestMatchFasta.sequence, 2, -1, -.5, -.1)
                if alignmentList > 1:
                    print "No one best mapping..."
                    print alignmentList
                alignment = alignmentList[0]
                alignedQuerySeq, alignedBestMatch = str(alignment.matrix[bestMatchFasta.title])

            # Get the sequences:
            alignedQuerySeq = str(alignment.matrix[queryFasta.title])
            alignedBestMatch = str(alignment.matrix[bestMatchFasta.title])

            # Get the left and right boundries:
            leftBoundary = len(re.search("^(-*)", alignedQuerySeq).groups()[0])
            rightBoundary = len(alignedQuerySeq) - len(re.search("(-*)$", alignedQuerySeq).groups()[0])

            # Identify the gaps introduced into the alignmnet by aligning the query:
            introducedGaps = []
            i = 0
            j = 0
            while j < len(alignedBestMatch):
                if alignedBestMatch[j] == '-' and bestMatchFasta.sequence[i] != '-':
                    introducedGaps.append(j)
                else:
                    assert alignedBestMatch[j] ==  bestMatchFasta.sequence[i]
                    i += 1
                j += 1

            # Remove flanking gaps from query:
            queryFasta.sequence = alignedQuerySeq

            # Remove corresponding sequence from best matching homologue:
            bestMatchFasta.sequence = alignedBestMatch
 
            # Remove flanking gaps from query:
            queryFasta.sequence = alignedQuerySeq[leftBoundary:rightBoundary]

            # Remove corresponding sequence from best matching homologue:
            bestMatchFasta.sequence = alignedBestMatch[leftBoundary:rightBoundary]

            # Queue the query and aligned best homologue for writing:
            seqList = []
            seqList.append(queryFasta)
            seqList.append(bestMatchFasta)

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

