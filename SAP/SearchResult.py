import os, sys, re

from SAP import UtilityFunctions as utils
from SAP import NCBIXML # locally hacked to better parse string info
from SAP import NW

def _listOrTuple(x):
    return isinstance(x, (list, tuple))

def _flatten(sequence, toExpand=_listOrTuple):
    for item in sequence:
        if toExpand(item):
            for subitem in _flatten(item, toExpand):
                yield subitem
        else:
            yield item

class ResultEntry(object):
    """
    id
    input_query_letters
    sequence
    score
    subject_start
    subject_length
    query_start
    query_length
    subject_strand
    significance   
    """
    def __init__(self, hitId, input_query_letters, hspList):

        # id:
        self.id = str(hitId)

        # input length
        self.input_query_letters = input_query_letters
 
        # strand:
        if hspList[0].query_start > hspList[0].query_end:
           self.subject_strand = -1
        else:
           self.subject_strand = 1

        self.hspAvIdentity = 0.0
        self.hspCoverage = 0.0
        for i, hsp in enumerate(hspList):
#             # adjust for possible overlap to nex hsp:
#             overlap = 0
#             if i < len(hspList):
#                 overlap = max(0, self.subject_strand * (hsp.query_end - hsp.query_start))

#            self.hspAvIdentity += hsp.identities / float(hsp.sbjct_end - hsp.sbjct_start + 1 - hsp.gaps)
            self.hspAvIdentity += hsp.identities
            if hsp.gaps == (None, None): # apparently the default...
                hsp.gaps = 0
            #self.hspCoverage += hsp.sbjct_end - hsp.sbjct_start + 1 - hsp.gaps
            self.hspCoverage += self.subject_strand * (hsp.query_end - hsp.query_start) + 1

        #self.hspAvIdentity /= float(len(hspList))
        self.hspAvIdentity /= self.hspCoverage

        if self.hspCoverage > self.input_query_letters:
            for h in hspList:
                print self.id, self.hspCoverage, self.input_query_letters, h.gaps, h.query_start, h.query_end, h.sbjct_start, h.sbjct_end
            sys.exit()

        self.hspCoverage /= float(self.input_query_letters)

            
        # query start and end coordinates
        queryCoords = [x for x in _flatten([(x.query_start, x.query_end) for x in hspList])]
        queryCoords.sort()
        self.query_start = queryCoords[0] - 1 # Coords are 1-based
        self.query_length = queryCoords[-1] - queryCoords[0] + 1

        # subject start and end coordinates
        subjectCoords = [x for x in _flatten([(x.sbjct_start, x.sbjct_end) for x in hspList])]
        subjectCoords.sort()
        self.subject_start = subjectCoords[0] - 1 # Coords are 1-based
        self.subject_length = subjectCoords[-1] - subjectCoords[0] + 1

        # none sequence so far:
        self.sequence = None

        # no alignment so far:
        self.alignment = None

        # score:
        self.score = 0
        for hsp in hspList:
           self.score += hsp.bits

        # evalue
        self.significance = min([x.expect for x in hspList]) # SHOULD BE SELF.EXPECT     <--- NB

    def __str__(self):
#         s = "hit id:\t%s\n" % self.id
#         s += "score:\t%f\n" % self.score
#         s += "expect:\t%.2e\n" % self.significance
#         s += "query:\n"
#         s += "\tstart\t%d\n" % self.query_start
#         s += "\tlength\t%d\n" % self.query_length
#         s += "subject:\n"
#         s += "\tstart\t%d\n" % self.subject_start
#         s += "\tlength\t%d\n" % self.subject_length
#         s += "\tstrand\t%d\n" % self.subject_strand

        s = "%s %d %d %d, score: %.2e, min e: %.2e cov: %.2f ident: %.2f query: %d %d" % (self.id, self.subject_start, self.subject_length, self.subject_strand, self.score, self.significance, self.hspCoverage, self.hspAvIdentity, self.query_start, self.query_length)

        if self.alignment:
            width = 100
            end = width
            alnLength = len(self.alignment[0])
            s += "\nAlignment:\n"
            for start in range(0, alnLength, width):
                s += self.alignment[0][start:min(end, alnLength)] + "\n"
                s += self.alignment[1][start:min(end, alnLength)] + "\n\n"
                end += width
            s += "\n"
            
        return s

    def _similarityScore(seq1, seq2):
        assert (len(seq1) == len(seq2))

        seq1LeftBoundary = len(re.search("^(-*)", seq1).groups()[0])
        seq1RightBoundary = len(seq1) - len(re.search("(-*)$", seq1).groups()[0])
        seq2LeftBoundary = len(re.search("^(-*)", seq2).groups()[0])
        seq2RightBoundary = len(seq2) - len(re.search("(-*)$", seq2).groups()[0])
        leftBoundary = max(seq1LeftBoundary, seq2LeftBoundary)
        rightboundary = min(seq1RightBoundary, seq2RightBoundary)

        seq1Trunc = seq1[leftBoundary:rightboundary]
        seq2Trunc = seq2[leftBoundary:rightboundary]

        length = len(seq1Trunc)

        extraGaps = 0
        ident = 0
        for i in range(0, length):
            if seq1Trunc[i] == "-" and seq2Trunc[i] == "-":
                extraGaps += 1
                continue
            if seq1Trunc[i] == seq2Trunc[i]:
                ident += 1
        score = ident/float(length - extraGaps)
        return score


    def _normalizedAlignmentScore(self, gapOpen=5, gapExt=2, gapFlank=0, mismatch=4, match=5):

        seq1, seq2 = self.alignment
        assert (len(seq1) == len(seq2))

        seq1LeftBoundary = len(re.search("^(-*)", seq1).groups()[0])
        seq1RightBoundary = len(seq1) - len(re.search("(-*)$", seq1).groups()[0])
        seq2LeftBoundary = len(re.search("^(-*)", seq2).groups()[0])
        seq2RightBoundary = len(seq2) - len(re.search("(-*)$", seq2).groups()[0])
        leftBoundary = max(seq1LeftBoundary, seq2LeftBoundary)
        rightboundary = max(seq1RightBoundary, seq2RightBoundary)

        seq1Trunc = seq1[leftBoundary:rightboundary]
        seq2Trunc = seq2[leftBoundary:rightboundary]

        length = len(seq1Trunc)

        extraGaps = 0
        ident = 0
        score = - gapFlank * (leftBoundary + rightboundary)
        for i in range(0, length):

            if seq1Trunc[i] == "-" and seq2Trunc[i] == "-":
                extraGaps += 1
                continue

            if seq1Trunc[i] == "-" or seq2Trunc[i] == "-":
               if inGap is False:
                  score -= gapOpen
                  inGap = True
               else:
                  score -= gapExt
            else:
               inGap = False

               if seq1Trunc[i] == seq2Trunc[i]:
                  score += match
               else:
                  score -= mismatch

        return score / float(length - extraGaps)

    def remap(self, queryRecord, dbIndex):

        gapOpen = 15
        gapExt = 6.66
        gapFlank = 0

        queryLength = len(queryRecord.sequence)

        hitRecordLength = dbIndex.index[self.id].length

        flank = int(1.2 * (queryLength - self.query_length))
        subjectStartIndex = max(0, self.subject_start - flank)
        subjectEndIndex = min(self.subject_start + self.subject_length + flank, hitRecordLength)

        subjectCutoutLength = subjectEndIndex - subjectStartIndex
        cutout = dbIndex.get(self.id, subjectStartIndex, subjectCutoutLength, self.subject_strand)

        print len(queryRecord.sequence), len(cutout)
        alignment = utils.pairwiseClustalw2('query', queryRecord.sequence, 'subject', cutout)
        alignedQuery = alignment.matrix['query'].tostring()
        alignedSubject = alignment.matrix['subject'].tostring()
#         status, alignedQuery, alignedSubject = NW.align(queryRecord.sequence, cutout, gapOpen, gapExt, gapFlank)


#         utils._printAlignment(alignedQuery, alignedSubject)

        self.sequence = alignedSubject.replace("-", "")

        queryStartGaps = len(re.search("^(-*)", alignedQuery).groups()[0])
        queryEndGaps = len(re.search("(-*)$", alignedQuery).groups()[0])

        subjectStartGaps = len(re.search("^(-*)", alignedSubject).groups()[0])
        subjectEndGaps = len(re.search("(-*)$", alignedSubject).groups()[0])

        self.alignment = [alignedQuery, alignedSubject]
#         if queryEndGaps: # if this is zero we get and empty slice using it with minus in the slice
#             self.alignment = [alignedQuery[queryStartGaps:-queryEndGaps], alignedSubject[queryStartGaps:-queryEndGaps]]
#         else:
#             self.alignment = [alignedQuery[queryStartGaps:], alignedSubject[queryStartGaps:]]

        #self.score = self._similarityScore()
        self.score = self._normalizedAlignmentScore(gapOpen=gapOpen, gapExt=gapExt, gapFlank=gapFlank)

        self.query_start = subjectStartGaps
        self.query_length = queryLength - subjectStartGaps - subjectEndGaps

        subject_start = subjectStartIndex + queryStartGaps
        subject_end = subjectEndIndex - queryEndGaps

        self.subject_start = subject_start
        self.subject_length = subject_end - subject_start

class BlastSearchResult(object):   
    """
    Generic class for search attributes.
    """
    def __init__(self, blastRecord, realign=False):
 
       if blastRecord is None:
          return

       self.input_query_letters = blastRecord.query_letters
       
       self.hits = []
 
       # loop over all hits:
       for i in range(len(blastRecord.alignments)):
 
          # hit title:
          hitId = blastRecord.alignments[i].title.strip()
 
          # make two lists with hsps maatching on the two different strands:
          posStrandHsps = []
          negStrandHsps = []
          for hsp in blastRecord.alignments[i].hsps:
             if hsp.query_start <= hsp.query_end:
                posStrandHsps.append(hsp)
             else:
                negStrandHsps.append(hsp)
 
          # turn the hsp lists into search entries:
          for lst in [posStrandHsps, negStrandHsps]:
             # split hsp list if there is more than one hit to the same strand:
             if lst:
                 for hspList in self._splitHit(lst):
                    searchEntry = ResultEntry(hitId, self.input_query_letters, hspList)
                    self.hits.append(searchEntry)
 
       # sort by score:
       self.hits.sort(key=lambda x: (x.score, x.id), reverse=True)
 
    def _splitHit(self, hspList):

        # here it does not matter if which strand we are on, just how far the hsps are
        # from each other:
        hspList.sort(key=lambda x: (min(x.sbjct_start, x.sbjct_end), max(x.sbjct_start, x.sbjct_end)))
        #hspList.sort(key=lambda x: (x.sbjct_start, x.sbjct_end))

        # list for split hsp lists:
        splitList = []

        tmpList = []
        # we require that hsps for same hit are in running order on both query and subject
        # and that all hsps are spanned by no more than three times the query length.
        for hsp in hspList:
            if not tmpList or (abs(max(hsp.sbjct_start, hsp.sbjct_end) - min(tmpList[0].sbjct_start, tmpList[0].sbjct_end)) < self.input_query_letters * 3 and max(tmpList[-1].query_start, tmpList[-1].query_end) <= min(hsp.query_start, hsp.query_end)):
                tmpList.append(hsp)
            else:
                # start a new list of hsps:
                splitList.append(tmpList)
                tmpList = [hsp]
        if tmpList:
            splitList.append(tmpList)

        return splitList

    def __str__(self):

        s = ""
        if self.hits:
            for entry in self.hits:
                s += str(entry) + "\n"
        else:
            s += "No hits in search\n"
        return s

    def __iter__(self):
        return self.next()

    def next(self):
        for record in self.hits:
            yield record
        raise StopIteration
