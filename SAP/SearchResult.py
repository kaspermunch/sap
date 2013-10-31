import os, sys

from SAP import UtilityFunctions as util
from SAP.Bio.Blast import NCBIXML


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
   score
   subject_start
   subject_length
   query_start
   query_length
   query_strand
   significance
   """
   def __init__(self, **keywords):
      self.__dict__.update(keywords)


class BlastSearchResult(object):   
   """
   Generic class for search attributes.
   """
   def __init__(self, blastRecord, realign=False):

      self.search_list = []

      if blastRecord is None:
         return

      for i in range(len(blastRecord.alignments)):

         searchEntry = ResultEntry()

         searchEntry.id = blastRecord.alignments[i].title.strip()



#          hspCoords = []
#          queryCoords = []
# 
#          searchEntry.score = 0
#          hsp = blastRecord.alignments[i].hsps[0]
#          searchEntry.score += hsp.bits
#          hspCoords += [hsp.sbjct_start, hsp.sbjct_end]
#          queryCoords.append([hsp.query_start, hsp.query_end])
# 
#          hspCoords.sort()
#          searchEntry.subject_start = hspCoords[0] - 1 # Coords are 1-based
#          searchEntry.subject_length = hspCoords[-1] - hspCoords[0] + 1
# 
#          newQueryCoords = []
#          for x in _flatten(queryCoords):
#             newQueryCoords.append(x)
#          queryCoords = newQueryCoords
#          queryCoords.sort()
#          searchEntry.query_start = queryCoords[0] - 1 # Coords are 1-based
#          searchEntry.query_length = queryCoords[-1] - queryCoords[0] + 1
# 
#          if blastRecord.alignments[i].hsps[0].query_start > blastRecord.alignments[i].hsps[0].query_end:
#             searchEntry.query_strand = -1
#          else:
#             searchEntry.query_strand = 1

         searchEntry.score = 0
         hsp = blastRecord.alignments[i].hsps[0]
         searchEntry.score = hsp.bits

         if hsp.query_start > hsp.query_end:
            searchEntry.query_start = hsp.query_end - 1 # Coords are 1-based
            searchEntry.query_length = hsp.query_start - hsp.query_end + 1
            searchEntry.query_strand = -1
         else:
            searchEntry.query_start = hsp.query_start - 1 # Coords are 1-based
            searchEntry.query_length = hsp.query_end - hsp.query_start + 1
            searchEntry.query_strand = 1

         if hsp.sbjct_start > hsp.sbjct_end:
            searchEntry.subject_start = hsp.sbjct_end - 1 # Coords are 1-based
            searchEntry.subject_length = hsp.sbjct_start - hsp.sbjct_end + 1
            searchEntry.subject_strand = -1
         else:
            searchEntry.subject_start = hsp.sbjct_start - 1 # Coords are 1-based
            searchEntry.subject_length = hsp.sbjct_end - hsp.sbjct_start + 1
            searchEntry.subject_strand = 1




         searchEntry.significance = blastRecord.descriptions[i].e

         self.search_list.append(searchEntry)
