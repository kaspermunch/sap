
try:
   import cPickle as pickle
except:
   import pickle
import os, sys, time, re, pickle

from SAP.Bio.EUtils.Datatypes import DBIds
from SAP.Bio.EUtils.ThinClient import ThinClient

from SAP import Fasta
from SAP import UtilityFunctions as utils
# from SAP import NCBIWWW # locally hacked to allow retrieval or more hits
from SAP import NCBIXML # locally hacked to better parse string info
from SAP import Taxonomy
from SAP import XML2Obj

from SAP.SearchResult import BlastSearchResult

from SAP import Homology # I should rename this to Homology

from SAP.ThreadedWorkerQueue import ThreadedWorkerQueue

# from SAP import NW

class SearchResult(BlastSearchResult):   
   """
   Generic class for search attributes.
   """
   def __init__(self, blastRecord):
       BlastSearchResult.__init__(self, blastRecord)

       # Seperate the gi from the full hit titles:
       for i in range(len(self.hits)):         
           self.hits[i].id = self.hits[i].id.split('|')[1]


class DB:

    def __init__(self, options):
        self.options = options
        self.prevExcludeList = None
        
    def search(self, fastaRecord, excludelist=[]):

        # Check that we get a new exclude list for every blast:
        if self.prevExcludeList is not None and excludelist == self.prevExcludeList:
            # Seems we are running in circles.
            return None
        else:
           self.prevExcludeList = excludelist
           
        # Get a blast record:
        useBlastCache = True
        for i in range(3):
            time.sleep(i * 2)
            #blastRecord = self.getBlastRecord(fastaRecord, excludelist, useBlastCache)
            blastRecord = self._netblast_search(fastaRecord, excludelist=excludelist, usecache=useBlastCache)
            if blastRecord is None:
                print "failed - trying again.\n\t\t\t",
                useBlastCache = False
            elif len(blastRecord.descriptions) == 0:
                print "no hits - strange - trying again.\n\t\t\t",
                useBlastCache = False
            else:
                break
        if blastRecord is None:
            print "WARNING: Blast failed - skipping entry\n\t\t\t",
            return None
        else:
           #return SearchResult(blastRecord)
           searchResult = SearchResult(blastRecord)

#           return searchResult

#            ###################################################################################################
#            # We can get the blast ranking form the list of homologues by sorting by evalue:
#            blastRanks = []
#            prevMax = searchResult.hits[i].score
#            rank = 1
#            for i in range(len(searchResult.hits)):
#                if searchResult.hits[i].score < prevMax:
#                    rank += 1
#                    prevMax = searchResult.hits[i].score
#                blastRanks.append((rank, searchResult.hits[i].id))
#            ###################################################################################################

           rerankedSearchResult = self._rerank_search_result(fastaRecord, searchResult)

#            ###################################################################################################
#            # We can get the blast ranking form the list of homologues by sorting by evalue:
#            rerankedRanks = []
#            prevMax = searchResult.hits[0].score
#            rank = 1
#            for i in range(len(searchResult.hits)):
#                if searchResult.hits[i].score < prevMax:
#                    rank += 1
#                    prevMax = searchResult.hits[i].score
#                rerankedRanks.append((rank, searchResult.hits[i].id))
# 
# 
#            distanceSum = 0
# 
#            for blastTup in blastRanks:
#               
#                found = False
#                for rerankedTup in rerankedRanks:
#                    if blastTup[1] == rerankedTup[1]:
#                        print blastTup[0], rerankedTup[0], blastTup[0] - rerankedTup[0]
#                        distanceSum += abs(blastTup[0]-rerankedTup[0])
#                        found = True
#                assert found, rerankedTup
# 
#            rerankedScore = float(distanceSum) / len(rerankedRanks)
# #            print len(rerankedRanks), rerankedScore
# # 
# #            sys.exit()
#            ###################################################################################################
           


           return rerankedSearchResult

    #fastaRecord, resultEntry, dbIdx
    def _rerank_search_result(self, fastaRecord, searchResult):
        """
        Retrieve all sequences and globally align them each to the query. Then rerank the
        list of hits based on alignment scores.
        """

        print 're-ranking results...',
        sys.stdout.flush()

        # get list of all hit ids:
        idList = [x.id for x in  searchResult.hits]

        # set up download stream:
        eutils = ThinClient()
        dbids = DBIds("nucleotide",  idList)
        fileob = eutils.efetch_using_dbids(dbids, retmode="text", rettype="fasta")    

        # make an iterator to read fasta entries as we get them:
        fastaIterator = Fasta.Iterator(fileob, Fasta.RecordParser())

        # set up a threaded queue for writing files:
        fileWriteQueue = ThreadedWorkerQueue(utils.writeFile)

        for resultIdx, resultEntry in enumerate(searchResult.hits):

            # Get next fasta record from the iterator:
            hitFastaEntry = fastaIterator.next()

            # Convert title to just the gi nr:
            hitFastaEntry.title = hitFastaEntry.title.split('|')[1]

            # make sure we are in sync:
            assert hitFastaEntry.title == resultEntry.id , hitFastaEntry.title + " " +  resultEntry.id

            # Queue the hit fasta entry for writing to the cache:
            fastaFileName = os.path.join(self.options.dbcache, hitFastaEntry.title + ".fasta")
            fileWriteQueue.enqueue([fastaFileName, str(hitFastaEntry)])

            # align and remap coordinates:
            resultEntry = utils.remap_db_hit(fastaRecord, hitFastaEntry, resultEntry)

            # populate it with the sequence:
            resultEntry.sequence = hitFastaEntry.sequence[resultEntry.subject_start:resultEntry.subject_start + resultEntry.subject_length]

            # updata search list:
            searchResult.hits[resultIdx] = resultEntry
       
        # sort list of hits based on new scores:
        searchResult.hits.sort(lambda a, b: cmp(b.score,a.score))

        # Wait for the remaining files to be written:
        fileWriteQueue.close()
        
        print "done.\n\t\t\t",
        sys.stdout.flush()

        return searchResult


    def _printAlignment(self, seq1, seq2):

        assert len(seq1) == len(seq2)

        width = 100
        end = width
        alnLength = len(seq1)
        print
        for start in range(0, alnLength, width):
            print seq1[start:min(end, alnLength)]
            print seq2[start:min(end, alnLength)]
            end += width
            print
        print

    def _netblast_search(self, fastaRecord, excludelist=[], usecache=True):
        """
        Blast against genbank over web
        """

        # Make a query to filter the returned results:
        if excludelist:
            entrezQuery = '(' + self.options.limitquery + ') NOT (uncultured[WORD] OR ' + '[ORGN] OR '.join(excludelist) + '[ORGN])'
        else:
            entrezQuery = '(' + self.options.limitquery + ') NOT uncultured[WORD]'

        fileSuffix = ''
        for name in excludelist:
            l = re.split(r'\s+', name)
            for n in l:
                fileSuffix += n[0]
        if fileSuffix:
           fileSuffix = '_' + fileSuffix

        # File name used for blast cache
        blastFileName = os.path.join(self.options.blastcache, "%s.%d_%s%s.xml" % (fastaRecord.title,
                                               self.options.maxblasthits, self.options.minsignificance, fileSuffix))
        
        if usecache and os.path.exists(blastFileName) and os.path.getsize(blastFileName)>0:
            # Use cached blast result
            if excludelist:
               print "\n\t\tUsing cached Blast results (excluding %s)..." % ', '.join(excludelist),
            else:
                print "\n\t\tUsing cached Blast results...", 
            sys.stdout.flush()
        else:
            # Make a query to filter the returned results:
            if excludelist:
                print "\n\t\tSearching database (excluding %s)..." % ', '.join(excludelist), 
            else:
                print "\n\t\tSearching database...", 
            sys.stdout.flush()

            fastaRecordFileName = os.path.join(self.options.project, utils.randomString(8))
            fastaRecordFile = open(fastaRecordFileName, 'w')
            fastaRecordFile.write(str(fastaRecord))
            fastaRecordFile.close()
            resultHandle = None
            if self.options.nolowcomplexfilter:
                filterOption = '-F F'
            else:
                filterOption = '-F T'

            if self.options.blastwordsize:
                wordSize = '-W %s' % self.options.blastwordsize
            else:
                wordSize = ''

            blastCmd = 'blastcl3 -p blastn -m 7 -d nr %s %s -e %s -v %d -b %d -u "%s" -i %s -o %s' \
                       % (wordSize, filterOption, self.options.minsignificance, self.options.maxblasthits, self.options.maxblasthits, \
                          entrezQuery, fastaRecordFileName, blastFileName)

            for i in range(20):
                time.sleep(2 * i)
                try:
                    #os.system(blastCmd)
                    retval = os.system(blastCmd)
                    #retval = utils.systemCall(blastCmd, stdout='IGNORE', stderr='IGNORE')
                    sys.stdout.flush()
                    if retval != 0:
                       print "Netblast failed with return value %d. Trying again..." % retval
                       continue
                    break
                except:
                    print "Netblast failed. Trying again..."
                    pass
            os.remove(fastaRecordFileName)

        # Read file from cache
        blastHandle = open(blastFileName, 'r')

        # Parse the result:
        blastParser = NCBIXML.BlastParser()
        try:
            blastRecord = blastParser.parse(blastHandle)
            print "done.\n\t\t\t",
            sys.stdout.flush()
        except:
            blastRecord = None
        blastHandle.close()

        return blastRecord
     
    # SEE IF YOU CAN AVOID DOWNLOADING THE FULL GENBANK XML SOMEHOW IF YOU HAVE THE
    # FASTA. IS THERE SOME UTILITY IN THE EUTILS FOR THE TAXONOMY DATABASE THAT LETS
    # ME GET THE TAXID FROM THE GI WITHOUT THE FULL XML?

#     def getBatchTaxonomies(self, giList):

       


  
    def get(self, gi):
        """
        Look up genbank records by their GI
        """

        taxonomyFileName = os.path.join(self.options.dbcache, gi + ".tax")
        fastaFileName = os.path.join(self.options.dbcache, gi + ".fasta")

        if (os.path.exists(taxonomyFileName) and os.path.getsize(taxonomyFileName) != 0 and
            os.path.exists(fastaFileName) and os.path.getsize(fastaFileName) != 0):
            retrievalStatus = "(c)"
            taxonomy = utils.safeReadTaxonomyCache(taxonomyFileName)
            sequence = utils.safeReadFastaCache(fastaFileName)
        else:
            retrievalStatus = "(d)"

            taxonXref = None
            seqLength = None

            successful = False
            for tries in range(10):
                try:
                    eutils = ThinClient()
                    dbids = DBIds("nucleotide",  [str(gi)])
                    fp = eutils.efetch_using_dbids(dbids, retmode="xml")

                    # Get the cross ref to the taxonomy database:
                    taxonXrefRE = re.compile("<GBQualifier_value>taxon:(\d+)</GBQualifier_value>")
                    seqLengthRE = re.compile("<GBSeq_length>(\d+)</GBSeq_length>")
                    sequenceRE = re.compile("<GBSeq_sequence>([a-zA-Z]+)</GBSeq_sequence>")

                    taxonXref = None
                    seqLength = None
                    sequence = None

                    while taxonXref is None or sequence is None:
                        line = fp.readline()
                        if not line:
                            break
                        taxonMatch = taxonXrefRE.search(line)
                        lengthMatch = seqLengthRE.search(line)
                        sequenceMatch = sequenceRE.search(line)
                        if taxonMatch:
                            if taxonXref is None:
                               taxonXref = taxonMatch.group(1)
#                             else:
#                                print "There was more than one taxon xref for %s. Picking the first one (%s)." % (gi, taxonXref)
                        if lengthMatch:
                            seqLength = lengthMatch.group(1)
                        if sequenceMatch:
                            sequence = sequenceMatch.group(1)

                    if not (taxonXref and sequence):
                       # Give it another try:
                       continue

                except KeyboardInterrupt:
                   sys.exit()
                except MemoryError:
                    # Write an empty file to cache to keep the script from
                    # trying to download the sequence next time.
                    utils.writeFile(fastaFileName, '')
                    return None, retrievalStatus.replace(")", "!M)")
                except:
                   ## print ' retrieving failed - retrying'
                   time.sleep(tries * 5)
                   continue
                else:
                   successful = True
                   fp.close()
                   break
                if not successful:
                   return None, retrievalStatus.replace(")", "!D2)")

            if not (taxonXref and gi and sequence):
                # The entry did not have a cross ref to the taxonomy database:
                return None, retrievalStatus.replace(")", "!T2)")

            # Make an object to hold the taxonomy:
            taxonomy = Taxonomy.Taxonomy()
            try:
               taxonomy.populateFromNCBI(dbid=taxonXref,
#                                          allow_unclassified=self.options.unclassified,
                                         minimaltaxonomy=self.options.minimaltaxonomy)
            except Taxonomy.NCBIPopulationError, X:
               return None, retrievalStatus.replace(")", " !%s)" % X.status)
               
            # Dump the taxonomy object to a file:
            fp = open(taxonomyFileName, 'w')
            pickle.dump(taxonomy, fp)
            fp.close()

            # Upcase the sequence:
            sequence = sequence.upper()

            # Cache the sequence:
            fastaEntry = ">%s\n%s\n" % (gi, sequence)
            utils.writeFile(fastaFileName, fastaEntry)

        # Create an object instance to hold the homologue data:
        return Homology.Homologue(gi=gi,
                                  sequence=sequence,
                                  taxonomy=taxonomy), retrievalStatus

#     def _safeReadFastaCache(self, fastaFileName):
#        """
#        Reads in a genbank fasta entry, and tries again if reading
#        failes because others are writing to it.
#        """
#        fastaEntry = None
#        for i in range(5):
#           fastaFile = open(fastaFileName, 'r')
#           fastaIterator = Fasta.Iterator(fastaFile)
#           try:
#              fastaEntry = fastaIterator.next()
#           except:             
#              print 'Fasta cache reading failed - retrying'
#           fastaFile.close()
#           if fastaEntry is not None:
#              break
#           time.sleep(i * 2)
#        assert fastaEntry, fastaFileName
#        return fastaEntry.sequence
# 
# #     def _safeReadFastaCache(self, fastaFileName):
# #        """
# #        Reads in a genbank fasta entry, and tries again if reading
# #        failes because others are writing to it.
# #        """
# #        fastaEntry = None
# #        for i in range(5):
# #           fastaFile = open(fastaFileName, 'r')
# #           fastaIterator = Fasta.Iterator(fastaFile)
# #           try:
# #              fastaEntry = fastaIterator.next()
# #           except:             
# #              print 'Fasta cache reading failed - retrying'
# #           fastaFile.close()
# #           if fastaEntry is not None:
# #              break
# #           time.sleep(i * 2)
# #        assert fastaEntry, fastaFileName
# #        seq = fastaEntry.seq.data
# #        return seq
# 
#     def _safeReadTaxonomyCache(self, taxonomyFileName):
#        """
#        Reads in a genbank taxonomy entry, and tries again if reading
#        failes because others are writing to it.
#        """
#        taxonomy = None
#        for i in range(5):
#           taxonomyFile = open(taxonomyFileName, 'r')
#           try:
#              taxonomy = pickle.load(taxonomyFile)
#           except:
#              print 'Taxonomy cache reading failed - retrying'
#           taxonomyFile.close()    
#           if taxonomy is not None:
#              break
#           time.sleep(i * 2)
#        assert taxonomy is not None, taxonomyFile
#        return taxonomy


if __name__ == "__main__":
    main()


