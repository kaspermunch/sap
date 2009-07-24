
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


class SearchResult(BlastSearchResult):   
   """
   Generic class for search attributes.
   """
   def __init__(self, blastRecord):
       BlastSearchResult.__init__(self, blastRecord)

       # Seperate the gi from the full hit titles:
       for i in range(len(self.search_list)):         
           self.search_list[i].id = self.search_list[i].id.split('|')[1]


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
           return SearchResult(blastRecord)


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


