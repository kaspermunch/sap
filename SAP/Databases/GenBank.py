
try:
   import cPickle as pickle
except:
   import pickle
import os, sys, time, re, pickle

from SAP.Bio.EUtils.Datatypes import DBIds
from SAP.Bio.EUtils.ThinClient import ThinClient

from SAP import Fasta
from SAP import UtilityFunctions as util
# from SAP import NCBIWWW # locally hacked to allow retrieval or more hits
from SAP import NCBIXML # locally hacked to better parse string info
from SAP import Taxonomy
from SAP import XML2Obj

from SAP import Homology # I should rename this to Homology

class DB:

    def __init__(self, options):
        self.options = options

    def search(self, fastaRecord, excludelist=[], usecache=True):
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
                                               self.options.maxblasthits, self.options.evaluecutoff, fileSuffix))
        
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

            fastaRecordFileName = os.path.join(self.options.project, util.randomString(8))
            fastaRecordFile = open(fastaRecordFileName, 'w')
            fastaRecordFile.write(str(fastaRecord))
            fastaRecordFile.close()
            resultHandle = None
            if self.options.nolowcomplexfilter:
                filterOption = '-F F'
            else:
                filterOption = '-F T'
            blastCmd = 'blastcl3 -p blastn -m 7 -d nr %s -e %s -v %d -b %d -u "%s" -i %s -o %s' \
                       % (filterOption, self.options.evaluecutoff, self.options.maxblasthits, self.options.maxblasthits, \
                          entrezQuery, fastaRecordFileName, blastFileName)

            for i in range(20):
                time.sleep(2 * i)
                try:
                    os.system(blastCmd)
                    break
                except:
                    print "Blast failed. Trying again..."
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


  
    def get(self, gi, evalue):
        """
        Look up genbank records by their GI
        """

        taxonomyFileName = os.path.join(self.options.genbankcache, gi + ".tax")
        fastaFileName = os.path.join(self.options.genbankcache, gi + ".fasta")

        if (os.path.exists(taxonomyFileName) and os.path.getsize(taxonomyFileName) != 0 and
            os.path.exists(fastaFileName) and os.path.getsize(fastaFileName) != 0):
            retrievalStatus = "(c)"
            taxonomy = self.safeReadTaxonomyCache(taxonomyFileName)
            sequence = self.safeReadFastaCache(fastaFileName)
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
                            taxonXref = taxonMatch.group(1)
                        if lengthMatch:
                            seqLength = lengthMatch.group(1)
                        if sequenceMatch:
                            sequence = sequenceMatch.group(1)
                except MemoryError:
                    # Write an empty file to cache to keep the script from
                    # trying to download the sequence next time.
                    util.writeFile(fastaFileName, '')
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
                                         minimaltaxonomy=self.options.minimaltaxonomy,
                                         subspecieslevel=self.options.subspecieslevel)
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
            util.writeFile(fastaFileName, fastaEntry)

        # Create an object instance to hold the homologue data:
        return Homology.Homologue(gi=gi,
                                  sequence=sequence,
                                  evalue = evalue, 
                                  taxonomy=taxonomy,
                                  options=self.options), retrievalStatus

    def safeReadFastaCache(self, fastaFileName):
       """
       Reads in a genbank fasta entry, and tries again if reading
       failes because others are writing to it.
       """
       fastaEntry = None
       for i in range(5):
          fastaFile = open(fastaFileName, 'r')
          fastaIterator = Fasta.Iterator(fastaFile)
          try:
             fastaEntry = fastaIterator.next()
          except:             
             print 'Fasta cache reading failed - retrying'
          fastaFile.close()
          if fastaEntry is not None:
             break
          time.sleep(i * 2)
       assert fastaEntry, fastaFileName
       return fastaEntry.sequence

#     def safeReadFastaCache(self, fastaFileName):
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
#        seq = fastaEntry.seq.data
#        return seq

    def safeReadTaxonomyCache(self, taxonomyFileName):
       """
       Reads in a genbank taxonomy entry, and tries again if reading
       failes because others are writing to it.
       """
       taxonomy = None
       for i in range(5):
          taxonomyFile = open(taxonomyFileName, 'r')
          try:
             taxonomy = pickle.load(taxonomyFile)
          except:
             print 'Taxonomy cache reading failed - retrying'
          taxonomyFile.close()    
          if taxonomy is not None:
             break
          time.sleep(i * 2)
       assert taxonomy is not None, taxonomyFile
       return taxonomy

if __name__ == "__main__":
    main()


