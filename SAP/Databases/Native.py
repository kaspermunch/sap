#!/usr/bin/python

import sys, os, re, shelve, tempfile, StringIO
from SAP import Fasta, Taxonomy, Homology
from SAP.Bio.Alphabet import IUPAC
from SAP.FastaIndex import FastaIndex
from SAP.InstallDependencies import assertBlastInstalled
from SAP.UtilityFunctions import *
from SAP.FindPlugins import *
from SAP.Bio.Blast import NCBIXML
try:
    import cPickle as pickle
except ImportError:
    import pickle

from SAP.SearchResult import BlastSearchResult

from SAP.Exceptions import AnalysisTerminated


class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass

class WrongFormatError(Error):
    """
    Exception raised for errors in alignment of taxonomic levels.
    """
    def __init__(self, msg=None):
        self.msg = "Wrong format of FASTA header line: " + msg

class SearchResult(BlastSearchResult):   
    """
    Generic class for search attributes.
    """
    pass

class DB(object):

    def __init__(self, fastaFileName, options, rebuild=False):
        """
        Opens a premade local data base and builds an index of taxon names.
        """

        self.options = options
        self.fastaFileName = fastaFileName

        baseName, ext = os.path.splitext(os.path.basename(self.fastaFileName))
        self.baseName = baseName
        self.dbDirName = os.path.splitext(self.fastaFileName)[0] + '.sapdb'
        self.blastSequenceFileName = os.path.join(self.dbDirName, baseName + '.fasta')
        self.blastDB = os.path.join(self.dbDirName, baseName)
        #self.dbFileName = os.path.join(self.dbDirName, baseName + '.db')
        self.indexFileName = os.path.join(self.dbDirName, baseName + '.taxidx')
        self.summaryTreeFileName = os.path.join(self.dbDirName, baseName + '.html')

        # Build shelve database if not present:
        if rebuild or not os.path.exists(self.dbDirName) or len(glob.glob(os.path.join("%s/*" % self.dbDirName))) != 6:
            self._build()

        print "Loading database...",
        sys.stdout.flush()
        # Just open the db:
        #self.db = shelve.open(self.dbFileName, 'r')
        with open(self.indexFileName, 'r') as f:
            self.index = pickle.load(f)

        self.sequence_index = FastaIndex(self.blastSequenceFileName)
        print 'done'
        sys.stdout.flush()

    @staticmethod
    def escape(s):
        return s.replace('(', '\(').replace(')', '\)')


    def _build(self):
        """
        Build a local database from a file of fasta entries.
        """

        print "Building native SAP database...",
        sys.stdout.flush()

        # Make a dir for the database files:
        if not os.path.exists(self.dbDirName):
            os.mkdir(self.dbDirName)

        # Write a newline char version to a tmp file:
        fileContent = readFile(self.fastaFileName)
        fileContent = re.sub(r'\r+', '\n', fileContent)
        tmpFile, tmpFileName = tempfile.mkstemp()
        writeFile(tmpFileName, fileContent)
        self.fastaFileName = tmpFileName

        # # Create a new shelve:
        # db = shelve.open(self.dbFileName, 'n')

        # Create an iterator over fasta entries:
        fastaFile = open(self.fastaFileName, 'r')
        fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        

        allowedLetters = IUPAC.IUPACAmbiguousDNA().letters

        # Initialize the index:
        index = {}

        # File to write sequences without gaps to for use as blast database:
        blastSequenceFile = open(self.blastSequenceFileName, 'w')

        taxonomySummary = Taxonomy.TaxonomySummary()

        for fastaRecord in fastaIterator:

            taxonomy = Taxonomy.Taxonomy()
            try:
                identifier, taxonomyString, organismName = re.split(r'\s*;\s*', fastaRecord.title.strip())
                taxonomy.populateFromString(taxonomyString)
                taxonomy.organism = organismName
            except (ValueError, Taxonomy.ParsingError) as exc:
                raise WrongFormatError(fastaRecord.title)
            if not (identifier and taxonomyString and organismName):
                raise WrongFormatError(fastaRecord.title)

            # Add an entry to the index:
            for taxonomyLevel in taxonomy:
                index.setdefault(taxonomyLevel.name, []).append(identifier)

            # taxonomySummary.addTaxonomy(taxonomy)

            # Uppercase and replace wildcard letters with Ns:
            fastaRecord.sequence = fastaRecord.sequence.upper()
            fastaRecord.sequence = re.sub('[^%s-]' % allowedLetters, 'N', fastaRecord.sequence)

            # Write the sequence without gaps to the blast file:
            fastaRecord.sequence = fastaRecord.sequence.replace('-', '')
            blastSequenceFile.write(str(fastaRecord))

        # remove tmp files:
        os.close(tmpFile)
        os.unlink(tmpFileName)

        # Close the fasta file:
        fastaFile.close()

        # Pickle the index:
        indexFile = open(self.indexFileName, 'w')
        pickle.dump(index, indexFile)
        indexFile.close()

# <h1>Taxonomic summary for local database: ''' + self.dbDirName + '</h1><div class="taxonomysummary"><pre>' + str(taxonomySummary) + "</pre></div></html>"
#
#         fh = open(self.summaryTreeFileName, 'w')
#         fh.write(html)
#         fh.close()

        # Format the blast database:
        blastSequenceFile.close()
        cmd = "makeblastdb -dbtype nucl -title %s -in %s" % (self.blastDB, self.blastSequenceFileName)

        cmd = self.escape(cmd)

        systemCall(cmd, stdout='IGNORE', stderr='IGNORE')

        # Create an index for the blastSequenceFile and pickle it:
        self.sequence_index = FastaIndex(self.blastSequenceFileName)

        # Make sure all the files have been written:
        import glob, time
        tries = 5
        while tries and len(glob.glob(os.path.join("%s/*" % self.dbDirName))) != 7:
            time.sleep(1)
            tries -= 1

        print "done"
        sys.stdout.flush()

    def search(self, fastaRecord, excludelist=[], usecache=True):

        # Write the query to a tmp file:
        tmpQueryFileIdent, tmpQueryFileName = tempfile.mkstemp()
        writeFile(tmpQueryFileName, str(fastaRecord))

        # File name used for blast cache
        fileSuffix = ''
        for name in excludelist:
            l = re.split(r'\s+', name)
            for n in l:
                fileSuffix += n[0]
        if fileSuffix:
           fileSuffix = '_' + fileSuffix

        blastFileName = os.path.join(self.options.blastcache, "%s.%d_%s%s.xml" % (fastaRecord.title,
                                               self.options.maxblasthits, self.options.minsignificance, fileSuffix))

        if usecache and os.path.exists(blastFileName) and os.path.getsize(blastFileName)>0:
            # Use cached blast result
            if excludelist:
               print "\n\t\tUsing cached Blast results (excluding %s)..." % ', '.join(excludelist),
            else:
                print "\n\t\tUsing cached Blast results...", 
            sys.stdout.flush()

            blastFile = open(blastFileName, 'r')
        else:
            if excludelist:
                print "\n\t\tSearching database (excluding %s)..." % ', '.join(excludelist), 
            else:
                print "\n\t\tSearching database...", 
            sys.stdout.flush()

            if excludelist:
                # Build a blast db of the relevant records:
                excludeIDs = []
                for taxon in excludelist:
                    excludeIDs.extend(self.index[taxon])
                includeIDs = self.sequence_index.keys.difference(set(excludeIDs))

                if not includeIDs:
                    print "done (database exhausted)\n\t\t\t",
                    sys.stdout.flush()       
                    return SearchResult(None)

                blastDBfileName = os.path.join(self.options.project, "tmpBlastDB.fasta")
                tmpFastaFile = open(blastDBfileName, 'w')

                for key in includeIDs:
                    tmpFastaFile.write(str(self.sequence_index.get_entry(key)))
                tmpFastaFile.close()
                cmd = "makeblastdb -dbtype nucl -title %s -in %s" % (blastDBfileName, blastDBfileName)
                #cmd = "xdformat -n -o %s %s" % (blastDBfileName, blastDBfileName)
                cmd = self.escape(cmd)
                systemCall(cmd, stdout='IGNORE', stderr='IGNORE')
            else:
                blastDBfileName = self.blastSequenceFileName
                #blastDBfileName = self.blastDB
            
            # Blast:
            if self.options.blastwordsize:
                wordSize = '-word_size %s' % self.options.blastwordsize
            else:
                wordSize = ''

            if self.options.nolowcomplexfilter:
                filterOption = '-dust no'
            else:
                filterOption = '-dust yes' # FIXME: Check that this is an ok default... It is not the defalut in blastn

            cmd = 'blastn -db %s -outfmt 5 %s %s -evalue %s -max_target_seqs %s -query %s' \
                       % (blastDBfileName, wordSize, filterOption, self.options.minsignificance, self.options.maxblasthits,
                          tmpQueryFileName)
               
#            cmd = "blastn %s -e %f -p blastn -d %s.fasta -i %s -m 7" % (wordSize, self.options.minsignificance, blastDBfileName, tmpQueryFileName)
            cmd = self.escape(cmd)

            STARTUPINFO = None
            if os.name == 'nt':
                STARTUPINFO = subprocess.STARTUPINFO()
                STARTUPINFO.dwFlags |= subprocess.STARTF_USESHOWWINDOW
            proc = subprocess.Popen(cmd, shell=True, startupinfo=STARTUPINFO,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout_value, stderr_value = proc.communicate()
            blastContent = str(stdout_value)

#            stdin, stdout, stderr = os.popen3(cmd)

            for f in glob.glob(os.path.join(self.options.project, 'tmpBlastDB.*')):
                os.remove(f)

#             blastContent = stdout.read()
# 
#             stdout.close()
#             stdin.close()
#             stderr.close()

            # This is a hack to remove an xml tag that just confuses things with blastn:
            blastContent = re.sub(r'<Hit_id>.*?</Hit_id>', '', blastContent)

            writeFile(blastFileName, blastContent)

            #blastFile = StringIO.StringIO(blastContent)
            blastFile = open(blastFileName, 'r')

        #blastRecord = NCBIXML.read(blastFile)
        try:
            blastRecord = NCBIXML.read(blastFile)
        except:
            blastRecord = None

        blastFile.close()

        os.close(tmpQueryFileIdent)
        os.unlink(tmpQueryFileName)

        print "done.\n\t\t\t",
        sys.stdout.flush()
        
        return SearchResult(blastRecord)
        
    def get(self, gi):

        try:
            record = self.sequence_index.get_entry(gi)
            taxonomy = Taxonomy.Taxonomy()
            try:
                identifier, taxonomyString, organismName = re.split(r'\s*;\s*', record.title)
                taxonomy.populateFromString(taxonomyString)
                taxonomy.organism = organismName
            except (ValueError, Taxonomy.ParsingError) as exc:
                raise WrongFormatError(record.title)
            if not (identifier and taxonomyString and organismName):
                raise WrongFormatError(record.title)

            retrievalStatus = '(l)'
        except:
            retrievalStatus = '(l!?)'
            return None, retrievalStatus

        return Homology.Homologue(gi=gi,
                                  sequence=record.sequence,
                                  taxonomy=taxonomy), retrievalStatus


if __name__ == "__main__":

    pass
    # import sys
    # inputfile = sys.argv[1]

    # db = DB(inputfile, rebuild=True)
    # #db = DB(inputfile)

    # print db.sequence_index.keys()
    # print
    # print db.index.keys()

    # fastaRecord = Fasta.Record(title='myfasta',
    #                            sequence='TTTAATGTTAGGAGCTCCAGACATGGCATTCCCTCGTATAAATAATATAAGATTTTGATTATTACCCCCTAAGAAGAATAGTTGAAAGAGGAGC')
    #
    # #blastRecord = db.search(fastaRecord, excludelist=['Tachinus brevipennis'])
    # blastRecord = db.search(fastaRecord)
    #
    # for i in range(len(blastRecord.alignments)):
    #
    #     gi = blastRecord.alignments[i].title.split(' ')[0]
    #
    #     description = blastRecord.descriptions[i]
    #
    #     print db.get(gi, description.e)

