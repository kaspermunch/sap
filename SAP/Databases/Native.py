#!/usr/bin/python

import sys, os, re, shelve, tempfile, StringIO
from sets import Set
from SAP import Fasta, Taxonomy, Homology
from SAP.Bio.Alphabet import IUPAC
from SAP.InstallDependencies import assertBlastInstalled
from SAP.UtilityFunctions import *
from SAP.FindPlugins import *
from SAP import NCBIXML # locally hacked to better parse string info
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

        self.got_wublast = False
        #self.got_wublast = findOnSystem('blastn') and findOnSystem('xdformat')

        self.options = options
        self.fastaFileName = fastaFileName

        baseName, ext = os.path.splitext(os.path.basename(self.fastaFileName))
        self.baseName = baseName
        self.dbDirName = os.path.splitext(self.fastaFileName)[0] + '.sapdb'
        self.blastSequenceFileName = os.path.join(self.dbDirName, baseName + '.fasta')
        self.blastDB = os.path.join(self.dbDirName, baseName)
        self.dbFileName = os.path.join(self.dbDirName, baseName + '.db')
        self.indexFileName = os.path.join(self.dbDirName, baseName + '.idx')
        self.summaryTreeFileName = os.path.join(self.dbDirName, baseName + '.html')

        # Build shelve database if not present:
        if rebuild or not os.path.exists(self.dbDirName) or len(glob.glob(os.path.join("%s/*" % self.dbDirName))) != 7:
            self._build()

        # Just open the db:
        self.db = shelve.open(self.dbFileName, 'r')
        indexFile = open(self.indexFileName, 'r')
        self.index = pickle.load(indexFile)


    def escape(self, s):
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

        # Create a new shelve:
        db = shelve.open(self.dbFileName, 'n')

        # Create an iterator over fasta entries:
        fastaFile = open(self.fastaFileName, 'r')
        fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        

        allowedLetters = IUPAC.IUPACAmbiguousDNA().letters
            
        # Initialize the index:
        index = {}

        # Length of alignment:
        alignmentLength = None

        # Boolean indicating if sequences seemm to be pre-aligned:
        isPreAligned = True

        # Attribute to hold GenBank plugin if necessary:
        self.genbank = None

        # File to write sequences without gaps to for use as blast database:
        blastSequenceFile = open(self.blastSequenceFileName, 'w')

        taxonomySummary = Taxonomy.TaxonomySummary()

        firstMissedTaxonomy = True
        for fastaRecord in fastaIterator:
            try:
                # See if we can get the taxonomy form the fasta title line as we are supposed to:
                taxonomy = Taxonomy.Taxonomy()                        
                try:
                    identifier, taxonomyString, organismName = re.split(r'\s*;\s*', fastaRecord.title)
                except ValueError:
                    raise WrongFormatError(fastaRecord.title)
                identifier.strip()
                taxonomyString.strip()
                organismName.strip()
                if not (identifier and taxonomyString and organismName):
                    raise WrongFormatError(fastaRecord.title)
                fastaRecord.title = identifier        
                taxonomy.populateFromString(taxonomyString)
                taxonomy.organism = organismName

            except WrongFormatError, exe:

                #print exe.msg

                if firstMissedTaxonomy:
                    print "\nRetrieving missing taxonomic annotation from NCBI"
                    sys.stdout.flush()
                    firstMissedTaxonomy = False

                # if not - load the genbank plugin and retrieve from genbank
                if self.genbank is None:
                    try:
                        plugin = findPlugin('GenBank', 'sap.database')
                        self.genbank = plugin.DB(self.options)
                    except PluginNotFoundError, X:
                        # If that fails it is either because we are running an osx application in which
                        # case we can not load dynamically. So we try to see if the plugin should be
                        # one of the ones that come with the distribution in which case we can just
                        # import it:
                        from SAP.Databases import GenBank as plugin
                        self.genbank = plugin.DB(self.options)

                # If the fasta file was downloaded form GenBank the header lines will look like this: 'gi|187481301|gb|EU154882.1| Turdus pilar...
                genBankGiMatch = re.match(r'gi\|(\d+)\|', fastaRecord.title)
                idMatch = re.match(r'([\d\w.]+)', fastaRecord.title)
                if genBankGiMatch:
                    dbid = genBankGiMatch.group(1)
                elif idMatch:
                    dbid = idMatch.group(1)
                else:
                    raise AnalysisTerminated(1, 'Taxonomic annotation in sequence header missing or wrongly formatted and sequence id is not a NCBI gi or accession number either.\nSee the manual pages for how to format a local dabase file.')              

                homologue, retrieval_status = self.genbank.get(dbid)

                fastaRecord.title = safeName(fastaRecord.title)
                # we remove the underscores because we need to be able to split on '-' to get the indentifier later...
                fastaRecord.title = fastaRecord.title.replace('_', '') 

                if homologue is None:
                    print "WARNING: Retrieval annotation for \"%s\" was not possible (error code %s)." % (dbid, retrieval_status)
                    continue
                    
                taxonomy = homologue.taxonomy

            taxonomySummary.addTaxonomy(taxonomy)

            # Make sure we don't have doubles:
            if db.has_key(fastaRecord.title):
                raise AnalysisTerminated(1, "WARNING: Database sequences must have unique names. Double found: %s" % fastaRecord.title)

            # Replace wildcard letters with Ns:
            fastaRecord.sequence = re.sub('[^%s-]' % allowedLetters, 'N', fastaRecord.sequence)

            # Create a database entry:
            db[str(fastaRecord.title)] = {'fastaRecord': fastaRecord,
                                   'taxonomy': taxonomy}

            # Write the sequence without gaps to the blast file:
            fastaRecord.sequence = fastaRecord.sequence.replace('-', '')
            blastSequenceFile.write(str(fastaRecord))

            ## # Print to an extra fasta file with the formatted header included:
            ## fastaRecord.title += " ; %s ; %s" (str(taxonomy), taxonomy.organism)
            ## fastaSequenceFile.write(str(fastaRecord))

            # Add an entry to the index:
            for taxonomyLevel in taxonomy:
                index.setdefault(taxonomyLevel.name, []).append(fastaRecord.title)

        # Syncronize shelve:
        db.sync()
        db.close()

        # remove tmp files:
        os.close(tmpFile)
        os.unlink(tmpFileName)

        # Close the fasta file:
        fastaFile.close()

        # Pickle the index:
        indexFile = open(self.indexFileName, 'w')
        pickle.dump(index, indexFile)
        indexFile.close()
        
        # Write the summary tree:
        html = '''<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<style type="text/css">

body {  
    font-family: verdana, geneva, arial, helvetica, sans-serif;
    font-size: 12px; 
    font-style: normal; 
    text-decoration: none; 
    font-weight: lighter; 
}

a {
    font-style: normal; 
    font-weight: lighter; 
}

.taxonomysummary {
    background-color: #ffffff;
    padding-left: 40px;
    padding-top: 15px;
}

.header {
    background-color: #999999;
}

.dark {
    background-color: #dddddd;
}

.light {
    background-color: #eeeeee;
}

p { 
    font-size: 12px;
    margin-top: 2px;
    margin-left: 5px;
}

h1 {
    font-size: 13pt;
    line-height: 13pt;
    background-color: #ddd;
    /* border:2px solid black; */
    padding-bottom: 0pt;
    margin-bottom: 0pt;
/*     margin-left: 3pt; */
    padding: 5pt;
}

h2 {
    font-size: 11pt;
    margin-bottom: 0pt;
    margin-left: 3pt;
    font-weight: bold
}

a:active blue {
    color: #666666;
}

a:active {
    color: #666666;
}

a:link {
    color: #000000;
    text-decoration: underline;
}

a:visited {
    color: #000000;
    text-decoration: underline;
}

a:hover {
    color: #000000;
    text-decoration: underline;
}
</style>
<link rel="stylesheet" type="text/css" href="tooltip.css" >
<title>Local database summary</title>
</head>
<h1>Taxonomic summary for local database: ''' + self.dbDirName + '</h1><div class="taxonomysummary"><pre>' + str(taxonomySummary) + "</pre></div></html>"

        fh = open(self.summaryTreeFileName, 'w')
        fh.write(html)
        fh.close()

        # Format the blast database:
        blastSequenceFile.close()
        if self.got_wublast:
            cmd = "xdformat -n -o %s %s" % (self.blastDB, self.blastSequenceFileName)
        else:
            cmd = "formatdb -p F -t %s -i %s" % (self.blastDB, self.blastSequenceFileName)

        cmd = self.escape(cmd)

        systemCall(cmd, stdout='IGNORE', stderr='IGNORE')
                
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
                includeIDs = Set(self.db.keys()).difference(excludeIDs)

                blastDBfileName = "tmpBlastDB.fasta"
                tmpFastaFile = open(blastDBfileName, 'w')
                for key in includeIDs:
                    tmpFastaFile.write(str(self.db[str(key)]['fastaRecord']))
                tmpFastaFile.close()
                cmd = "xdformat -n -o %s %s" % (blastDBfileName, blastDBfileName)
                cmd = self.escape(cmd)
                systemCall(cmd, stdout='IGNORE', stderr='IGNORE')
            else:
                blastDBfileName = self.blastDB
            
            # Blast:
            if self.options.blastwordsize:
                wordSize = '-W %s' % self.options.blastwordsize
            else:
                wordSize = ''
                
            if self.got_wublast:
                cmd = "blastn %s -E %s %s -mformat 7" % (wordSize, self.options.minsignificance, blastDBfileName, tmpQueryFileName)
            else:
                cmd = "blastall %s -e %f -p blastn -d %s.fasta -i %s -m 7" % (wordSize, self.options.minsignificance, blastDBfileName, tmpQueryFileName)
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

            for f in glob.glob('tmpBlastDB.*'):
                os.remove(f)

#             blastContent = stdout.read()
# 
#             stdout.close()
#             stdin.close()
#             stderr.close()

            # This is a hack to remove an xml tag that just confuses things with blastn:
            if not self.got_wublast:
                blastContent = re.sub(r'<Hit_id>.*?</Hit_id>', '', blastContent)

            writeFile(blastFileName, blastContent)

            blastFile = StringIO.StringIO(blastContent)

        blastParser = NCBIXML.BlastParser()
        try:
            blastRecord = blastParser.parse(blastFile) 
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
            fastaRecord = self.db[str(gi)]['fastaRecord']
            taxonomy = self.db[str(gi)]['taxonomy']
            retrievalStatus = '(l)'
        except:
            retrievalStatus = '(l!?)'
            return None, retrievalStatus

        return Homology.Homologue(gi=gi,
                                  sequence=fastaRecord.sequence,
                                  taxonomy=taxonomy), retrievalStatus


if __name__ == "__main__":

    import sys
    inputfile = sys.argv[1]

    #db = DB(inputfile, rebuild=True)
    db = DB(inputfile)

    print db.db.items()
    print
    print db.index.items()

    fastaRecord = Fasta.Record(title='myfasta',
                               sequence='TTTAATGTTAGGAGCTCCAGACATGGCATTCCCTCGTATAAATAATATAAGATTTTGATTATTACCCCCTAAGAAGAATAGTTGAAAGAGGAGC')

    #blastRecord = db.search(fastaRecord, excludelist=['Tachinus brevipennis'])
    blastRecord = db.search(fastaRecord)

    for i in range(len(blastRecord.alignments)):

        gi = blastRecord.alignments[i].title.split(' ')[0]     

        description = blastRecord.descriptions[i]

        print db.get(gi, description.e)

