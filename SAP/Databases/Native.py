#!/usr/bin/python

import sys, os, re, shelve, tempfile, StringIO
from sets import Set
from SAP import Fasta, Taxonomy, Homology
from SAP.UtilityFunctions import *
from SAP import NCBIXML # locally hacked to better parse string info
try:
    import cPickle as pickle
except ImportError:
    import pickle

class DB:

    def __init__(self, fastaFileName, options, rebuild=False):
        """
        Opens a premade local data base and builds an index of taxon names.
        """
        if not findOnSystem('blastn'):
            print 'blastn is not available'
            sys.exit()
        if not findOnSystem('xdformat'):
            print 'xdformat is not available'
            sys.exit()

        self.options = options
        self.fastaFileName = fastaFileName

        baseName, ext = os.path.splitext(os.path.basename(self.fastaFileName))
        self.baseName = baseName
        self.dbDirName = baseName + '.sapdb'
        self.blastSequenceFileName = os.path.join(self.dbDirName, baseName + '.fasta')
        self.blastDB = os.path.join(self.dbDirName, baseName)
        self.dbFileName = os.path.join(self.dbDirName, baseName + '.db')
        self.indexFileName = os.path.join(self.dbDirName, baseName + '.idx')

        # Build shelve database if not present:
        if rebuild or not os.path.exists(self.dbDirName):
            self._build()
        else:
            # Just open the db:
            self.db = shelve.open(self.dbFileName, 'r')
            indexFile = open(self.indexFileName, 'r')
            self.index = pickle.load(indexFile)


    def _build(self):
        """
        Build a local database from a file of fasta entries.
        """
        # Make a dir for the database files:
        if not os.path.exists(self.dbDirName):
            os.mkdir(self.dbDirName)

        # Create a new shelve:
        db = shelve.open(self.dbFileName, 'n')

        # Create an iterator over fasta entries:
        fastaFile = open(self.fastaFileName, 'r')
        fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        

#         # Initialize sequenceMatrix, sequenceMatrixNames and bootstrapList for creating cached bootstraps:
#         matrixNames = []
#         matrix = []
#         bootstrapList = []

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

        firstMissedTaxonomy = True
        for fastaRecord in fastaIterator:
            try:
                # See if we can get the taxonomy form the fasta title line as we are supposed to:
                taxonomy = Taxonomy.Taxonomy()
                identifier, taxonomyString, organismName = re.split(r'\s*;\s*', fastaRecord.title)
                taxonomyString.strip()
                fastaRecord.title = identifier        
                taxonomy.populateFromString(taxonomyString)
            except ValueError:                
                # if not - load the genbank plugin and retrieve from genbank
                if self.genbank is None:
                    plugin = findPlugin('GenBank', 'sap.database')
                    genbank = plugin.DB(self.options)
                identifier = fastaRecord.title
                if not re.match(r'\d+', identifier):
                    raise Exception
                if firstMissedTaxonomy:
                    print "Retrieving taxonomy from NCBI's Taxonomy Browser:\n\t%s" % identifier
                    firstMissedTaxonomy = False
                else:
                    print "\t%s" % identifier

                homologue, ret = genbank.get(identifier, None)
                taxonomy = homologue.taxonomy
                organismName = taxonomy.organism

            # Make sure we don't have doubles:
            assert not db.has_key(identifier)

#             # Add the sequence to the sequence matrix:
#             if alignmentLength == None or len(fastaRecord.sequence) == alignmentLength:
#                 matrixNames.append(fastaRecord.title)
#                 matrix.append(fastaRecord.sequence)
#                 alignmentLength = len(fastaRecord.sequence)
#             else:
#                 isPreAligned = False

            # Create a database entry:
            db[str(identifier)] = {'fastaRecord': fastaRecord,
                                   'taxonomy': taxonomy,
                                   'organismName': organismName}

            # Write the sequence without gaps to the blast file:
            fastaRecord.sequence = fastaRecord.sequence.replace('-', '')
            blastSequenceFile.write(str(fastaRecord))

            # Add an entry to the index:
            for taxonomyLevel in taxonomy:
                index.setdefault(taxonomyLevel.name, []).append(identifier)

        # Syncronize shelve:
        db.sync()
        self.db = db

        # Pickle the index:
        indexFile = open(self.indexFileName, 'w')
        pickle.dump(index, indexFile)
        self.index = index

        # Format the blast database:
        blastSequenceFile.close()
        cmd = "xdformat -n -o %s %s" % (self.blastDB, self.blastSequenceFileName)
        stdin, stdout, stderr = os.popen3(cmd)            

#         if isPreAligned:
# 
#             # Generate a list of index lists for bootstrapping:
#             bootstrapIndexList = []
#             alignmentLength = sequenceMatrix
#             for i in range(1000):
#                 tmpList = []
#                 for j in range(alignmentLength):            
#                     tmpList.append(random.randint(0, alignmentLength-1))
#                 bootstrapIndexList.append(tmpList)
# 
#             # Generate 1000 distance matrices:
#             for idxList in bootstrapIndexList:
#                 distanceMatrix = []
#                 for i in len(sequenceMatrix):
#                     distanceMatrix.append([])
#                     for j in alignmentLength:
#                         distanceMatrix[i].append(distance(i, j))
# 
#             # Dump the bootstrapList and matrixNames:




    def search(self, fastaRecord, excludelist=[], usecache=True):

        print "\n\t\tSearching database...", 
        sys.stdout.flush()

        # Write the query to a tmp file:
        tmpDirName = tempfile.mkdtemp()
        tmpQueryFile, tmpQueryFileName = tempfile.mkstemp(dir=tmpDirName)
        writeFile(tmpQueryFileName, str(fastaRecord))

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
            stdin, stdout, stderr = os.popen3(cmd)
        else:
            blastDBfileName = self.blastDB
            
        # Blast:
        cmd = "blastn %s %s -mformat 7" % (blastDBfileName, tmpQueryFileName)
        stdin, stdout, stderr = os.popen3(cmd)
        blastParser = NCBIXML.BlastParser()
        try:
            blastRecord = blastParser.parse(StringIO.StringIO(stdout.read())) 
        except:
            blastRecord = None        

        print "done.\n\t\t\t",
        sys.stdout.flush()

        return blastRecord
        
    def get(self, gi, evalue):

        fastaRecord = self.db[str(gi)]['fastaRecord']
        taxonomy = self.db[str(gi)]['taxonomy']

        retrievalStatus = '(n)'

        return Homology.Homologue(gi=gi,
                                  sequence=fastaRecord.sequence,
                                  evalue = evalue, 
                                  taxonomy=taxonomy,
                                  options=None), retrievalStatus


if __name__ == "__main__":

    import sys
    inputfile = sys.argv[1]

    #db = DB(inputfile, rebuild=True)
    db = DB(inputfile)

#     print db.db.items()
#     print
#     print db.index.items()

    fastaRecord = Fasta.Record(title='myfasta', sequence='TTTAATGTTAGGAGCTCCAGACATGGCATTCCCTCGTATAAATAATATAAGATTTTGATTATTACCCCCTTCATTATCTTTATTATTAATAAGAAGAATAGTTGAAAGAGGAGC')

    #db.search(fastaRecord, excludelist=['Tachinus brevipennis'])
    blastRecord = db.search(fastaRecord)

    for i in range(len(blastRecord.alignments)):

        # There are some problems here: the gi is atracted differently than what is hardwired into Homology.py ... We need to parse it out in the plugin and not in Homology.py

        gi = blastRecord.alignments[i].title.split(' ')[0]     

        description = blastRecord.descriptions[i]

        #(homologue, retrievalStatus) = db.get(gi, description.e)
        print db.get(gi, description.e)

