

from Bio import Entrez

from SAP import Fasta

from SAP.UtilityFunctions import *



genus = "Turdus"

speciesList = ["amaurochalinus", "migratorius", "rufiventris"]

tmpAlignmentFileName = 'kasper'
outputTmpFileName = "tmpkasper"

sequenceLists = {}

fastaFileName = "kasper.fasta"
fastaFile = open(fastaFileName, 'w')

for species in speciesList:

    handle = Entrez.esearch(db="nucleotide", retmax=10, term="%s %s[ORGN] AND barcode[keyword]" % (genus, species))
    record = Entrez.read(handle)
    
    #print record["Count"]
    
    handle = Entrez.efetch(db="nucleotide", id=','.join(record["IdList"]), rettype="fasta", retmax=10)
    
    fastaIterator = Fasta.Iterator(handle, Fasta.RecordParser())
    for entry in fastaIterator:
        entry.title = entry.title.split('|')[1]
        sequenceLists.setdefault(species, []).append(entry)
        
        fastaFile.write(str(entry))

fastaFile.close()

alignmentoptions = "-gapopen=50"

commandLine = "clustalw2 -infile=%s -output=NEXUS -outfile=%s %s" % (fastaFileName, outputTmpFileName, alignmentoptions)
systemCall(commandLine, stdout='IGNORE', stderr='IGNORE')

writeFile(tmpAlignmentFileName, readFile(outputTmpFileName))
os.unlink(outputTmpFileName)

alignment = Nexus.Nexus(tmpAlignmentFileName)

alignmentLength = None
alignmentString = ""
for gi in alignment.matrix.keys():
    
    alignmentString += "%-10s %s\n" % (gi, alignment.matrix[gi].tostring())
    alignmentLength = len(alignment.matrix[gi].tostring())

tree = "((0,1):3,2):4"
        
s = genus + "\n"
s += "# For finding appropriate priors\n"
s += "3\n"
s += " ".join(speciesList) + "\n"
s += str(tree) + "\n"
s += "1\n"
s += "%s 0" % genus
for species in speciesList[1:]:
    s += " %d" % len(sequenceLists[species])
s += " %d H 1.000000 A1\n%s" % (alignmentLength, alignmentString)


print s



