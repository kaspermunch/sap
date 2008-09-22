
import re, random, string, time, sys, os, glob
from pkg_resources import load_entry_point, get_entry_info, get_entry_map, iter_entry_points, EntryPoint
from SAP.Bio.Nexus import Nexus

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass
    
class PluginNotFoundError(Error):
    """
    Raised when a plugin is not found.
    """
    def __init__(self, plugin):
        self.plugin = plugin
        
def pairwiseClustalw2(id1, sequence1, id2, sequence2):
    """
    Make a pairwise alignment using a system call to clustalw2 and
    returns a Nexus alignment object.
    """
    tmpAlnFileName = os.path.abspath(randomString(6) + ".tmp")
    tmpFastaFileName = tmpAlnFileName + '.fasta'
    tmpFastaFileContents = ">%s\n%s\n>%s\n%s\n" % (id1, sequence1, id2, sequence2)
    writeFile(tmpFastaFileName, tmpFastaFileContents)
    if os.name == 'nt':
        commandLine = "clustalw2 %s -output=NEXUS -gapopen=50 -outfile=%s > %s.log" \
                      % (os.path.basename(tmpFastaFileName), os.path.basename(tmpAlnFileName), os.path.basename(tmpAlnFileName))
    else:
        commandLine = "clustalw2 -infile=%s -output=NEXUS -gapopen=50 -outfile=%s &> %s.log" \
                      % (os.path.basename(tmpFastaFileName), os.path.basename(tmpAlnFileName), os.path.basename(tmpAlnFileName))
    os.system(commandLine)

    # small bug fix:
    nexusContents = readFile(tmpAlnFileName)
    nexusContents = re.sub("(symbols=\"[^\"]*\")", r"[\1]", nexusContents)
    alignment = Nexus.Nexus(nexusContents)
    # Remove tmp alignment files:
    for f in glob.glob(tmpAlnFileName + '*'):
        if os.path.exists(f):
            os.remove(f)
    return alignment

def findOnSystem(filename):
    """
    Finds specified filename in system path.
    """
    fileFound = 0
    paths = []
    if os.environ.has_key('PATH'):
        paths = os.environ['PATH'].split(os.pathsep)

    for path in paths:
        if os.path.exists(os.path.join(path, filename)):
            fileFound = 1
            break

    # Look in sensible places:
    if not fileFound:

        if os.name == 'posix':
            # In OSX applications the PATH set by user shell is not available:
            for path in [os.path.join(os.environ['HOME'], 'bin'), os.path.join(os.environ['HOME'], 'usr/local/bin'),'/usr/local/bin/']:
                if os.path.exists(os.path.join(path, filename)):
                    paths.append(path)
                    ## os.putenv('PATH', os.pathsep.join(paths))
                    os.environ['PATH'] = os.pathsep.join(paths)
                    fileFound = 1
                    break

        if os.name == 'nt':
            # On windows look in directories under "Program Files":
            for path in glob.glob(r'C:\Program Files\*'):
                if os.path.exists(os.path.join(path, filename)):                    
                    paths.append(path)
                    ## os.putenv('PATH', os.pathsep.join(paths))
                    os.environ['PATH'] = os.pathsep.join(paths)
                    fileFound = 1
                    break
                
    if fileFound:
        return os.path.abspath(path)
    else:
        return None

def findPlugin(name, entry_point_name):
    plugin = None
    for entryp in iter_entry_points(entry_point_name):
       if entryp.name == name:
          plugin = entryp.load()
          break
    if not plugin:
       raise PluginNotFoundError, name
    return plugin

def similarityScore(seq1, seq2):
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


def poolStatus(pool):

    output = []
    error = []
    for q in pool.getQueuedOutput():
        output.append(q)
    for q in pool.getQueuedError():
        error.append(q)
    if output or error:
        print "\nUpdate on parallel jobs."
        if output:
            print "\tOutput:"
            for o in output:
                print "\t\t" + o.replace('\n', '\n\t\t').strip()
        if error:
            print "\tError:"
            for e in error:
                print "\t\t" + e.replace('\n', '\n\t\t').strip()
        sys.stdout.flush()

def writeNexusFile(fileName, seqList):
    """
    Takes a list if Fasta objects and writes a phylip formatted
    entry to a file.
    """
    fp = open(fileName, 'w')
    alnLength = len(seqList[0].sequence)
    width = 100
    #width = alnLength
    end = width
    entry = "#NEXUS\nbegin data;\n        dimensions ntax=%d nchar=%d;\n        format datatype=dna missing=? gap=-;\nmatrix\n"  % (len(seqList), alnLength)

    # Find the longest name:
    names = []
    seqs = []
    longestName = 0
    for fastaRecord in seqList:
        assert len(fastaRecord.sequence) == alnLength, "SEQUENCES NOT SAME LENGTH"
        names.append(fastaRecord.title)
        seqs.append(fastaRecord.sequence)
        longestName = max(len(fastaRecord.title), longestName)
        
    for start in range(0, alnLength, width):
        for i in range(len(names)):
            s = seqs[i][start:min(end, alnLength)]
            entry += "%- *s %s\n" % (longestName, names[i], s)         
        end += width
        entry += "\n"
    entry += "end;\n"
    fp.write(entry)
    fp.close()


def writePhylipFile(fileName, seqList):
    """
    Takes a list of (name, Seq object) tuples and writes a phylip formatted
    entry to a file.
    """
    fp = open(fileName, 'w')
    alnLength = len(seqList[0][1].data)
    #width = 100
    width = alnLength
    end = width
    entry = "%d %d\n" % (len(seqList), alnLength)

    # Find the longest name:
    names = []
    seqs = []
    longestName = 0
    for name, seq in seqList:
        seq = seq.tostring()
        assert len(seq) == alnLength, "SEQUENCES NOT SAME LENGTH"
        names.append(name)
        seqs.append(seq)
        longestName = max(len(name), longestName)

    for start in range(0, alnLength, width):
        for i in range(len(seqList)):
            s = seqs[i][start:min(end, alnLength)]
            if start == 0:
                entry += "%- *s %s\n" % (longestName, names[i], s)
            else:
                entry += "%- *s %s\n" % (longestName, '', s)
        end += width
        entry += "\n"
    fp.write(entry)
    fp.close()

# def writePhylipFile(fileName, seqList):
#     """
#     Takes a list of (name, Seq object) tuples and writes a phylip formatted
#     entry to a file.
#     """
#     fp = open(fileName, 'w')
#     alnLength = len(seqList[0][1].data)
#     #width = 100
#     width = alnLength
#     end = width
#     entry = "%d %d\n" % (len(seqList), alnLength)
#     for start in range(0, alnLength, width):
#         for name, seq in seqList:
#             s = seq[start:min(end, alnLength)].tostring()
#             if start == 0:
#                 entry += "%- 30s %s\n" % (name, s)
#             else:
#                 entry += "%- 30s %s\n" % ('', s)
#         end += width
#         entry += "\n"
#     fp.write(entry)
#     fp.close()


def randomString(length):
    chars = string.ascii_lowercase + string.digits
    rand = ''
    for i in range(length):
        rand = rand + random.choice(chars)
    return rand

def listOrTuple(x):
    return isinstance(x, (list, tuple))


def flatten(sequence, toExpand=listOrTuple):
    for item in sequence:
        if toExpand(item):
            for subitem in flatten(item, toExpand):
                yield subitem
        else:
            yield item

def safeName(name):
    """
    Change the fasta title to something safe. E.g. MrBayes does not like dashes.
    """    
    name = name.strip()
    name = re.sub(r'[-/\\\(\)|<>.,:;\[\] ]', '_', name)
    name = re.sub(r'[!@#$%^&*+?~`"\']', '', name)
    return name

def readFile(fileName):
    """
    Read in contents of the file given by fileName
    """
    fp = None
    for tries in range(10):
        try:
            fp = open(fileName, "r")
        except IOError:
            time.sleep(tries * 5)
            continue
        break
    if fp is None:
        print "Could not find/read file: %s" % fileName
        sys.exit(1)
    else:
        contents = fp.read()
        fp.close()
        return contents

def readFileLines(fileName):
    """
    Read in contents of the file given by fileName
    """
    fp = None
    for tries in range(10):
        try:
            fp = open(fileName, "r")
        except IOError:
            time.sleep(tries * 5)
            continue
        break
    if fp is None:
        print "Could not find/read file: %s" % fileName
        sys.exit(1)
    else:
        contents = fp.readlines()
        fp.close()
        return contents

def writeFile(fileName, contents):
    """
    Write contents to the file given by fileName
    """
    try:
        file = open(fileName, 'w')
        file.write(contents)
        file.close()
    except:
        if os.path.exists(fileName):
            os.remove(fileName)

def copyCacheForSequenceDoubles(dict, options):
    """
    Copies the cached result files for the sample-sequences that
    occour more than once, because we only the analysis on one
    representative.
    """
    for copyFrom in dict.keys():
        # Find all the cache files related to copyFrom:
        for cacheDir in [options.blastcache, options.homologcache,
                         options.alignmentcache, options.treescache, options.treestatscache]:
            for copyFile in glob.glob(os.path.join(cacheDir, copyFrom + '.*')):
                path, fileName = os.path.split(copyFile)
                baseName, suffix = os.path.splitext(fileName)
                # Write a copy of each file:
                for copyTo in dict[copyFrom]:
                    #newFileName = os.path.join(path, copyTo + suffix)
                    newFileName = os.path.join(path, fileName.replace(copyFrom + '.', copyTo + '.'))
                    s = readFile(copyFile)
                    s = s.replace(copyFrom, copyTo)
                    writeFile(newFileName, s)


if __name__ == "__main__":


    import Fasta
    from SAP.Bio import Seq
    fastaFileName = sys.argv[1]
    if not os.path.exists(fastaFileName):
        raise Exception

    fastaFile = open(fastaFileName, 'r')
    fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        
    seqList = []
    for r in fastaIterator:
        seqList.append(r)

    writeNexusFile('tmp.nex', seqList)
# 
#     fastaFile = open(fastaFileName, 'r')
#     fastaIterator = Fasta.Iterator(fastaFile, parser=Fasta.RecordParser())        
#     seqList = []
#     for r in fastaIterator:
#         seqList.append((r.title, Seq.Seq(r.sequence)))
# 
#     writePhylipFile('tmp.phylip', seqList)
