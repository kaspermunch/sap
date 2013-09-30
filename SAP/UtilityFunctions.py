
try:
   import cPickle as pickle
except:
   import pickle
import re, random, string, time, sys, os, glob, subprocess
# from pkg_resources import load_entry_point, get_entry_info, get_entry_map, iter_entry_points, EntryPoint
from SAP.Bio.Nexus import Nexus
from SAP import Fasta

from SAP import NW

from SAP.Exceptions import AnalysisTerminated

# homedir = os.path.expanduser('~')
# 
# # ...works on at least windows and linux. 
# # In windows it points to the user's folder 
# #  (the one directly under Documents and Settings, not My Documents)
# 
# 
# # In windows, you can choose to care about local versus roaming profiles.
# # You can fetch the current user's through PyWin32.
# #
# # For example, to ask for the roaming 'Application Data' directory:
# #  (CSIDL_APPDATA asks for the roaming, CSIDL_LOCAL_APPDATA for the local one)
# #  (See microsoft references for further CSIDL constants)
# try:
#     from win32com.shell import shellcon, shell            
#     homedir = shell.SHGetFolderPath(0, shellcon.CSIDL_APPDATA, 0, 0)
#  
# except ImportError: # quick semi-nasty fallback for non-windows/win32com case
#     homedir = os.path.expanduser("~")
    
# class Error(Exception):
#     """
#     Base class for exceptions in this module.
#     """
#     pass
#     
# class PluginNotFoundError(Error):
#     """
#     Raised when a plugin is not found.
#     """
#     def __init__(self, plugin):
#         self.plugin = plugin

def _printAlignment(seq1, seq2):

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

# def remap_db_hit(queryRecord, hitRecord, dbHit):
#    
#     # Truncate sequences with flanks of same length as the query on each
#     # side to make sure the query sequence is covered in the alignment
#     # no matter where the homologue maches. 
#     queryLength = len(queryRecord.sequence)
# 
#     leftFlank = 2 * dbHit.query_start
#     rightFlank = 2 * (queryLength - dbHit.query_length - dbHit.query_start)
# 
#     # Check for reverse complementation macth:                            
#     strandMatch = 1
#     if dbHit.query_strand == -1:
#        strandMatch = -1
#        leftFlank, rightFlank = rightFlank, leftFlank
# 
#     startIndex = max(0, dbHit.subject_start - leftFlank)
#     endIndex = min(dbHit.subject_start + dbHit.subject_length + rightFlank, len(hitRecord.sequence))
# 
#     # reverse complement sequence if needed:
#     if strandMatch == -1:
#         hitRecord.sequence = string.translate(hitRecord.sequence[::-1], string.maketrans('AGTC', 'TCAG'))
# 
#     status, alignedQuery, alignedHomol = NW.align(queryRecord.sequence, hitRecord.sequence[startIndex:endIndex])
# 
#     _printAlignment(alignedQuery, alignedHomol)
#     print queryRecord
#     print hitRecord
# 
# #       print dbHit.score/float(len(alignedHomol)), utils.alignmentScore(alignedQuery, alignedHomol)
# 
#     dbHit.score = similarityScore(alignedQuery, alignedHomol)
#     #dbHit.score = alignmentScore(alignedQuery, alignedHomol)
# 
#     leftBoundary = len(re.search("^(-*)", alignedQuery).groups()[0])
#     rightBoundary = len(re.search("(-*)$", alignedQuery).groups()[0])
# 
#     dbHit.query_start -= leftFlank - leftBoundary
#     dbHit.query_length += leftFlank - leftBoundary + rightFlank - rightBoundary
# 
#     dbHit.subject_start -= leftFlank - leftBoundary
#     dbHit.subject_length += leftFlank - leftBoundary + rightFlank - rightBoundary
# 
# #       print leftFlank, leftBoundary, rightFlank, rightBoundary, dbHit.subject_start, dbHit.subject_length
# 
#     return dbHit

# def remap_db_hit(queryRecord, hitRecord, dbHit):
#    
#     # Truncate sequences with flanks of same length as the query on each
#     # side to make sure the query sequence is covered in the alignment
#     # no matter where the homologue maches. 
#     queryLength = len(queryRecord.sequence)
# 
#     flank = 2 * (queryLength - dbHit.query_length)
# 
#     subjectStartIndex = max(0, dbHit.subject_start - flank)
#     subjectEndIndex = min(dbHit.subject_start + dbHit.subject_length + flank, len(hitRecord.sequence))
# 
#     # reverse complement sequence if needed:
#     if dbHit.subject_strand == -1:
#         hitRecord.sequence = string.translate(hitRecord.sequence[::-1], string.maketrans('AGTC', 'TCAG'))
# 
#     status, alignedQuery, alignedSubject = NW.align(queryRecord.sequence, hitRecord.sequence[subjectStartIndex:subjectEndIndex])
#     #_printAlignment(alignedQuery, alignedSubject)
# 
#     dbHit.sequence = alignedSubject.replace("-", "")
# 
#     dbHit.score = similarityScore(alignedQuery, alignedSubject)
#     #dbHit.score = alignmentScore(alignedQuery, alignedSubject)
# 
#     queryStartGaps = len(re.search("^(-*)", alignedQuery).groups()[0])
#     queryEndGaps = len(re.search("(-*)$", alignedQuery).groups()[0])
# 
#     subjectStartGaps = len(re.search("^(-*)", alignedSubject).groups()[0])
#     subjectEndGaps = len(re.search("(-*)$", alignedSubject).groups()[0])
# 
#     dbHit.query_start = subjectStartGaps
#     dbHit.query_length = queryLength - subjectStartGaps - subjectEndGaps
# 
#     subject_start = subjectStartIndex - queryStartGaps
#     subject_end = subjectEndIndex - queryEndGaps
# 
#     dbHit.subject_start = subject_start
#     dbHit.subject_length = subject_end - subject_start
# 
#     return dbHit

def remap_db_hit(queryRecord, dbHit, dbIndex):
   
    # Truncate sequences with flanks of same length as the query on each
    # side to make sure the query sequence is covered in the alignment
    # no matter where the homologue maches. 
    queryLength = len(queryRecord.sequence)
    hitRecordLength = dbIndex.index[dbHit.id].length

    flank = int(1.2 * (queryLength - dbHit.query_length))
    subjectStartIndex = max(0, dbHit.subject_start - flank)
    subjectEndIndex = min(dbHit.subject_start + dbHit.subject_length + flank, hitRecordLength)

    subjectCutoutLength = subjectEndIndex - subjectStartIndex
    cutout = dbIndex.get(dbHit.id, subjectStartIndex, subjectCutoutLength, dbHit.subject_strand)

    status, alignedQuery, alignedSubject = NW.align(queryRecord.sequence, cutout, 5, 2, 1)

    _printAlignment(alignedQuery, alignedSubject)
#     sys.exit()
    
    dbHit.sequence = alignedSubject.replace("-", "")

    dbHit.score = similarityScore(alignedQuery, alignedSubject)
    #dbHit.score = alignmentScore(alignedQuery, alignedSubject)

    queryStartGaps = len(re.search("^(-*)", alignedQuery).groups()[0])
    queryEndGaps = len(re.search("(-*)$", alignedQuery).groups()[0])

    subjectStartGaps = len(re.search("^(-*)", alignedSubject).groups()[0])
    subjectEndGaps = len(re.search("(-*)$", alignedSubject).groups()[0])

    dbHit.query_start = subjectStartGaps
    dbHit.query_length = queryLength - subjectStartGaps - subjectEndGaps

    subject_start = subjectStartIndex - queryStartGaps
    subject_end = subjectEndIndex - queryEndGaps

    dbHit.subject_start = subject_start
    dbHit.subject_length = subject_end - subject_start

    return dbHit


def systemCall(cmd, stdout=None, stderr=None):
    """
    Executes a system call without a shell/CMD. Takes only single comands with no
    redirection or ';' chars. Used where we don't want CMD to pop up on windows.
    """
    if os.name == 'nt':
        null = 'nul'
    else:
        null = '/dev/null'

    if stdout == 'IGNORE':
        stdout = open(null, 'w')

    if stderr == 'IGNORE':
        stderr = open(null, 'w')

    STARTUPINFO = None
    if os.name == 'nt':
       STARTUPINFO = subprocess.STARTUPINFO()
       STARTUPINFO.dwFlags |= subprocess.STARTF_USESHOWWINDOW

    cmdlist = re.split(r'\s+', cmd)
    
    try:
        retcode = subprocess.call(cmdlist, shell=False, env=os.environ, stdout=stdout, stderr=stderr, startupinfo=STARTUPINFO)
        if retcode < 0:
            print 'System call "%s" was terminated with return code' % cmd, -retcode
            return False
        else:
            return retcode
    except OSError, e:
        print 'System call "%s" failed:' % cmd, e

    stdout.close()
    stderr.close()

def pairwiseClustalw2(id1, sequence1, id2, sequence2):
    """
    Make a pairwise alignment using a system call to clustalw2 and
    returns a Nexus alignment object.
    """
    tmpAlnFileName = os.path.abspath(randomString(6) + ".tmp")
    tmpFastaFileName = tmpAlnFileName + '.fasta'
    tmpFastaFileContents = ">%s\n%s\n>%s\n%s\n" % (id1, sequence1, id2, sequence2)
    writeFile(tmpFastaFileName, tmpFastaFileContents)
#     if os.name == 'nt':
#         commandLine = "clustalw2 %s -output=NEXUS -gapopen=50 -outfile=%s > %s.log" \
#                       % (os.path.basename(tmpFastaFileName), os.path.basename(tmpAlnFileName), os.path.basename(tmpAlnFileName))
#     else:
#         commandLine = "clustalw2 -infile=%s -output=NEXUS -gapopen=50 -outfile=%s &> %s.log" \
#                       % (os.path.basename(tmpFastaFileName), os.path.basename(tmpAlnFileName), os.path.basename(tmpAlnFileName))
#     os.system(commandLine)
    if os.name == 'nt':
        commandLine = "clustalw2 %s -output=NEXUS -gapopen=50 -outfile=%s" \
                      % (os.path.basename(tmpFastaFileName), os.path.basename(tmpAlnFileName))
    else:
        commandLine = "clustalw2 -infile=%s -output=NEXUS -gapopen=50 -outfile=%s" \
                      % (os.path.basename(tmpFastaFileName), os.path.basename(tmpAlnFileName))
    systemCall(commandLine, stdout='IGNORE', stderr='IGNORE')
    #os.system(commandLine)

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

# def findPlugin(name, entry_point_name):
#     plugin = None
#     for entryp in iter_entry_points(entry_point_name):
#        if entryp.name == name:
#           plugin = entryp.load()
#           break
#     if not plugin:
#        raise PluginNotFoundError, name
#     return plugin

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
    #width = 100
    width = alnLength
    end = width
    entry = "#NEXUS\nbegin data;\n        dimensions ntax=%d nchar=%d;\n        format datatype=dna missing=? gap=-;\nmatrix\n"  % (len(seqList), alnLength)

    # Find the longest name:
    names = []
    seqs = []
    longestName = 0
    for fastaRecord in seqList:
        assert len(fastaRecord.sequence) == alnLength, "SEQUENCES NOT SAME LENGTH: %s != %s" % (len(fastaRecord.sequence), alnLength)
        names.append(fastaRecord.title)
        seqs.append(fastaRecord.sequence)
        longestName = max(len(fastaRecord.title), longestName)
        
    for start in range(0, alnLength, width):
        for i in range(len(names)):
            s = seqs[i][start:min(end, alnLength)]
            entry += "%- *s %s\n" % (longestName, names[i], s)         
        end += width
        entry += "\n"
    entry += ";\nend;\n"
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
    chars = string.ascii_lowercase + string.digits #  We use only lowercase to make names cross-platform compatible...
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
        raise AnalysisTerminated(1, "Could not find/read file: %s" % fileName)
#         print "Could not find/read file: %s" % fileName
#         sys.exit(1)
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
        raise AnalysisTerminated(1, "Could not find/read file: %s" % fileName)
#         print "Could not find/read file: %s" % fileName
#         sys.exit(1)
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

def safeReadFastaCache(fastaFileName):
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

def safeReadTaxonomyCache(taxonomyFileName):
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
