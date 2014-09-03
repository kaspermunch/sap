import re
import sys, os, string
import cPickle as pickle
from SAP.Fasta import Record

from SAP.ProgressBar import ProgressBar

class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass


class IndexingError(Error):
    """
    Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which
                      the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        print self.expression, ": ", self.message


def fastaGenerator(fileob):

    separator = "\n>"
    seplen = len(separator)
    block = fileob.read(8192)

    while block:
        where = block.find(separator)

        if where < 0:
            moredata = fileob.read(8192)
            if moredata:
                block += moredata
                continue
            yield block
            return

        where += seplen / 2
        yield block[:where]
        block = block[where:]

class FastaIndex(object):

    def __init__(self, inputFileName, rebuild=False, progress=False):

        self.indexFileName = inputFileName + ".idx"
        self.inputFile = open(inputFileName, 'r')

        # Check for existence of an index file:
        if rebuild or not os.path.exists(self.indexFileName):

            self.index = {}
            regex = re.compile(r'>(\S+)')
            pos = 0
            self.keys = set()
            for chunk in fastaGenerator(self.inputFile):
                key = regex.match(chunk).group(1)
                self.keys.add(key)
                if key in self.index:
                    raise IndexingError(key, "Sequence ids must be unique")
                self.index[key] = (pos, len(chunk))
                pos += len(chunk)

            indexFile = open(self.indexFileName, 'w')
            pickle.dump(self.index, indexFile)
            indexFile.close()

        else:
            indexFile = open(self.indexFileName, 'r')
            self.index = pickle.load(indexFile)
            self.keys = set(self.index.keys())
            indexFile.close()

    def get_entry(self, seq_id):

        start, length = self.index[seq_id]

        self.inputFile.seek(start)
        block = self.inputFile.read(length)
        end_title = block.find("\n")

        return Record(block[1:end_title], block[end_title:].replace("\n","").replace(" ", ""))

    def keys(self):
        return self.keys

    def __in__(self, key):
        return key in self.keys


if __name__ == "__main__":

    inputFileName = sys.argv[1]

    idx = FastaIndex(inputFileName)

    if len(sys.argv) == 3:
        print idx.get_entry(sys.argv[2])
