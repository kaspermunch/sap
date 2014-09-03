
import re

class Record(object):
    """
    Simple Record class to bundle attributes.
    """
    def __init__(self, title=None, sequence=None):

        self.title = title
        self.sequence = sequence
        self.lineLength = 100
        
    def __str__(self):
        length = len(self.sequence)
        s =  ">%s\n" % self.title
        idx = 0

        while length > idx:
            s += self.sequence[idx:min(length,idx+self.lineLength)] + "\n"
            idx += self.lineLength
        return s

class RecordParser(object):
    """
    Dummy class to make the interface the same as that of Fasta in Biopython.
    """
    def __init__(self):
        pass

def _makeRecord(block):
    """
    Parses a fasta entry block and returns a Record object.
    """
    end_title = block.find("\n")
    return Record(block[1:end_title], block[end_title:].replace("\n","").replace(" ", ""))
    
def Iterator(fileob, parser=None):
    """
    The parser argument is only used to give the module the same interface as Fasta in
    Biopython
    """
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
            yield _makeRecord(block)
            return

        where += seplen / 2
        yield _makeRecord(block[:where])
        block = block[where:]


if __name__ == '__main__':

    fob = open('../../testing/release.fasta', 'r')

    it = Iterator(fob)
    record = it.next()
    print record.title
    print record.sequence

#     for record in Iterator(fob):
#         print record
    
