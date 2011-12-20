import sys,os

class Table:
    def __init__(self):
        self.data = dict()
        self.rows = 0;
            
    def add_columns(self, **vals):
        for k in vals:
            assert not self.rows or self.rows == len(vals.values()[0])
            assert k not in self.data
            self.data[k] = vals[k]
        return self

    def load_kv(self, file):
        self.data = dict()
        self.rows = 1;

        fin = open(file, "r")
        for line in fin :
            k = line[0:line.find("=")].strip()
            v = line[line.find("=")+1:len(line)]
            self.data[k] = [v]
        
        fin.close()

        return self

    def load(self, file):
        self.data = dict()
        self.rows = 0
        header = True

        if not os.path.exists(file):
            return self


        cols = None
        
        fin = open(file, "r")
        for line in fin :
            if header:
                cols = line.strip().split('\t')
                for c in cols:
                    self.data[c.strip()] = []

                header = False
            else:
                d = line.strip().split('\t')
                for i in range(0, len(cols)):
                    self.data[cols[i].strip()].append(d[i])

                self.rows = self.rows + 1

        return self
            

    def has(self, **vals):
        for r in range(0, self.rows):
            good = True
            for k in vals:
                if k not in self.data:
                    return False

                if str(self.data[k][r]).strip() != str(vals[k]).strip():
                    good = False

            if good:
                return True

        return False

    def merge(self, other):
        res = Table()
        #copy self
        for k in self.data:
            res.data[k] = []
            for d in self.data[k]:
                res.data[k].append(d)

        for k in other.data:
            #fill out missing columns
            if k not in self.data:
                res.data[k] = []
                for i in range(0, self.rows):
                    res.data[k].append("NULL")

            #append other table
            for d in other.data[k]:
                res.data[k].append(d)

        for k in self.data:
            #fill other missing columns
            if k not in other.data:
                for i in range(0, other.rows):
                    res.data[k].append("NULL")

        res.rows = self.rows + other.rows

        return res

    def add_row(self, **vals):
        tmp = Table()
        for k in vals:
            tmp.data[k] = [vals[k]]

        tmp.rows = 1

        return self.merge(tmp)

    def join(self, other):
        if self.rows != other.rows:
            raise Exception("does not compute")

        res = Table()

        for k in self.data:
            res.data[k] = self.data[k]

        for k in other.data:
            if k not in res.data:
                res.data[k] = other.data[k]


        res.rows = self.rows
        return res
            
    def get_row(self, i):
        res = []
        for k in self.data:
            res.append(str(self.data[k][i]).strip())

        return res
        

    def get_header(self):
        res = []
        for k in self.data:
            res.append(k.strip())

        return res

    def __str__(self):
        res = "\t".join(self.get_header()) + "\n"
        #print res
        for i in range(0, self.rows):
            r = "\t".join(self.get_row(i))
            #print r
            res = res + r + "\n"
        return res

    def write(self, filename):
        f = open(filename, "w")
        f.write(str(self))
        f.close()
