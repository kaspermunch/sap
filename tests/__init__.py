import unittest, sys, shutil

from SAP import ConsoleScripts
from SAP import Options

class TestSAP(unittest.TestCase):

    def setUp(self):
        """
        make a temp dir for test files. write a query and a database file.
        """
        self.queryFileName = 'tests/query.fasta'
        self.projectDir = "testproject"
        self.databaseFileName = "tests/testdb.fasta"

    def test_1_dependencies(self):
	sys.argv = ['', '--compile', 'COI[Gene Name] AND Aves[ORGN]', '--database', self.databaseFileName]
        ConsoleScripts.sap()

    def test_2_dependencies(self):
        sys.argv = ['', '--installdependencies']
        ConsoleScripts.sap()

    def test_3_native_with_barcoder(self):
        sys.argv = ['', "-d", self.projectDir, "--email", "test@test.com", "-S", "Barcoder", "--database", self.databaseFileName, self.queryFileName]
        ConsoleScripts.sap()

    def test_4_genbank_with_cnj(self):
        sys.argv = ['', "--project", self.projectDir, "--email", "test@test.com", "-S", "ConstrainedNJ", "--database", "GenBank", self.queryFileName]
        ConsoleScripts.sap()
    
#     def test_4_ima2(self):
# #         sys.argv = self.args + ['', "--ghostpopulation"]
#         pass

    def tearDown(self):
        """
        delete any test files.
        """
#        shutil.rmtree(self.projectDir)

if __name__ == '__main__':
    unittest.main()
