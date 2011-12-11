import unittest, sys

from SAP import ConsoleScripts
from SAP import Options

class TestSAP(unittest.TestCase):

    def setUp(self):
        """
        make a temp dir for test files. write a query and a database file.
        """
        self.projectDir = ""
        self.databaseFileName = ""
        self.args = ["--project", self.projectDir]
        pass

    def test_Default(self):
        ConsoleScripts.sap()
        pass

    def test_NativeDB(self):
#         sys.argv = self.args + ["--database", self.databaseFileName, queryFileName]
#         ConsoleScripts.sap()
        pass
    
    def test_ConstrainedNJ(self):
#         sys.argv = self.args + ["--assignmnet", "ConstrainedNJ"]
        pass

    def test_IMa2(self):
#         sys.argv = self.args + ["--ghostpopulation"]
        pass

    def tearDown(self):
        """
        delete any test files.
        """
        pass

if __name__ == '__main__':
    unittest.main()
