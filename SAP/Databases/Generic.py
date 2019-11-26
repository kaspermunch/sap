

# # An external plugin package should include a class named DB supplying
# # the methods shown below. The package should contain a file with the
# # code e.g. MyPlugInModule.py and a file named setup.py with the
# # following content:
#  
# import ez_setup
# ez_setup.use_setuptools()
# from setuptools import setup, find_packages
# setup(name='MyPlugInPackage>',
#       version='1.0',
#       description='',
#       long_description='',
#       author='Your Name',
#       author_email='your@email.com',
#       packages = find_packages(exclude=['ez_setup']),
#       entry_points = { 'sap.database': [ 'MyPlugInName = MyPlugInModule', ], },
#       )
# 
# # Then the develper of the plugin just uploads to PyPI:
# setup.py bdist_egg upload         # create an egg and upload it
# setup.py sdist upload             # create a source distro and upload it
# setup.py sdist bdist_egg upload   # create and upload both
# 
# # Now everyone can install the plugin like this:
# easy_install MyPlugInPackage


try:
   import cPickle as pickle
except:
   import pickle
import os, sys, time, re, pickle

from SAP import Fasta
from Bio.EUtils.Datatypes import DBIds
from Bio.EUtils.ThinClient import ThinClient

from SAP import Homology # I should rename this to Homology

class DB(object):

    def __init__(self):
        """
        Instantiate the database, build it if necessary.
        """
        pass


    def search(self, fastaRecord, excludelist=[], usecache=True):
        """
        Called when the database is searched. Must return a blast record.
        """

        return SAP.Bio.Blast.Record()

    def get(self, gi):
       """
       Called to retrieve a homologue and its taxonomic annotation
       from the database. Must return a fully populated
       HomolCompiler.Homologue object.
       """

       return HomolCompiler.Homologue(), retrievalStatus
  
