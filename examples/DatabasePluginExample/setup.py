
import ez_setup
ez_setup.use_setuptools()
from setuptools import setup, find_packages
setup(name='MyPlugInPackage>',
      version='1.0',
#       description='',
#       long_description='',
#       author='Kasper Munch',
#       author_email='kasmunch@bi.ku.dk',
#       install_requires = ['biopython >= 1.3'],
#       dependency_links = ["http://biopython.org/DIST"],
      packages = find_packages(exclude=['ez_setup']),
#       package_dir = {'Taxoniphy': 'Taxoniphy'},
#       include_package_data = True,
      entry_points = { 'taxoniphy.database': [ 'MyPlugIn = MyPlugInModule', ], },
      )
