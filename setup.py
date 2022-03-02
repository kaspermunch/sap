
# # You're right, the environment variables weren't getting imported.  So I
# # just added it to my python code instead:
# import os; os.environ['CC'] = 'g++'; os.environ['CXX'] = 'g++';
# os.environ['CPP'] = 'g++'; os.environ['LDSHARED'] = 'g++'
# 
# # I found that CC determined the 1st command and LDSHARED determined the
# # second one (gcc vs. g++).  But strangely enough, neither (CC=g++ and
# # LDSHARED=g++) or (CC=g++ and LDSHARED=gcc) worked.  But if I do (CC=gcc
# # and LDSHARED=gcc) then manually redo the 2nd g++, then it works.



# bdist_mpkg (to create a osx gui installer) does not currently work on Leopard - but may at some point.
import ez_setup
ez_setup.use_setuptools()

import sys, glob, os
import platform
from setuptools import setup, Extension, find_packages

if sys.version < '2.4' or sys.version >= '3.0':
    print "This package requires Python version 2.4 or higher and will not run with python 3. Aborting."
    print "Visit http://www.python.org/download for a suitable version of Python."
    print "Aborting."
    sys.exit()

guiscript = 'SAP/GUI.py'

# grep for r'version="%prog\s+[\d.]+" in SAP/Options.py
version = "1.9.10"

if sys.platform == 'darwin':
    # Cross-platform applications generally expect sys.argv to be used for opening files:
    extra_options = dict(app=[guiscript],
                         # setup_requires=['py2app'],
                         options = dict(py2app = dict( argv_emulation=True,
                                                       iconfile='icons/app.icns'),                                       
                                        plist = dict( CFBundleName               = "SAP",
                                                      CFBundleShortVersionString = version,     # must be in X.X.X format
                                                      CFBundleGetInfoString      = "SAP " + version,
                                                      CFBundleExecutable         = "SAP",
                                                      CFBundleIdentifier         = "com.example.myappname",
                                                      CFBundleDocumentTypes=[dict(CFBundleTypeExtensions=['.fasta'],
                                                                                  #dict(CFBundleTypeIconFile='doc.icns'),
                                                                                  CFBundleTypeName='Fasta file',
                                                                                  CFBundleTypeRole="Viewer"), ],
                                                      #LSPrefersPPC=True,
                                                      ),
                                        )
                         )

    data_files = []
elif sys.platform == 'win32':
    try:
        import py2exe
        # NB: Only do"python setup.py bdist_wininst --install-script winpostinstall.py" when using py2exe on gui
        extra_options = dict(windows=[{"script": guiscript }],
                             # Look here: http://wiki.wxpython.org/index.cgi/DistributingYourApplication
                             #app=[guiscript],
#                              scripts=['winpostinstall.py'],
                             #setup_requires=['py2exe'],
                             )
    except ImportError:
        extra_options = {}

#     import fnmatch
#     import os
#     rootPath = '/'
#     pattern = '(Microsoft.VC90.CRT.manifest)|(msvcm90.dll)|(msvcp90.dll)|(msvcr90.dll)'
#     for root, dirs, files in os.walk(rootPath):
#         for filename in fnmatch.filter(files, pattern):
#             data_file_list.append(os.path.join(root, filename))
    data_files = [("Microsoft.VC90.CRT", glob.glob('Microsoft.VC90.CRT/*'))]
#     assert len(data_file_list) == 4
#     print data_file_list
#     data_files = [("Microsoft.VC90.CRT", tuple(data_file_list))]
else:
    extra_options = {}
    data_files = []

EXTRA_COMPILE_ARGS = []
EXTRA_LINK_ARGS = []

# This is because EPD is built against the 10.5 SDK to run on all more recent versions of OSX. If you use the 
if sys.platform=='darwin':# and sys.executable != '/usr/bin/python' and 'anaconda' not in sys.executable:
   import platform
   v, _, _ = platform.mac_ver()
   osx_version = float('.'.join(v.split('.')[:2]))
#    if v >= 10.9:
#       EXTRA_LINK_ARGS.append('-lstdc++')#'-L/Developer/SDKs/MacOSX10.5.sdk/usr/lib']
   if v >= 10.10:
      EXTRA_COMPILE_ARGS.append("-stdlib=libc++")
#      EXTRA_COMPILE_ARGS.append("-mmacosx-version-min=10.14")
      os.environ["MACOSX_DEPLOYMENT_TARGET"] = platform.mac_ver()[0]
      
    #   EXTRA_COMPILE_ARGS.append("-stdlib=libstdc++")
# EXTRA_COMPILE_ARGS = ['-fvisibility=hidden']

setup(name='SAP',
      test_suite='tests',
      version=version,
      description='Statistical Assignment Package (SAP)',
      long_description='SAP does statistial assignment of unknown DNA to estabilish what taxononomic groups the DNA sample originates from. Itx uses a Baysian approach to calculate a probability distribution over all taxa represented in a sequence database. The probability of assignment to eaach taxa serves as a measure of confidence in the assignment.',
      author='Kasper Munch',
      author_email='kasmunch@bi.ku.dk',
      url='http://www.binf.ku.dk/~kasper/wiki/SAP.html',
      packages = find_packages(),
      package_dir = {'SAP': 'SAP'},
      include_package_data = True,      
      entry_points = { 'console_scripts': [ 'sap = SAP.ConsoleScripts:sap', ],
                       #'gui_scripts': [ 'sap_gui = SAP.GUI:start_gui', ],
                       'sap.database': [ 'Native = SAP.Databases.Native',
                                         'GenBank = SAP.Databases.GenBank' ],
                       'sap.alignment': [ 'Clustalw2 = SAP.Alignment.Clustalw2',
#                                           'MapToPreAligned = SAP.Alignment.MapToPreAligned',
                                          'PreAligned = SAP.Alignment.PreAligned' ],
                       'sap.assignment': [ 'Barcoder = SAP.Assignment.Barcoder',
                                           'ConstrainedNJ = SAP.Assignment.ConstrainedNJ' ],
                       },
      ext_modules=[Extension('SAP.Assignment.Barcoder.wrapBarcoder', # Has to be prefixed with the underscore because this is the name generated by swig.
                             glob.glob('ext/barcoder/*.cpp'),
                             include_dirs=['ext/barcoder'],
                             extra_link_args = EXTRA_LINK_ARGS,
                             export_symbols = ['runprogram'],
                             extra_compile_args = EXTRA_COMPILE_ARGS,
                             ),
                   Extension('SAP.Assignment.ConstrainedNJ.cConstrainedNJlib', # Has to be prefixed with the underscore because this is the name generated by swig.
#                             glob.glob('ext/constrnj/*.c[px][px]'),
                             glob.glob('ext/constrnj/*.cpp'),
                             include_dirs=['ext/constrnj'],
                             extra_link_args = EXTRA_LINK_ARGS,
                             extra_compile_args = EXTRA_COMPILE_ARGS,
                             ),
#                    Extension('SAP.Bio.Nexus.cnexus',
#                              ['SAP/Bio/Nexus/cnexus.c'],
#                              extra_link_args = EXTRA_LINK_ARGS,
#                              ),
                   Extension('SAP.PostAnalysis.IMa.wrapIMa', # Has to be prefixed with the underscore because this is the name generated by swig.
                             glob.glob('ext/ima/*.c') + glob.glob('ext/ima/ima/src/*.c'),
                             include_dirs=['ext/ima', 'ext/ima/ima', 'ext/ima/ima/src'],
                             extra_link_args = EXTRA_LINK_ARGS,
                             extra_compile_args = EXTRA_COMPILE_ARGS,
                             ),
                   #Extension('SAP.Bio.Nexus.cnexus',
                   #          ['SAP/Bio/Nexus/cnexus.c']
                   #          ),
                   ],
      data_files = data_files,
      **extra_options
      )
