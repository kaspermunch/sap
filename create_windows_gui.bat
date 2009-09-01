REM Install python 2.5 if you don't have it. I highly recommend the Enthought Python Distribution (EPD) 
REM that comes with a lot of things that makes compiling SAP a breeze.
REM 
REM Add ";C:\Python25;" (without the quotes) to the enf of your "path" environmental variable
REM 
REM Install ansi wxpython for 2.5
REM 
REM Install py2exe for 2.5
REM 
REM If have mingw installed (you do if you are using the Enthought Python Distribution) then just skip to where it says INSTALL USING MINGW COMPILER.
REM If not, read below:
REM
REM INSTALL USING BORLAND COMPILER:
REM
REM Install Borland compiler 5.5 and make the necessary changes to Path and environmental
REM variables (form the supplemental information on the website):
REM 
REM  Configuring the system environment:
REM  
REM      Open a console box.
REM      1.  Start | Run...
REM      2.  Type "command" into the field [Enter]
REM  
REM      * If Windows 95/98:
REM      Navigate to the root in order to modify the PATH reference in the autoexec.bat file.
REM      3.  Type "cd" [Enter]
REM      4.  Type "edit autoexec.bat" [Enter]
REM      5.  Insert a line and type "PATH=C:\BORLAND\BCC55\BIN;%PATH%"
REM      6.  Save the changes (Alt-F then hit S).
REM      7.  Exit edit. (Alt+F then press X).
REM  
REM      * If Windows NT:
REM      Add a path reference to the Environment variables:
REM      3.  Using the mouse, right-click on the "My Computer" icon (on your desktop) and choose "Properties".
REM      4.  Click on the "Environment" tab.
REM      5.  Click on "Path" in the "System Variables" field.
REM      6.  Highlight the "Path" System variable (bottom).
REM      7.  Click in the "Value" field.
REM      8.  Append the line with ";C:\BORLAND\BCC55\BIN;" (exactly 1 semi-colon between references)
REM      9.  Click on the "Set" button.
REM      10.  Click OK (in the "System Properties" window)
REM  
REM      * Or, if Windows 2000/XP:
REM      Add a path reference to the Environment variables:
REM      3.  Using the mouse, right-click on the "My Computer" icon (on your desktop) and choose "Properties".
REM      4.  Click on the "Advanced" tab.
REM      5.  Click on the "Environment Variables..." button.
REM      6.  Highlight the "Path" System variable (bottom).
REM      7.  Click on the "Edit..." button.
REM      8.  Append the line with ";C:\BORLAND\BCC55\BIN;"
REM      9.  Click OK (in the "Edit System Variables")
REM      10. Click OK (in the "Environment Variables" window) and click OK (in the "System Properties" window) Navigating to the directory, "c:\Borland\bcc55\bin"
REM      11. cd borland [Enter]
REM      12. cd bcc55 [Enter]
REM      13. cd bin [Enter]
REM  
REM  ______________________________
REM  
REM  Creating the configuration files:
REM  
REM      Note: The command line should read:  C:\BORLAND\BCC55\BIN
REM  
REM      Part 1: Creating BCC32.CFG.
REM      1.  Type "edit bcc32.cfg" [Enter]  (This creates the file and opens a blank window in the editor).
REM      2. Add these lines:
REM  
REM          -I"c:\Borland\Bcc55\include"
REM          -L"c:\Borland\Bcc55\lib"
REM  
REM      3. Save the changes (Alt-F then hit S).
REM      4. Exit edit. (Alt+F then press X).
REM  
REM      Part 2: Creating ILINK32.CFG
REM      5. Type "edit ilink32.cfg"  (This creates the file and opens a blank window in the editor).
REM      6. Add these lines:
REM  
REM          -L"c:\Borland\Bcc55\lib"
REM  
REM      7. Save the changes (Alt-F then hit S).
REM      8. Exit edit. (Alt+F then press X).
REM      9. Type "exit" [Enter]
REM      10. Restart Windows.
REM   
REM You also have to convert Python's library python24.lib into the Borland format. You can do this as follows:
REM  
REM cd C:\Python25\libs
REM coff2omf python24.lib python24_bcpp.lib
REM 
REM The converted files have to reside in the same directories as the normal libraries.
REM 
REM How does Distutils manage to use these libraries with their changed names? If the
REM extension needs a library (eg. foo) Distutils checks first if it finds a library with
REM suffix _bcpp (eg. foo_bcpp.lib) and then uses this library. In the case it doesn't find
REM such a special library it uses the default name (foo.lib.)
REM 
REM To install the command line version of SAP do:
REM 
REM python setup.py setopt --command=build --option=compiler --set-value=bcpp
REM python setup.py install
REM 
REM Before you can run sap from the command line you need to do the follwoing
REM
REM INSTALL USING MINGW COMPILER:
REM
REM python setup.py build -c mingw32 install
REM
REM AFTER YOU INSTALL:
REM 
REM   Add a path reference to the Environment variables:
REM   3.  Using the mouse, right-click on the "My Computer" and choose "Properties".
REM   4.  Click on the "Advanced" tab or "Advanced system settings".
REM   5.  Click on the "Environment Variables..." button.
REM   6.  Highlight the "Path" System variable (bottom).
REM   7.  Click on the "Edit..." button.
REM   8.  Append the line with ";C:Python25\Scripts;"
REM   9.  Click OK on your way out.
REM 
REM Locate "Set Associations" using the search field in "Control Panels". Hithlight the ".py"
REM file extension and associate it with the python executable in C:\Python25. Do the same for
REM the ".pyc" extension.

REM If you want to create a windows distribution for an unsupported version of Windows (to make SAP available to more users) install NSIS.
REM You can then use the following commands:
REM 
REM Using mingw compiler:       
REM python setup.py setopt --command=build --option=compiler --set-value=mingw32
REM
REM Using borland compiler:       
REM python setup.py setopt --command=build --option=compiler --set-value=bcpp
REM 
REM python setup.py bdist_wininst
REM 
REM REM python setup.py bdist_msi
REM 
REM python setup.py py2exe
REM 
REM copy C:\Python25\lib\site-packages\wx-2.8-msw-ansi\wx\gdiplus.dll dist
REM copy C:\Python25\lib\site-packages\wx-2.8-msw-ansi\wx\MSVCP71.dll dist
REM 
REM "C:\Program Files\NSIS\makensis.exe" win_installer_script.nsi

python setup.py setopt --command=build --option=compiler --set-value=mingw32

python setup.py bdist_wininst

REM python setup.py bdist_msi

python setup.py py2exe

copy C:\Python25\lib\site-packages\wx-2.8-msw-ansi\wx\gdiplus.dll dist
copy C:\Python25\lib\site-packages\wx-2.8-msw-ansi\wx\MSVCP71.dll dist

"C:\Program Files\NSIS\makensis.exe" win_installer_script.nsi
