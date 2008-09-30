import sys, os

if sys.argv[1] == '-install':

    try:
        program_folder = get_special_folder_path(CSIDL_PROGRAMS)        
    except OSError:
        pass

    try:
        start_menu_folder = get_special_folder_path(CSIDL_STARTMENU)
    except OSError:
        pass

    #########################################################
    for s in ["CSIDL_APPDATA", "CSIDL_COMMON_STARTMENU",  "CSIDL_STARTMENU",
              "CSIDL_COMMON_DESKTOPDIRECTORY", "CSIDL_DESKTOPDIRECTORY",
              "CSIDL_COMMON_STARTUP", "CSIDL_STARTUP", "CSIDL_COMMON_PROGRAMS",
              "CSIDL_PROGRAMS",  "CSIDL_FONTS"]:
        try:
            path = get_special_folder_path(s)
            print "%s: %s", (s, path)
        except OSError:
            print "Could not find", s
    
    # Do we find our executable here?
    import glob
    print glob.glob(os.path.join(program_folder, 'SAP', '*'))
    print tglob.glob(os.path.join(program_folder, '*'))
    #########################################################


    # Make a dir for the shortcut:
    shortcut_folder = os.path.join(start_menu_folder, 'SAP')

    print "would create dir:", shortcut_folder
    #os.mkdirs(os.path.join(start_menu_folder, 'SAP'))
    #directory_created(shortcut_folder)

    # make the shortcut:
    shortcut = os.path.join(shortcut_folder, 'SAP', 'SAP.lnk')
    description = "..."
    target = os.path.join(program_folder, 'sap.exe')    
    print "would create shortcut", target, description, shortcut
    #create_shortcut(target, description, shortcut)
    #file_created(shortcut)

if sys.argv[1] == '-remove':
    pass

else:
    pass


create_shortcut( target, description, filename[, arguments[, workdir[, iconpath[, iconindex]]]]) 


The first link try to launch python with a script as argument, so values are :

target=os.path.join(sys.exec_prefix,"python.exe")
description="..."
shortcut=os.path.join(some_dir,"mylink.lnk")
arguments="path_to_myscript.py"

and finally 
create_shortcut(target,description,shortcut,arguments)

This creates a link but when it is activated, it tells 
"C:\Python23\python.exe is not a valid Win32 application."
This sounds like the problem is coming from the file C:\Python23\python.exe and not from the shortcut I just create. Is something wrong with the way I make the link ?


My second link starts Acrobat reader with a file as argument. This time, values are :

target=os.path.join(
       "C:\\Program Files\\Adobe\\Acrobat 5.\\Acrobat",
       "Acrobat.exe") 
description="..."
shortcut=os.path.join(some_dir,"Mydocumentpdf.lnk")
arguments="path_to_doc.pdf"

create_shortcut(target,description,shortcut,arguments)


# 
# const wincerapi.CSIDL_STARTMENU;
# 
# File system directory containing Start menu items.
# 
# 
# CSIDL_PROGRAMS (FOLDERID_Programs)
# 
#     The file system directory that contains the users program groups (which are themselves file system directories). A typical path is C:\Documents and Settings\username\Start Menu\Programs.
# 
# 
#     Which folders are available depends on the exact Windows version, and probably also the configuration. For details refer to Microsoft's documentation of the SHGetSpecialFolderPath() function. 
# 
# create_shortcut( 	target, description, filename[, arguments[, workdir[, iconpath[, iconindex]]]])
#     This function creates a shortcut. target is the path to the program to be started by the shortcut. description is the description of the shortcut. filename is the title of the shortcut that the user will see. arguments specifies the command line arguments, if any. workdir is the working directory for the program. iconpath is the file containing the icon for the shortcut, and iconindex is the index of the icon in the file iconpath. Again, for details consult the Microsoft documentation for the IShellLink interface.
# 
#     
#     print "installing"
