import sys, os

if sys.argv[1] == '-install':

    try:
        start_menu_programs_folder = get_special_folder_path('CSIDL_PROGRAMS')
    except OSError:
        pass

    # Make a dir for the shortcut:
    shortcut_folder = os.path.join(start_menu_programs_folder, 'SAP')

    if not os.path.exists(shortcut_folder):
        os.makedirs(os.path.join(start_menu_programs_folder, 'SAP'))
    directory_created(shortcut_folder)

    # make the shortcut:
    shortcut = os.path.join(shortcut_folder, 'SAP.lnk')
    description = "..."
    target = os.path.join(sys.prefix, 'Scripts', 'sap_gui.exe')    
    if os.path.exists(shortcut):
        os.unlink(shortcut)
    create_shortcut(target, "Statistical Assignmnet Package", shortcut)
    file_created(shortcut)

if sys.argv[1] == '-remove':
    pass

else:
    pass


