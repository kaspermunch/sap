Page license
Page instfiles
UninstPage uninstConfirm
UninstPage instfiles

Name SAP
LicenseText "SAP is distributed under the GNU Public License" "Install SAP"


!define UninstLog "uninstall.log"
Var UninstLog
 
; Uninstall log file missing.
LangString UninstLogMissing ${LANG_ENGLISH} "${UninstLog} not found!$\r$\nUninstallation cannot proceed!"
 
; AddItem macro
!macro AddItem Path
 FileWrite $UninstLog "${Path}$\r$\n"
!macroend
!define AddItem "!insertmacro AddItem"
 
; File macro
!macro File FilePath FileName
 IfFileExists "$OUTDIR\${FileName}" +2
  FileWrite $UninstLog "$OUTDIR\${FileName}$\r$\n"
 File "${FilePath}${FileName}"
!macroend
!define File "!insertmacro File"
 
; Copy files macro
!macro CopyFiles SourcePath DestPath
 IfFileExists "${DestPath}" +2
  FileWrite $UninstLog "${DestPath}$\r$\n"
 CopyFiles "${SourcePath}" "${DestPath}"
!macroend
!define CopyFiles "!insertmacro CopyFiles"
 
; Rename macro
!macro Rename SourcePath DestPath
 IfFileExists "${DestPath}" +2
  FileWrite $UninstLog "${DestPath}$\r$\n"
 Rename "${SourcePath}" "${DestPath}"
!macroend
!define Rename "!insertmacro Rename"
 
; CreateDirectory macro
!macro CreateDirectory Path
 CreateDirectory "${Path}"
 FileWrite $UninstLog "${Path}$\r$\n"
!macroend
!define CreateDirectory "!insertmacro CreateDirectory"
 
; SetOutPath macro
!macro SetOutPath Path
 SetOutPath "${Path}"
 FileWrite $UninstLog "${Path}$\r$\n"
!macroend
!define SetOutPath "!insertmacro SetOutPath"
 
; WriteUninstaller macro
!macro WriteUninstaller Path
 WriteUninstaller "${Path}"
 FileWrite $UninstLog "${Path}$\r$\n"
!macroend
!define WriteUninstaller "!insertmacro WriteUninstaller"
 
Section -openlogfile
 CreateDirectory "$INSTDIR"
 IfFileExists "$INSTDIR\${UninstLog}" +3
  FileOpen $UninstLog "$INSTDIR\${UninstLog}" w
 Goto +4
  SetFileAttributes "$INSTDIR\${UninstLog}" NORMAL
  FileOpen $UninstLog "$INSTDIR\${UninstLog}" a
  FileSeek $UninstLog 0 END
SectionEnd


############################################################################

# define name of installer
outFile ".\dist\SAP installer.exe"

# define installation directory
#installDir $DESKTOP\test_installation
installDir "C:\Program Files\SAP"

Section "Install Main"
SectionIn RO
 
 ${SetOutPath} $INSTDIR
 ${WriteUninstaller} "$INSTDIR\uninstall SAP.exe"

 ${File} "dist\" "GUI.exe" 
 ${File} "dist\" "MSVCR71.dll"
 ${File} "dist\" "SAP.Assignment.Barcoder._Barcoder.pyd"
 ${File} "dist\" "SAP.Assignment.ConstrainedNJ._cConstrainedNJlib.pyd"
 ${File} "dist\" "SAP.Bio.Nexus.cnexus.pyd"
 ${File} "dist\" "_hashlib.pyd"
 ${File} "dist\" "_socket.pyd"
 ${File} "dist\" "_ssl.pyd"
 ${File} "dist\" "bz2.pyd"
 ${File} "dist\" "gdiplus.dll"
 ${File} "dist\" "library.zip"
 ${File} "dist\" "msvcp71.dll"
 ${File} "dist\" "pyexpat.pyd"
 ${File} "dist\" "python25.dll"
 ${File} "dist\" "select.pyd"
 ${File} "dist\" "unicodedata.pyd"
 ${File} "dist\" "w9xpopen.exe"
 ${File} "dist\" "wx._controls_.pyd"
 ${File} "dist\" "wx._core_.pyd"
 ${File} "dist\" "wx._gdi_.pyd"
 ${File} "dist\" "wx._misc_.pyd"
 ${File} "dist\" "wx._windows_.pyd"
 ${File} "dist\" "wxbase28uh_net_vc.dll"
 ${File} "dist\" "wxbase28uh_vc.dll"
 ${File} "dist\" "wxmsw28uh_adv_vc.dll"
 ${File} "dist\" "wxmsw28uh_core_vc.dll"
 ${File} "dist\" "wxmsw28uh_html_vc.dll"

 # create a shortcut named "new shortcut" in the start menu programs directory
 # point the new shortcut at the program uninstaller
 SetShellVarContext all
 createShortCut "$SMPROGRAMS\SAP.lnk" "$INSTDIR\GUI.exe"


SectionEnd

############################################################################

Section -closelogfile
 FileClose $UninstLog
 SetFileAttributes "$INSTDIR\${UninstLog}" READONLY|SYSTEM|HIDDEN
SectionEnd
 
Section Uninstall
 
 ; Can't uninstall if uninstall log is missing!
 IfFileExists "$INSTDIR\${UninstLog}" +3
  MessageBox MB_OK|MB_ICONSTOP "$(UninstLogMissing)"
   Abort
 
 Push $R0
 Push $R1
 Push $R2
 SetFileAttributes "$INSTDIR\${UninstLog}" NORMAL
 FileOpen $UninstLog "$INSTDIR\${UninstLog}" r
 StrCpy $R1 0
 
 GetLineCount:
  ClearErrors
   FileRead $UninstLog $R0
   IntOp $R1 $R1 + 1
   IfErrors 0 GetLineCount
 
 LoopRead:
  FileSeek $UninstLog 0 SET
  StrCpy $R2 0
  FindLine:
   FileRead $UninstLog $R0
   IntOp $R2 $R2 + 1
   StrCmp $R1 $R2 0 FindLine
 
   StrCpy $R0 $R0 -2
   IfFileExists "$R0\*.*" 0 +3
    RMDir $R0  #is dir
   Goto +3
   IfFileExists $R0 0 +2
    Delete $R0 #is file
 
  IntOp $R1 $R1 - 1
  StrCmp $R1 0 LoopDone
  Goto LoopRead
 LoopDone:
 FileClose $UninstLog
 Delete "$INSTDIR\${UninstLog}"
 RMDir "$INSTDIR"
 Pop $R2
 Pop $R1
 Pop $R0

 SetShellVarContext all
 delete "$SMPROGRAMS\SAP.lnk"

SectionEnd
