!include "MUI2.nsh"
Name "modhmm-${PACKAGE_VERSION}"
outfile "${CMAKE_CURRENT_BINARY_DIR}/modhmm-${PACKAGE_VERSION}-win32.exe"
installDir $PROGRAMFILES\modhmm

!define MUI_ABORTWARNING

!insertmacro MUI_PAGE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/license-win32-statically-built"
; !insertmacro MUI_PAGE_COMPONENTS 
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
 
!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES
  
!insertmacro MUI_LANGUAGE "English"

Section "Dummy Section" SecDummy
  SetOutPath "$INSTDIR"
  file modhmmt.exe
  file modhmms.exe
  file modhmmseqalign.exe
  file add_alphabet.exe
  WriteUninstaller "$INSTDIR\Uninstall.exe"
SectionEnd

Section "Uninstall"
  Delete "$INSTDIR\Uninstall.exe"
  Delete "$INSTDIR\modhmmt.exe"
  Delete "$INSTDIR\modhmms.exe"
  Delete "$INSTDIR\modhmmseqalign.exe"
  Delete "$INSTDIR\add_alphabet.exe"
  RMDir "$INSTDIR"
SectionEnd