; MAFigate Windows Installer - Inno Setup Script
;
; Prerequisites to build this installer:
;   1. Install Inno Setup: https://jrsoftware.org/isinfo.php
;   2. Run this script from Inno Setup Compiler
;
; The installer will:
;   - Copy MAFigate app files to Program Files
;   - Create Start Menu and Desktop shortcuts
;   - On first launch, the app auto-creates a Python venv and installs deps

[Setup]
AppName=MAFigate
AppVersion=1.0.0
AppPublisher=VarianTalker
AppPublisherURL=https://github.com/your-org/variantalker
DefaultDirName={autopf}\MAFigate
DefaultGroupName=MAFigate
OutputDir=output
OutputBaseFilename=MAFigate-1.0.0-Windows-Setup
Compression=lzma2
SolidCompression=yes
SetupIconFile=icon.ico
UninstallDisplayIcon={app}\icon.ico
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
WizardStyle=modern
PrivilegesRequired=lowest
LicenseFile=..\..\LICENSE

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
; Main launcher
Source: "launch.bat"; DestDir: "{app}"; Flags: ignoreversion
; Streamlit app source
Source: "..\..\MAFigate.py"; DestDir: "{app}\streamlit_app"; Flags: ignoreversion
Source: "..\..\requirements.txt"; DestDir: "{app}\streamlit_app"; Flags: ignoreversion
Source: "..\..\components\*"; DestDir: "{app}\streamlit_app\components"; Flags: ignoreversion recursesubdirs
Source: "..\..\config\*"; DestDir: "{app}\streamlit_app\config"; Flags: ignoreversion recursesubdirs
Source: "..\..\filters\*"; DestDir: "{app}\streamlit_app\filters"; Flags: ignoreversion recursesubdirs
Source: "..\..\page_modules\*"; DestDir: "{app}\streamlit_app\page_modules"; Flags: ignoreversion recursesubdirs
Source: "..\..\utils\*"; DestDir: "{app}\streamlit_app\utils"; Flags: ignoreversion recursesubdirs
; Vendored copies of the pipeline's filter code (see vendor\README.md). Without this
; the installed app has no filter code at all. _sync.py is the developer-side repair
; tool and needs bin\, which the installed app does not carry, so it is excluded.
; The drift guard itself lives in tests\, which this script never copies.
Source: "..\..\vendor\*"; DestDir: "{app}\streamlit_app\vendor"; Flags: ignoreversion recursesubdirs; Excludes: "_sync.py"
Source: "..\..\resources\*"; DestDir: "{app}\streamlit_app\resources"; Flags: ignoreversion recursesubdirs skipifsourcedoesntexist
; Icon
Source: "icon.ico"; DestDir: "{app}"; Flags: ignoreversion skipifsourcedoesntexist

[Icons]
Name: "{group}\MAFigate"; Filename: "{app}\launch.bat"; IconFilename: "{app}\icon.ico"; Comment: "Launch MAFigate"
Name: "{group}\Uninstall MAFigate"; Filename: "{uninstallexe}"
Name: "{autodesktop}\MAFigate"; Filename: "{app}\launch.bat"; IconFilename: "{app}\icon.ico"; Tasks: desktopicon; Comment: "Launch MAFigate"

[Run]
Filename: "{app}\launch.bat"; Description: "Launch MAFigate now"; Flags: nowait postinstall skipifsilent shellexec

[UninstallDelete]
Type: filesandordirs; Name: "{userprofile}\.mafigate"

[Messages]
WelcomeLabel2=This will install MAFigate on your computer.%n%nMAFigate is a local application for exploring and annotating MAF variant files. All data stays on your computer - nothing is uploaded to external servers.%n%nRequirement: Python 3.9 or later must be installed.%nDownload from: https://www.python.org/downloads/

[Code]
// Check if Python is available before installation
function InitializeSetup(): Boolean;
var
  ResultCode: Integer;
begin
  Result := True;
  if not Exec('python', '--version', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then
  begin
    if not Exec('python3', '--version', '', SW_HIDE, ewWaitUntilTerminated, ResultCode) then
    begin
      if MsgBox('Python 3 was not detected on your system.' + #13#10 +
                'MAFigate requires Python 3.9 or later.' + #13#10#13#10 +
                'Would you like to continue anyway?' + #13#10 +
                '(You can install Python later from python.org)',
                mbConfirmation, MB_YESNO) = IDNO then
        Result := False;
    end;
  end;
end;
