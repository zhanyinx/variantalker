; MAFigate Windows Installer - Inno Setup Script
;
; Prerequisites to build this installer:
;   1. Install Inno Setup: https://jrsoftware.org/isinfo.php
;   2. Run build_installer.bat — NOT this script on its own. It does two things this
;      script cannot do for itself: it stages the bundled interpreter into `python\`
;      (compiling straight from the Inno Setup Compiler fails on the missing Source:
;      until that has been done once) and it derives the version (see below, where
;      compiling without that step stops with an #error).
;
; The installer will:
;   - Copy MAFigate app files, and a portable CPython, to the application directory
;   - Create Start Menu and Desktop shortcuts
;   - On first launch, the app auto-creates a Python venv and installs deps
;
; The user needs no Python of their own (issue #261). The release is pinned in
; build\python_release.env, shared with the macOS build.
;
; THE VERSION IS NOT WRITTEN IN THIS FILE (issue #260). It was, twice — AppVersion and
; OutputBaseFilename — and it had drifted from the app's own APP_VERSION by a major
; number while nothing reconciled them. build/version.py reads config/constants.py and
; writes version.iss, which is generated, untracked, and included below.
;
; There is exactly one way in, and deliberately so. An earlier draft of this file also
; honoured `ISCC /DAppVersion=1.2.3` for the release workflow's benefit, which is the same
; externally-typed override the Makefile's `VERSION ?=` was deleted for: a way to stamp an
; installer with a number the app does not agree with. The workflow can run
; `version.py --write-iss` exactly as build_installer.bat does, and then the constant wins
; on every route. Including unconditionally is part of that — the include is the authority,
; so a /D on the command line cannot outrank it.
;
; With no include, the compile stops. Inno Setup expands an undefined symbol to nothing
; rather than failing, so an unversioned installer is the silent alternative.
#ifexist "version.iss"
  #include "version.iss"
#endif
#ifndef AppVersion
  ; Kept short because an ISPP string cannot be wrapped; the comment above is the long form.
  #error "No AppVersion. Build via build_installer.bat (it runs build/version.py --write-iss)."
#endif

[Setup]
AppName=MAFigate
AppVersion={#AppVersion}
AppPublisher=VarianTalker
; Rendered in Windows' Add/Remove Programs, so it is a user-facing link rather than a
; comment. It carried a scaffolding placeholder where the account name belongs — broken
; from the day it was written, and invisible to the sweep in tests/test_public_repo_name.py
; because it named neither repository correctly. tests/test_unsigned_artifact_copy.py now
; reads this key and forbids the placeholder token tree-wide.
AppPublisherURL=https://github.com/zhanyinx/variantalker
DefaultDirName={autopf}\MAFigate
DefaultGroupName=MAFigate
OutputDir=output
OutputBaseFilename=MAFigate-{#AppVersion}-Windows-Setup
Compression=lzma2
SolidCompression=yes
; No icon.ico is tracked in this repository, and the [Files] and [Icons] entries below were
; already written to tolerate that — but this directive is read by the *compiler*, so an
; absent icon aborted the build rather than falling back to Inno's own (issue #265, tag
; mafigate-v1.0.0-rc2, the first Windows compile that ever got this far). #ifexist keeps the
; intent and makes the tolerance uniform: drop an icon.ico beside this file and all three
; start using it. UninstallDisplayIcon points into {app} rather than at a source file, so it
; is resolved at install time and needs no guard.
#ifexist "icon.ico"
SetupIconFile=icon.ico
#endif
UninstallDisplayIcon={app}\icon.ico
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
WizardStyle=modern
PrivilegesRequired=lowest
; Three levels up, not two: this file is build/windows/, so ..\..\ is streamlit_app/ and the
; LICENSE is at the repository root beside it. Wrong since it was written, and invisible until
; a compiler read it (issue #265).
LicenseFile=..\..\..\LICENSE

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
; Streamlit's own config, which is what turns its usage reporting off — and so what makes
; this installer's welcome screen true. It must land beside MAFigate.py: Streamlit
; resolves a `.streamlit` directory relative to the main script, and launch.bat runs from
; the installation root one level above the app (issue #259).
Source: "..\..\.streamlit\*"; DestDir: "{app}\streamlit_app\.streamlit"; Flags: ignoreversion recursesubdirs
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
; The bundled interpreter, staged here by build_installer.bat from the release pinned in
; build\python_release.env. launch.bat runs {app}\python\python.exe and looks nowhere else,
; so this line is what makes the app launchable on a machine with no Python at all.
; No skipifsourcedoesntexist: an .exe built without this is an .exe that cannot start.
Source: "python\*"; DestDir: "{app}\python"; Flags: ignoreversion recursesubdirs createallsubdirs
; Icon
Source: "icon.ico"; DestDir: "{app}"; Flags: ignoreversion skipifsourcedoesntexist

[Icons]
Name: "{group}\MAFigate"; Filename: "{app}\launch.bat"; IconFilename: "{app}\icon.ico"; Comment: "Launch MAFigate"
Name: "{group}\Uninstall MAFigate"; Filename: "{uninstallexe}"
Name: "{autodesktop}\MAFigate"; Filename: "{app}\launch.bat"; IconFilename: "{app}\icon.ico"; Tasks: desktopicon; Comment: "Launch MAFigate"

[Run]
Filename: "{app}\launch.bat"; Description: "Launch MAFigate now"; Flags: nowait postinstall skipifsilent shellexec

[UninstallDelete]
; Uninstalling reclaims what the app generated, and nothing the user made.
;
; This used to be one line naming `{userprofile}\.mafigate` - the whole directory. That is
; not just the venv: the app keeps the parameter cache and the file history in there, so
; uninstalling silently destroyed a clinician's saved filter setups, with no dialog saying
; so, and a reinstall came back empty (issue #261).
;
; The venv and the logs are bytes this app wrote and can write again, so they go. Anything
; else in that directory stays, which is also why these name the two subdirectories rather
; than listing what to keep: something the app starts keeping between sessions tomorrow
; survives an uninstall without anyone remembering to add it here.
;
; {%USERPROFILE}, not {userprofile}: Inno has no constant of that name and aborts the
; compile with `Unknown constant "userprofile"`. {%NAME} is how it reads an environment
; variable, and %USERPROFILE% is what launch.bat builds these same two paths from, so the
; uninstaller and the launcher now agree by construction. The line above records the
; directory the old rule named and keeps its original spelling, which is why it still
; reads {userprofile} (issue #265, tag mafigate-v1.0.0-rc3 - the first compile to reach
; this section).
Type: filesandordirs; Name: "{%USERPROFILE}\.mafigate\venv"
Type: filesandordirs; Name: "{%USERPROFILE}\.mafigate\logs"

[Messages]
WelcomeLabel2=This will install MAFigate on your computer.%n%nMAFigate is a local application for exploring and annotating MAF variant files. All data stays on your computer - nothing is uploaded to external servers.%n%nEverything MAFigate needs is included - there is nothing to install first.%n%nThe first launch needs an internet connection while it builds its environment.
