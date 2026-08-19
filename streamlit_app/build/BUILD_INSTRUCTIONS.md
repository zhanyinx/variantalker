# MAFigate Desktop Installer Build Guide

Build standalone `.dmg` (macOS) and `.exe` (Windows) installers for MAFigate.

**Key principle**: MAFigate runs 100% locally — no data leaves the user's machine.

---

## macOS (.dmg)

The DMG bundles a portable Python 3.11 — **users need zero prerequisites** (no Homebrew, no Python install). Works on both Apple Silicon (M1/M2/M3/M4) and Intel Macs.

### Prerequisites (build machine only)
- macOS with Xcode Command Line Tools (`xcode-select --install`)
- Internet connection (to download Python on first build; cached for subsequent builds)
- Optional: `brew install create-dmg` for a polished DMG with drag-to-Applications layout

### Build

```bash
cd streamlit_app/build/mac

# Build a universal DMG (Intel + Apple Silicon) — default, recommended for sharing
./build_dmg.sh 1.0.0

# Build for a specific architecture only
./build_dmg.sh 1.0.0 --arch arm64      # Apple Silicon only (~33 MB)
./build_dmg.sh 1.0.0 --arch x86_64     # Intel only (~35 MB)
./build_dmg.sh 1.0.0 --arch universal  # Both architectures (~65 MB, == default)
```

> The default is now **universal** so the DMG runs on both Apple Silicon and
> Intel Macs; a single-arch build crashes on the other architecture when shared.
> The build also compiles a small native launcher (via `swiftc`, included with
> the Xcode Command Line Tools) so the app stays responsive and quits cleanly
> (Cmd-Q / Dock → Quit) instead of appearing as "not responding".

### Output
`MAFigate-1.0.0-macOS-universal.dmg` (or `-arm64` / `-x86_64`)

### What the user does
1. Open the `.dmg` file
2. Drag `MAFigate.app` to Applications
3. Double-click MAFigate from Applications
4. First launch: automatically installs pip dependencies (~1 min, needs internet)
5. Browser opens with MAFigate at `http://127.0.0.1:8501`

**No Python installation required.** Python 3.11 is bundled inside the app.

### App icon (optional)
Place an `icon.icns` file at `build/mac/MAFigate.app/Contents/Resources/icon.icns`.
Generate from PNG: `sips -s format icns icon.png --out icon.icns`

---

## Windows (.exe installer)

### Prerequisites
- Windows PC (cross-compilation from Mac is not supported)
- [Inno Setup 6](https://jrsoftware.org/isinfo.php) installed
- **Note**: Windows installer currently requires users to have Python 3.9+ installed

### Build

```cmd
cd streamlit_app\build\windows
build_installer.bat
```

### Output
`output\MAFigate-1.0.0-Windows-Setup.exe`

### What the user does
1. Install Python 3.9+ from https://www.python.org/downloads/ (check "Add Python to PATH")
2. Run `MAFigate-1.0.0-Windows-Setup.exe`
3. Follow the install wizard
4. Launch from Start Menu or Desktop shortcut
5. First launch: automatically installs Python dependencies (~1-2 min)
6. Browser opens with MAFigate at `http://127.0.0.1:8501`

### App icon (optional)
Place an `icon.ico` file at `build/windows/icon.ico`.

---

## How it works

### macOS
1. App bundle includes a full portable Python 3.11 (~33 MB)
2. On first launch, creates a virtual environment in `~/.mafigate/venv`
3. Installs pip dependencies from `requirements.txt`
4. Starts Streamlit on `127.0.0.1` (localhost only — not network-accessible)
5. Opens the default browser
6. Subsequent launches skip steps 2-3 and start in seconds

### Windows
1. Finds Python 3.9+ on the system PATH
2. Same steps 2-6 as macOS (venv at `%USERPROFILE%\.mafigate\venv`)

---

## File structure

```
build/
├── BUILD_INSTRUCTIONS.md     ← this file
├── launcher.py               ← cross-platform Python launcher (alternative)
├── mac/
│   ├── build_dmg.sh          ← builds the .dmg (downloads + bundles Python)
│   ├── .python_cache/        ← cached Python downloads (auto-created)
│   └── MAFigate.app/         ← .app bundle template
│       └── Contents/
│           ├── Info.plist
│           ├── MacOS/
│           │   └── launch.sh
│           └── Resources/
│               └── (icon.icns — add your own)
└── windows/
    ├── build_installer.bat    ← builds the .exe installer
    ├── installer.iss          ← Inno Setup script
    └── launch.bat             ← Windows launcher
```

---

## Distributing to users

Share the `.dmg` or `.exe` file via:
- USB drive
- Internal file share / NAS
- Institutional cloud storage (OneDrive, Google Drive, etc.)
- Email (if file size allows — ~33 MB for macOS)

Internet is needed only on first launch to install pip dependencies.

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| macOS "app is damaged" warning | Run: `xattr -cr /Applications/MAFigate.app` |
| macOS Gatekeeper blocks app | System Settings → Privacy & Security → Open Anyway |
| Windows "Python not found" | Install Python 3.9+ with "Add to PATH" checked |
| Windows SmartScreen warning | Click "More info" → "Run anyway" |
| Dependencies fail to install | Check internet connection; check `~/.mafigate/logs/` |
| Port 8501 in use | App automatically finds the next free port |
| App update: stale venv | Delete `~/.mafigate/venv` and relaunch |
