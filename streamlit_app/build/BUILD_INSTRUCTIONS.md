# MAFigate Desktop Installer Build Guide

Build standalone `.dmg` (macOS) and `.exe` (Windows) installers for MAFigate.

**Key principle**: MAFigate runs 100% locally — no data leaves the user's machine.

**For a release, do not build by hand.** Push a tag and CI builds both, on GitHub-hosted
macOS and Windows runners, from the one commit — see [Cutting a release](#cutting-a-release).
The instructions below are for building locally: a change to the packaging you want to try,
or a one-off artifact for a single machine. Nothing on this page is a prerequisite for a
release any more, which is the point of it.

## The version is not typed anywhere

It lives in `streamlit_app/config/constants.py` as `APP_VERSION`, and every artifact takes
it from there: the DMG's filename, the `.app` bundle's two `Info.plist` version keys, the
Windows installer's version and output filename, and the tag a release is cut from. Nothing
below asks you for a version, and neither build script accepts one (issue #260).

```bash
cd streamlit_app
make version        # what this tree builds as
make release-tag    # what to hand `git tag` — mafigate-v<version>
```

Bumping a release is therefore a one-line edit to `APP_VERSION`.
`streamlit_app/tests/test_installer_version.py` fails if any installer input grows a
literal of its own, which is the state this tree was in before #260: the app said one
number and every installer said another.

## Each build stamps which build it is

A version number cannot identify a build: the same `APP_VERSION` reaches a user as a `.dmg`,
as a Windows `.exe`, or as a clone of the source. So each build script writes a **build
stamp** before it packages anything — its channel, and the short commit it was cut from — and
the app's About dialog reports both beside the version (issue #263). With **no update check
ever** (#229), that dialog is the only thing that tells a user, and therefore you reading
their bug report, what they are running.

```bash
python3 build/build_stamp.py --channel macos-dmg          # what build_dmg.sh runs
python3 build/build_stamp.py --channel windows-installer  # what build_installer.bat runs
```

Both write `streamlit_app/config/build_stamp.py`, which is **gitignored and never
committed** — its absence is exactly how the app knows it is running from a source checkout,
so a tracked copy would make every clone claim to be an installer build. You do not need to
run either command by hand; both build scripts do it, and neither takes a channel from you.

A build leaves that stamp behind in your working tree, where it outlives the build that wrote
it: run the app from source afterwards and About will call it a `.dmg`. `make clean` removes
it, which restores the source-checkout answer.

Unlike the version, a missing `git` does not stop a build: the stamp is written with an empty
commit, the build warns, and About names the channel and the version without a commit. That
is deliberate — someone building from a downloaded archive has no `.git` to read, and an
installer that refuses to exist is worse than one that cannot name its commit. A build from a
tree with uncommitted changes to tracked files is stamped `-dirty`, which is the honest answer
to *can I check this out?*

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
./build_dmg.sh

# Build for a specific architecture only
./build_dmg.sh --arch arm64      # Apple Silicon only (~33 MB)
./build_dmg.sh --arch x86_64     # Intel only (~35 MB)
./build_dmg.sh --arch universal  # Both architectures (~65 MB, == default)
```

> There is no version argument. Passing one — how this was called before #260 — now fails
> rather than being ignored, wherever in the command line it appears, because an ignored
> positional would have built a DMG labelled with a number nobody chose.

> The default is now **universal** so the DMG runs on both Apple Silicon and
> Intel Macs; a single-arch build crashes on the other architecture when shared.
> The build also compiles a small native launcher (via `swiftc`, included with
> the Xcode Command Line Tools) so the app stays responsive and quits cleanly
> (Cmd-Q / Dock → Quit) instead of appearing as "not responding".

### Output
`MAFigate-<version>-macOS-universal.dmg` (or `-arm64` / `-x86_64`), where `<version>` is
what `make version` prints

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

The installer bundles the same portable Python the DMG does — **users need zero
prerequisites**. x64 only; ARM Windows runs it under emulation.

### Prerequisites (build machine only)
- Windows PC (cross-compilation from Mac is not supported — but a release does not need one,
  see [Cutting a release](#cutting-a-release): `windows-latest` is the machine)
- [Inno Setup 6](https://jrsoftware.org/isinfo.php) installed
- Python 3 on `PATH` on the **build** machine, only so the version and the build stamp can be
  written
- `curl` and `tar`, both shipped with Windows 10 1803 and later
- Optional: `git` on `PATH`, only so the stamp can name the commit. Without it the build
  succeeds and About reports the channel and version with no commit
- Internet connection (to download Python on first build; cached for subsequent builds)

### Build

```cmd
cd streamlit_app\build\windows
build_installer.bat
```

> Run `build_installer.bat`, not `installer.iss` directly. It does two things Inno Setup
> cannot do for itself, and the compile fails without either: it downloads and extracts the
> interpreter into `windows\python\` (a missing `Source:` otherwise), and it generates
> `version.iss` from `APP_VERSION` (a compile with no version supplied stops with an
> `#error` rather than producing an unversioned installer). If you do want the IDE, run
> `python ..\version.py --write-iss` and one full batch build beforehand.

### Output
`output\MAFigate-<version>-Windows-Setup.exe`

### What the user does
1. Run `MAFigate-<version>-Windows-Setup.exe`
2. Follow the install wizard
3. Launch from Start Menu or Desktop shortcut
4. First launch: automatically installs Python dependencies (~1-2 min, needs internet)
5. Browser opens with MAFigate at `http://127.0.0.1:8501`

**No Python installation required.** Python 3.11 is bundled inside the installer.

### Uninstalling
Uninstall removes the app and the environment it built (`~\.mafigate\venv`) and its logs.
It **keeps** the parameter cache, also in `~\.mafigate`, because that is the user's own
work — so a reinstall opens with their saved filter setups intact. To
start genuinely clean, delete `%USERPROFILE%\.mafigate` by hand.

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
1. The installed app includes a portable Python 3.11 at `<install dir>\python` (~56 MB on
   disk, after the same trim the DMG applies — 73 MB untrimmed), and `launch.bat` runs that
   one only: it does not look at `PATH` or at any Python the machine happens to have
2. Same steps 2-6 as macOS (venv at `%USERPROFILE%\.mafigate\venv`)

Both platforms bundle the **same** release, pinned once in `build/python_release.env`. Move
the pin there and both builds follow; there is no second place to edit.

---

## File structure

```
build/
├── BUILD_INSTRUCTIONS.md     ← this file
├── OPENING_MAFIGATE.md       ← what recipients see when the alert fires, said once
├── RELEASE_PAGE.md           ← the release page's copy; CI fills its two placeholders
├── RELEASES.md               ← which releases have actually been published; a guard reads it
├── version.py                ← reads APP_VERSION; the only route from the constant to a build
├── build_stamp.py            ← writes ../config/build_stamp.py (channel + commit), gitignored
├── python_release.env        ← the bundled CPython release, pinned once for both builds
├── mac/
│   ├── build_dmg.sh          ← builds the .dmg (downloads + bundles Python)
│   ├── MAFigateLauncher.swift ← native launcher compiled into the .app
│   ├── .python_cache/        ← cached Python downloads (auto-created)
│   └── MAFigate.app/         ← .app bundle template
│       └── Contents/
│           ├── Info.plist    ← version keys are placeholders; build_dmg.sh stamps them
│           ├── MacOS/
│           │   └── launch.sh
│           └── Resources/
│               └── (icon.icns — add your own)
└── windows/
    ├── build_installer.bat    ← builds the .exe (downloads + stages Python, then compiles)
    ├── .python_cache/         ← cached Python downloads (auto-created)
    ├── python/                ← extracted interpreter, staged for Inno Setup (auto-created)
    ├── installer.iss          ← Inno Setup script
    ├── version.iss            ← generated by version.py, gitignored (never edit)
    └── launch.bat             ← Windows launcher
```

---

## Cutting a release

One tag, two runners, one commit. `.github/workflows/mafigate-release.yml` (issue #264) builds
the DMG on `macos-latest` and the installer on `windows-latest` and attaches both to a **draft**
release. Before it existed, a release needed a Mac *and* a physical Windows machine, both driven
by hand — and when the Windows box was not to hand, Mac users got the new version and Windows
users silently kept the old one.

```bash
cd streamlit_app
make version        # what this tree builds as; edit APP_VERSION first if it needs bumping
make release-tag    # the tag to push — mafigate-v<version>
git tag "$(make -s release-tag)"
git push origin "$(make -s release-tag)"
```

Both artifacts carry the build stamp described above, and a runner's checkout is clean, so About in
a released app names the tag's commit with no `-dirty` suffix. The workflow does not write the stamp
itself — the two build scripts do, which is why it does not have to (issue #263).

Then:

1. Watch the run. `gh run watch <id> --exit-status`, or the Actions tab.
2. Check the draft. Both files attached, both named from `<version>`, and the notes naming the
   commit and each file's sha256.
3. Read the release page. **You do not write it here** — its copy is
   [RELEASE_PAGE.md](RELEASE_PAGE.md) and the workflow has already assembled it into the draft,
   with the downloads block and [OPENING_MAFIGATE.md](OPENING_MAFIGATE.md) filled in. Changing
   what it says means editing that file and cutting again, which is the point: a page typed into
   the web form once per release is a page nobody reviewed. Before publishing, read
   [One thing measured but not settled](#one-thing-measured-but-not-settled-worth-doing-before-the-release-page-goes-out)
   below — it is about what this page promises macOS users, and it is not settled.
4. **`make release-preflight`.** It refuses if the draft was built from a tree the repository
   has since moved on from. Do not skip it because you only cut the tag an hour ago — that is
   exactly the state it was written for: `mafigate-v1.0.0` was drafted from the tip of `main`,
   sat while twenty-two commits landed, four of which changed code that ships, and came within
   one click of publishing a known wrong-variant bug under the fixed version's number (#328).
   The draft says nothing about this. It names its commit correctly and reads exactly like a
   fresh one.
5. Publish it. **CI never does this.** It drafts and stops, so a human decides when a download
   exists.
6. Record it in [RELEASES.md](RELEASES.md), and update the channel table in
   `streamlit_app/README.md`. That ledger is what the README's installer copy is held against,
   in both directions: until you record the release, a guard refuses a README that links to a
   download; once you do, it refuses one that still calls the installers unavailable.

### Rehearsing it

The workflow's first run must not be the one that has to produce a real artifact for clinicians.
Push the same version with a suffix — `mafigate-v<version>-rc1`, or `-rehearsal1` — and the
version check accepts it as a **pre-release** instead of a release. Everything else runs
identically. Delete the draft afterwards; it never was a download.

Rehearse in the **public** repository, where both runners are free. macOS runners bill at a
multiplier on private repositories, which is also why a tag push is this workflow's only
trigger — no ordinary commit spends a macOS minute on it.

### What the workflow refuses

- **A tag that names a version this tree does not hold.** The tag is a human typing a number;
  `build/version.py` says what it should have been, and a mismatch fails the run before either
  build starts. `mafigate-v<version>-<suffix>` is the one accepted variation.
- **Two commits.** Each build job reports the commit it checked out, and the release job
  refuses to draft anything unless both match the tag's — including refusing an *empty* report,
  since three empty strings are all equal.
- **An artifact whose name does not carry the version.** Each job reads the name off what its
  build script produced — the naming convention is not repeated in the workflow — and checks the
  derived version is in it. Two artifacts of one kind, or none, is also a refusal.
- **A rebuild of a release someone already published.** `gh release edit --draft` would retract
  it, so a re-run against a published tag stops instead. Delete the release first if a rebuild
  is really what you want.
- **An installer with no version in it.** Before building, the Windows job compiles
  `installer.iss` with no `version.iss` present and requires ISCC to stop with its `#error`.
  Inno Setup expands an undefined symbol to nothing, so the silent alternative is an
  unversioned installer; that compile is the one thing `tests/test_installer_version.py` cannot
  do for itself.
- **Drifted vendor copies.** `vendor/_sync.py --check` runs first, on both platforms' behalf.

---

## Distributing to users

The release page is the channel; the list below is for a build you made by hand.

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
| macOS blocks the app, or calls it damaged | [OPENING_MAFIGATE.md](OPENING_MAFIGATE.md) — the wording lives there, not here |
| Windows SmartScreen warning | [OPENING_MAFIGATE.md](OPENING_MAFIGATE.md) — same file, Windows section |
| Windows "bundled Python is missing" | The install is incomplete — reinstall MAFigate. Installing a Python of your own will not help; the app does not look for one. |
| Dependencies fail to install | Check internet connection; check `~/.mafigate/logs/` |
| Port 8501 in use | App automatically finds the next free port |
| App update: stale venv | Delete `~/.mafigate/venv` and relaunch |

Those first two rows used to carry their own copy of the instructions, and the two copies
disagreed — this table gave one `xattr` incantation, the DMG build script's closing message
another. The wording now lives once, in [OPENING_MAFIGATE.md](OPENING_MAFIGATE.md), and
`tests/test_unsigned_artifact_copy.py` fails if it acquires a second home, alert quotations
included. Link to it; do not paste it.

Where each audience meets that file is decided by *when* the alert fires:

- **macOS** — Gatekeeper fires when the app is opened, after the DMG window is already on
  screen. `build_dmg.sh` therefore stages the note into the mounted window as `Opening
  MAFigate.txt`, beside the app, with a position of its own in the create-dmg layout.
- **Windows** — SmartScreen fires on the downloaded installer *before* any surface of ours
  can be drawn, so nothing we ship arrives in time. The release page is the only thing
  that reaches a Windows recipient first, and it must quote this file.

### One thing measured but not settled, worth doing before the release page goes out

The note leads macOS users through **System Settings → Privacy & Security**, which is the
route Apple documents, and which is where macOS 15 moved the only override (Control-click →
Open no longer works for software that is not correctly signed *or* notarized —
[Apple Developer news](https://developer.apple.com/news/?id=saqachfa)). There is a credible
secondary claim that for an **ad-hoc-signed** app — which is what `build_dmg.sh` produces,
with `codesign --force --deep --sign -` — the alert is a dead end and no override appears in
that pane at all.

What was measured here, on macOS 26, with two probe bundles signed exactly as this script
signs (one flat, one carrying a Mach-O binary under `Contents/Resources/`, the shape the
bundled Python arrives in):

- `codesign --verify --deep --strict` reports **valid on disk** and *satisfies its Designated
  Requirement*. The ad-hoc signature is valid, so the bundle is **not** in the
  invalid-signature class that produces the "damaged" alert.
- `syspolicy_check distribution` grades *Adhoc Signed App* a **Warning** and *Notary Ticket
  Missing* the only **Fatal** — the same fatal condition a Developer-ID-signed but
  unnotarized app has, which is the class the documented override exists for.
- `spctl -a -t exec` rejects it either way, quarantined or not, as expected.

None of that proves the pane offers the button, so the note names **both** triggers for its
Terminal fallback — damaged, or no button offered — and deliberately does not say how often
either happens. **The decisive test costs one afternoon and no session here can run it:**
build the DMG, copy the app out of it, `xattr -w com.apple.quarantine` the copy, open it on
a Sequoia-or-later Mac, and record which alert fires and whether Privacy & Security offers a
line for it. Do that before the release page tells people the GUI route works — and if it
turns out it does not, the answer is notarization, not a rewritten note.
