Opening MAFigate the first time
===============================

MAFigate is not signed with a paid Apple or Microsoft developer certificate. So the first
time you open it, your computer warns you that it cannot check who made it. The warning is
expected, and it is not a report about MAFigate: we have not bought a certificate, so
there is nothing for your computer to check the app against.

Find your system below and follow the steps. They are done by clicking; only the last
section needs a Terminal, and only if the steps for your system did not work.

macOS
-----

Double-click MAFigate in your Applications folder. An alert appears, and it says:

    Apple could not verify "MAFigate" is free of malware that may harm your Mac or
    compromise your privacy.

Its only buttons are "Move to Trash" and "Done". There is no "Open" button, and that is
normal — the permission is given in Settings instead.

1. Click "Done". Do not click "Move to Trash".
2. Open the Apple menu, then System Settings, then Privacy & Security.
3. Scroll down to Security. It now shows a line saying MAFigate was blocked, with an
   "Open Anyway" button beside it. Click that button.
4. A box asks you to confirm, and its wording is worth recognising:
   "Opening MAFigate will always allow it to run on this Mac."
   Click "Open", then enter your login password.
5. MAFigate opens. Every launch after this one opens normally, with no alert.

If your Mac runs macOS 14 or older, there is a shorter route: Control-click (right-click)
MAFigate in your Applications folder, choose "Open", then click "Open" in the alert that
appears. macOS 15 removed that shortcut, so on a current Mac use the five steps above.

Windows
-------

Windows blocks the installer before it starts, on a blue full-screen page.
It is headed "Windows protected your PC", and the only button you can see is "Don't run".

1. Click "More info", the small link under the message.
2. Click the "Run anyway" button that appears.
3. The installer continues normally.

You see this once, for the installer itself. Launching MAFigate afterwards shows nothing.

If the steps above did not work — the one step that needs Terminal
------------------------------------------------------------------

Two things end up here: macOS calls MAFigate damaged rather than unverified, or Privacy &
Security never offers you a button to click. Either way it takes one command, and this is
the only instruction on this page that needs the Terminal.

Open Terminal (Applications, then Utilities), paste the line below, and press Return:

    xattr -dr com.apple.quarantine /Applications/MAFigate.app

Then open MAFigate again.

Why the warnings happen at all
------------------------------

Removing them means buying a Developer ID certificate from Apple and a code-signing
certificate for Windows, and sending every build to Apple to be notarized. Until that
happens, both downloads stay unsigned and both warnings stay. Neither warning means a scan
found something — each one says only that the signature is missing.
