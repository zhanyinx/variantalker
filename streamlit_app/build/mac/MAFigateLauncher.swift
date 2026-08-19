// MAFigate native launcher
// =============================================================================
// A minimal AppKit application used as the .app's main executable. Its only job
// is to run a proper macOS event loop (so the OS does not flag the app as "not
// responding") while a child process (launch.sh) sets up the environment and
// runs the Streamlit server. When the user quits the app (Cmd-Q, Dock -> Quit,
// or the server exiting), the child is terminated cleanly so the Streamlit
// server is not left orphaned.
//
// Previously launch.sh was the bundle executable directly; because it blocks on
// `wait` and never services the AppKit run loop, macOS marked the app as not
// responding and it had to be force-quit.

import AppKit
import Foundation

final class AppDelegate: NSObject, NSApplicationDelegate {
    private var child: Process?
    private var pendingFile: String?

    // Document open (e.g. double-clicking a .maf file). For launch-time opens
    // this is delivered before applicationDidFinishLaunching, so we stash the
    // path and hand it to launch.sh when we start it.
    func application(_ sender: NSApplication, openFile filename: String) -> Bool {
        if child == nil {
            pendingFile = filename
        }
        return true
    }

    func applicationDidFinishLaunching(_ notification: Notification) {
        let scriptPath = Bundle.main.bundlePath + "/Contents/MacOS/launch.sh"

        let proc = Process()
        proc.executableURL = URL(fileURLWithPath: "/bin/bash")
        var args = [scriptPath]
        if let file = pendingFile {
            args.append(file)
        }
        proc.arguments = args

        // If the server process exits on its own, quit the wrapper too.
        proc.terminationHandler = { _ in
            DispatchQueue.main.async { NSApp.terminate(nil) }
        }

        do {
            try proc.run()
            child = proc
        } catch {
            NSApp.terminate(nil)
        }
    }

    func applicationShouldTerminate(_ sender: NSApplication) -> NSApplication.TerminateReply {
        // SIGTERM the launcher script; its trap stops the Streamlit server.
        if let c = child, c.isRunning {
            c.terminationHandler = nil
            c.terminate()
            c.waitUntilExit()
        }
        return .terminateNow
    }
}

// Build a minimal main menu so Cmd-Q (and Dock -> Quit) work.
func makeMainMenu() -> NSMenu {
    let mainMenu = NSMenu()
    let appItem = NSMenuItem()
    mainMenu.addItem(appItem)

    let appMenu = NSMenu()
    let name = ProcessInfo.processInfo.processName
    appMenu.addItem(
        withTitle: "Quit \(name)",
        action: #selector(NSApplication.terminate(_:)),
        keyEquivalent: "q"
    )
    appItem.submenu = appMenu
    return mainMenu
}

let app = NSApplication.shared
let delegate = AppDelegate()
app.delegate = delegate
app.setActivationPolicy(.regular)
app.mainMenu = makeMainMenu()
app.run()
