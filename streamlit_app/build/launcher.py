#!/usr/bin/env python3
"""
MAFigate Desktop Launcher

Cross-platform launcher that:
1. Finds or creates a Python virtual environment
2. Installs dependencies if needed
3. Starts the Streamlit server on localhost
4. Opens the browser automatically
5. Shuts down cleanly on exit
"""

import os
import sys
import subprocess
import webbrowser
import time
import signal
import socket
import threading
from pathlib import Path


APP_NAME = "MAFigate"
DEFAULT_PORT = 8501
MAX_PORT_ATTEMPTS = 10
STARTUP_TIMEOUT = 30  # seconds


def get_app_dir():
    """Get the directory containing the Streamlit app."""
    # When bundled, the app is next to the launcher or in a known relative path
    launcher_dir = Path(__file__).resolve().parent
    # Check if we're in build/ — app is one level up
    app_dir = launcher_dir.parent
    if (app_dir / "MAFigate.py").exists():
        return app_dir
    # Check if app is bundled alongside
    if (launcher_dir / "MAFigate.py").exists():
        return launcher_dir
    # Check common bundle layouts
    for candidate in [launcher_dir / "app", launcher_dir / "streamlit_app"]:
        if (candidate / "MAFigate.py").exists():
            return candidate
    print(f"ERROR: Could not find MAFigate.py. Looked in: {launcher_dir}, {app_dir}")
    sys.exit(1)


def find_free_port(start_port=DEFAULT_PORT):
    """Find an available port starting from start_port."""
    for port in range(start_port, start_port + MAX_PORT_ATTEMPTS):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(("127.0.0.1", port))
                return port
        except OSError:
            continue
    print(f"ERROR: No free port found in range {start_port}-{start_port + MAX_PORT_ATTEMPTS}")
    sys.exit(1)


def get_venv_dir(app_dir):
    """Get the virtual environment directory."""
    return app_dir / ".mafigate_venv"


def ensure_venv(app_dir):
    """Create virtual environment and install dependencies if needed."""
    venv_dir = get_venv_dir(app_dir)
    requirements = app_dir / "requirements.txt"

    if sys.platform == "win32":
        python_bin = venv_dir / "Scripts" / "python.exe"
        pip_bin = venv_dir / "Scripts" / "pip.exe"
        streamlit_bin = venv_dir / "Scripts" / "streamlit.exe"
    else:
        python_bin = venv_dir / "bin" / "python"
        pip_bin = venv_dir / "bin" / "pip"
        streamlit_bin = venv_dir / "bin" / "streamlit"

    # Create venv if it doesn't exist
    if not python_bin.exists():
        print(f"Setting up {APP_NAME} environment (first run, may take a minute)...")
        subprocess.run(
            [sys.executable, "-m", "venv", str(venv_dir)],
            check=True,
        )

    # Install/update dependencies
    marker = venv_dir / ".deps_installed"
    req_mtime = requirements.stat().st_mtime if requirements.exists() else 0
    marker_mtime = marker.stat().st_mtime if marker.exists() else 0

    if req_mtime > marker_mtime:
        print("Installing dependencies...")
        subprocess.run(
            [str(pip_bin), "install", "-q", "-r", str(requirements)],
            check=True,
        )
        marker.touch()

    return streamlit_bin


def wait_for_server(port, timeout=STARTUP_TIMEOUT):
    """Wait until the Streamlit server is accepting connections."""
    start = time.time()
    while time.time() - start < timeout:
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.settimeout(1)
                s.connect(("127.0.0.1", port))
                return True
        except (ConnectionRefusedError, OSError):
            time.sleep(0.5)
    return False


def main():
    app_dir = get_app_dir()
    port = find_free_port()

    print(f"Starting {APP_NAME}...")

    # Ensure virtual environment and dependencies
    streamlit_bin = ensure_venv(app_dir)

    # Start Streamlit server
    env = os.environ.copy()
    env["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"
    env["STREAMLIT_SERVER_HEADLESS"] = "true"

    process = subprocess.Popen(
        [
            str(streamlit_bin),
            "run",
            str(app_dir / "MAFigate.py"),
            f"--server.port={port}",
            "--server.address=127.0.0.1",
            "--server.headless=true",
        ],
        env=env,
        cwd=str(app_dir),
    )

    url = f"http://127.0.0.1:{port}"

    # Handle clean shutdown
    def shutdown(signum=None, frame=None):
        print(f"\nShutting down {APP_NAME}...")
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
        sys.exit(0)

    signal.signal(signal.SIGINT, shutdown)
    signal.signal(signal.SIGTERM, shutdown)
    if sys.platform != "win32":
        signal.signal(signal.SIGHUP, shutdown)

    # Wait for server and open browser
    print(f"Waiting for server on {url} ...")
    if wait_for_server(port):
        print(f"{APP_NAME} is running at {url}")
        webbrowser.open(url)
    else:
        print(f"WARNING: Server may not have started properly. Try opening {url} manually.")

    print(f"Press Ctrl+C to stop {APP_NAME}")

    # Keep running until the process exits or user interrupts
    try:
        process.wait()
    except KeyboardInterrupt:
        shutdown()


if __name__ == "__main__":
    main()
