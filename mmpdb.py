from mmpdblib import commandline
import signal
import os

if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal.SIG_DFL)  # Allow ^C to kill the process.
    if not os.name == "nt":
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)  # Allow the output pipe to be closed
    commandline.main()
