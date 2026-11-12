#!/usr/bin/env python3
"""Detached launcher: run a command in its own session; return at once.

Usage: python3 scripts/detach.py CMD [ARGS...]    (env DETACH_LOG=path.log optional)

Learned pattern (2026-11-11/12): never run long compute as lane foreground --
fork + setsid + dup2(log) + execvp; the child gets PYTHONUNBUFFERED=1 so log
tails stay live, and it survives terminal SIGINTs (the stray-SIGINT lesson).
"""
import os
import sys


def main():
    cmd = sys.argv[1:]
    if not cmd:
        print("usage: detach.py CMD [ARGS...]")
        return 2
    log = os.environ.get("DETACH_LOG")
    pid = os.fork()
    if pid == 0:
        os.setsid()
        fd = (os.open(log, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
              if log else os.open(os.devnull, os.O_WRONLY))
        os.dup2(fd, 1)
        os.dup2(fd, 2)
        nul = os.open(os.devnull, os.O_RDONLY)
        os.dup2(nul, 0)
        os.environ["PYTHONUNBUFFERED"] = "1"
        os.execvp(cmd[0], cmd)
        os._exit(127)  # unreachable
    sys.stdout.write("DETACHED pid=%d cmd=%s log=%s\n" % (pid, " ".join(cmd), log or "-"))
    sys.stdout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
