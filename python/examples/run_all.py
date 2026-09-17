#!/usr/bin/env python3
"""Run every example in this directory and report which ones pass.

    python3 run_all.py            # all of them
    python3 run_all.py obstacle   # just the ones whose name contains "obstacle"

Each example prints its own result and exits non-zero if it did not match its
reference, so this is a thin wrapper: it runs them in separate processes, keeps
the last few lines of each, and prints a table. Nothing here knows anything about
any particular example.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP = {"_common.py", "run_all.py", "rv2oe_casadi.py"}   # helpers, not examples


def main(argv):
    names = sorted(f for f in os.listdir(HERE)
                   if f.endswith(".py") and f not in SKIP
                   and (not argv or any(a in f for a in argv)))
    if not names:
        print("no examples matched")
        return 1

    width = max(len(n) for n in names)
    failures = []
    print("%-*s %8s %10s   %s" % (width, "example", "result", "seconds", "objective"))
    for n in names:
        t0 = time.time()
        r = subprocess.run([sys.executable, n], cwd=HERE,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        dt = time.time() - t0
        out = r.stdout.decode("utf-8", "replace").splitlines()
        obj = next((l.split(":", 1)[1].strip() for l in reversed(out)
                    if l.startswith("objective") and ":" in l), "")
        print("%-*s %8s %10.1f   %s"
              % (width, n, "pass" if r.returncode == 0 else "FAIL", dt, obj))
        if r.returncode != 0:
            failures.append((n, out[-12:]))

    if failures:
        for n, tail in failures:
            print("\n---- %s, last lines:" % n)
            for l in tail:
                print("   ", l)
    print("\n%d of %d passed" % (len(names) - len(failures), len(names)))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
