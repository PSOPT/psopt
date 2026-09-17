"""Shared bootstrap for the examples in this directory.

Every example begins with ``from _common import psopt``. That import puts the
source tree's ``python/`` directory on sys.path, so an example runs straight from
a checkout with nothing installed:

    cd python/examples && python3 brachistochrone.py

If PSOPT's Python package is installed instead, the installed one is used. The
examples used to open with a hard-coded ``sys.path.insert(0, "/tmp/psopt_py")``,
which worked on exactly one machine.
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PKG_PARENT = os.path.dirname(_HERE)          # .../psopt/python

try:
    import psopt                               # installed
except ImportError:
    if _PKG_PARENT not in sys.path:
        sys.path.insert(0, _PKG_PARENT)
    import psopt                               # from the source tree


def report(sol, label, reference=None, tol=None):
    """Print the result and, when a reference is given, check it.

    Returns True when the solve succeeded and matched. An example that prints a
    number without saying whether the solve converged is not much of a check,
    which is what sol.status is for.
    """
    ok = sol.status.success
    print("\n---- %s" % label)
    print("objective      : %.12g" % sol.objective)
    if reference is not None:
        rel = abs(sol.objective - reference) / max(abs(reference), 1.0)
        print("reference      : %.12g   (relative difference %.2e)" % (reference, rel))
        if tol is not None:
            ok = ok and rel <= tol
    print("status         : %s" % sol.status)
    if not sol.status.success:
        print("*** the solve did NOT converge:", sol.status.error_msg or
              "NLP return code %d" % sol.status.nlp_return_code)
    return ok
