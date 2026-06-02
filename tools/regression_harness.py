#!/usr/bin/env python3
"""
PSOPT example regression harness.

Builds and runs a representative subset of the PSOPT examples and compares
their numerical output against recorded baselines. Unlike a pure build/convergence
check, this catches *silent* numerical regressions: a change that still compiles
and still reports "Optimal Solution Found" but quietly alters the computed solution.

Per example it records up to three signals, whichever the example prints:
  * optimal_cost   - PSOPT's "Optimal (unscaled) cost function value:"
  * final_objective- IPOPT's final "Objective..........:" value (universal)
  * max_rel_error  - the largest "maximum relative local error" across phases
plus `solved` = number of "EXIT: Optimal Solution Found." lines.

Numeric signals are compared with a relative tolerance (default 1e-4) so that
last-digit, platform-dependent noise does not cause spurious failures, while a
real change in the solution (which moves digits well above 1e-4) is caught.

Usage
-----
  # 1. Configure + build the example subset (needs ADOL-C w/ ColPack, IPOPT, Eigen):
  python3 tools/regression_harness.py --build --record

  # 2. Later, after making changes, re-build and verify nothing moved:
  python3 tools/regression_harness.py --build --check

  # Run the whole catalogue instead of the default subset:
  python3 tools/regression_harness.py --build --examples all --check

Exit code is non-zero if any example fails in --check mode, so it can gate CI.
"""
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# A representative subset spanning the solver's feature space. Each entry lists
# any data files the example needs at run time (copied next to the binary).
DEFAULT_SUBSET = {
    "bryson_denham":     {"data": [], "note": "state-constrained, classic"},
    "breakwell":         {"data": [], "note": "state constraint"},
    "hypersensitive":    {"data": [], "note": "stiff, global pseudospectral"},
    "goddard":           {"data": [], "note": "singular arc, multi-phase"},
    "shuttle_reentry":   {"data": [], "note": "path constraints, aerospace"},
    "obstacle":          {"data": [], "note": "path constraints"},
    "twophase_schwartz": {"data": [], "note": "multi-phase + linkages"},
    "catmix":            {"data": [], "note": "larger NLP, chemical"},
    "predator":          {"data": ["predator.dat"], "note": "parameter estimation + data loader"},
    "param2":            {"data": ["param2.dat"], "note": "parameter estimation + data loader"},
    "dae_i3":            {"data": ["dae_i3.dat"], "note": "index-3 DAE + data loader"},
    "low_thrust":        {"data": ["T0.dat", "U0.dat", "X0.dat"], "note": "data loaders, larger"},
    "twoburn":           {"data": ["twoburn.txt"], "note": "multi-phase orbit transfer + data loader"},
}

# Examples excluded from the 'all' set. They still run if named explicitly.
EXCLUDE_FROM_ALL = {
    "climb",  # minimum-time-to-climb: does not converge within a practical CI timeout (>120s)
}

FLOAT = r"[-+]?\d+\.?\d*(?:[eE][-+]?\d+)?"
RE_COST = re.compile(r"Optimal \(unscaled\) cost function value:\s*(" + FLOAT + r")")
RE_OBJ = re.compile(r"Objective\.*:\s*(" + FLOAT + r")")
RE_ERR = re.compile(r"maximum relative local error:\s*(" + FLOAT + r")")
RE_SOLVED = re.compile(r"EXIT:\s*Optimal Solution Found")


def parse_signals(text):
    sig = {}
    costs = RE_COST.findall(text)
    objs = RE_OBJ.findall(text)
    errs = RE_ERR.findall(text)
    if costs:
        sig["optimal_cost"] = float(costs[-1])          # final mesh-refinement value
    if objs:
        sig["final_objective"] = float(objs[-1])        # IPOPT's final objective
    if errs:
        sig["max_rel_error"] = max(float(e) for e in errs)
    sig["solved"] = len(RE_SOLVED.findall(text))
    return sig


def make_gnuplot_stub(workdir):
    """A no-op 'gnuplot' so plot() calls return instantly and never block headless."""
    stub_dir = workdir / "_stubbin"
    stub_dir.mkdir(exist_ok=True)
    stub = stub_dir / "gnuplot"
    stub.write_text("#!/bin/sh\nexit 0\n")
    stub.chmod(0o755)
    return str(stub_dir)


def configure_and_build(build_dir, source_dir, targets, jobs):
    build_dir.mkdir(parents=True, exist_ok=True)
    if not (build_dir / "CMakeCache.txt").exists():
        subprocess.run(
            ["cmake", str(source_dir), "-DCMAKE_BUILD_TYPE=Release",
             "-DBUILD_EXAMPLES=ON", "-DHEADLESS=ON"],
            cwd=build_dir, check=True,
        )
    subprocess.run(["make", "PSOPT", "-j", str(jobs)], cwd=build_dir, check=True)
    subprocess.run(["make", *targets, "-j", str(jobs)], cwd=build_dir, check=True)


def find_binary(build_dir, name):
    for p in build_dir.rglob(name):
        if p.is_file() and os.access(p, os.X_OK):
            return p
    return None


def run_example(name, spec, build_dir, source_dir, timeout, stub_path):
    binary = find_binary(build_dir, name)
    if binary is None:
        return None, "binary not found (build it first)"
    rundir = binary.parent
    for d in spec.get("data", []):
        src = source_dir / "examples" / name / d
        if src.exists():
            shutil.copy(src, rundir / d)
    env = dict(os.environ, PATH=stub_path + os.pathsep + os.environ.get("PATH", ""))
    try:
        proc = subprocess.run([f"./{name}"], cwd=rundir, env=env, timeout=timeout,
                              capture_output=True, text=True)
    except subprocess.TimeoutExpired:
        return None, f"timeout after {timeout}s"
    sig = parse_signals(proc.stdout + proc.stderr)
    if sig.get("solved", 0) < 1:
        return sig, "did not report Optimal Solution Found"
    return sig, None


def compare(name, baseline, current, rtol):
    msgs = []
    ok = True
    if current.get("solved", 0) < baseline.get("solved", 0):
        ok = False
        msgs.append(f"solved {current.get('solved',0)} < baseline {baseline['solved']}")
    for key in ("optimal_cost", "final_objective", "max_rel_error"):
        if key in baseline:
            if key not in current:
                ok = False
                msgs.append(f"{key} missing")
                continue
            b, c = baseline[key], current[key]
            denom = max(abs(b), 1e-30)
            rel = abs(c - b) / denom
            # max_rel_error is a diagnostic that legitimately varies more; only
            # flag order-of-magnitude moves on it.
            tol = rtol if key != "max_rel_error" else max(rtol, 0.5)
            if rel > tol:
                ok = False
                msgs.append(f"{key}: {c:.8g} vs {b:.8g} (rel {rel:.1e} > {tol:.0e})")
    return ok, "; ".join(msgs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-dir", default=str(Path(__file__).resolve().parent.parent))
    ap.add_argument("--build-dir", default=None,
                    help="cmake build directory (default: <source>/build_regression)")
    ap.add_argument("--baselines", default=None,
                    help="baseline json (default: <source>/tools/regression_baselines.json)")
    ap.add_argument("--examples", default="subset",
                    help="'subset' (default), 'all', or comma-separated names")
    ap.add_argument("--build", action="store_true", help="configure + build before running")
    ap.add_argument("--record", action="store_true", help="write baselines from this run")
    ap.add_argument("--check", action="store_true", help="compare this run to baselines")
    ap.add_argument("--rtol", type=float, default=1e-4)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    args = ap.parse_args()

    source_dir = Path(args.source_dir).resolve()
    build_dir = Path(args.build_dir).resolve() if args.build_dir else source_dir / "build_regression"
    baselines_path = Path(args.baselines) if args.baselines else source_dir / "tools" / "regression_baselines.json"

    if args.examples == "subset":
        specs = DEFAULT_SUBSET
    elif args.examples == "all":
        ex_root = source_dir / "examples"
        specs = {p.name: {"data": [f.name for f in p.glob("*.dat")] + [f.name for f in p.glob("*.txt")]}
                 for p in sorted(ex_root.iterdir())
                 if p.is_dir() and p.name not in EXCLUDE_FROM_ALL}
    else:
        ex_root = source_dir / "examples"
        names = [n.strip() for n in args.examples.split(",") if n.strip()]
        specs = {}
        for n in names:
            if n in DEFAULT_SUBSET:
                specs[n] = DEFAULT_SUBSET[n]
            else:
                d = ex_root / n
                specs[n] = {"data": [f.name for f in d.glob("*.dat")] + [f.name for f in d.glob("*.txt")]}

    if args.build:
        configure_and_build(build_dir, source_dir, list(specs.keys()), args.jobs)

    stub_path = make_gnuplot_stub(build_dir)

    results = {}
    errors = {}
    print(f"{'example':<20} {'cost/objective':>20} {'max_rel_err':>14}  status")
    print("-" * 72)
    for name, spec in specs.items():
        sig, err = run_example(name, spec, build_dir, source_dir, args.timeout, stub_path)
        if sig is None:
            errors[name] = err
            print(f"{name:<20} {'-':>20} {'-':>14}  ERROR: {err}")
            continue
        results[name] = sig
        primary = sig.get("optimal_cost", sig.get("final_objective"))
        perr = sig.get("max_rel_error")
        note = "" if err is None else f"  ({err})"
        print(f"{name:<20} {primary:>20.8g} "
              f"{(perr if perr is not None else float('nan')):>14.3g}  ok{note}")

    if args.record:
        baselines_path.parent.mkdir(parents=True, exist_ok=True)
        baselines_path.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
        print(f"\nWrote {len(results)} baselines to {baselines_path}")
        return 0

    if args.check:
        if not baselines_path.exists():
            print(f"\nNo baselines at {baselines_path}; run with --record first.", file=sys.stderr)
            return 2
        baselines = json.loads(baselines_path.read_text())
        print("\n" + "=" * 72)
        n_fail = 0
        for name in specs:
            if name in errors:
                n_fail += 1
                print(f"FAIL {name}: {errors[name]}")
                continue
            if name not in baselines:
                print(f"NEW  {name}: no baseline (run --record to add)")
                continue
            ok, msg = compare(name, baselines[name], results[name], args.rtol)
            if ok:
                print(f"PASS {name}")
            else:
                n_fail += 1
                print(f"FAIL {name}: {msg}")
        print("=" * 72)
        print(f"{len(specs) - n_fail}/{len(specs)} passed")
        return 1 if n_fail else 0

    return 0


if __name__ == "__main__":
    sys.exit(main())
