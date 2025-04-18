#!/usr/bin/env python3
"""
Run selected PSOPT examples, parse each *psopt_solution_<name>.txt* file,
and decide pass / fail based on

  • IPOPT success string          → “NLP solver reports: … The problem has been solved!”
  • Optimal cost function value   → line starts with
                                    “Optimal (unscaled) cost function value:”

The list of examples and their reference costs are provided as two
comma‑separated environment variables or CLI flags:

  EXAMPLES="twoburn,low_thrust,zpm,shuttle_reentry,launch"
  REF_COSTS="-2.367249e-01,-2.203380e-01,6.680110e+06,-3.414119e+01,-7.529661e+03"

Tolerance is relative  (`abs(cost-ref)/abs(cost) <= tol`) where the default
is 0.01 
"""

import argparse, os, pathlib, subprocess, json, sys, re, math, textwrap
from datetime import datetime

def parse_args():
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent(__doc__))
    p.add_argument("--exe-dir", required=True,
                   help="folder that holds example sub‑dirs with executables")
    p.add_argument("--summary", required=True,
                   help="write a JSON summary here")
    p.add_argument("--examples", default=os.getenv("EXAMPLES",""))
    p.add_argument("--ref-costs", default=os.getenv("REF_COSTS",""))
    p.add_argument("--tol", type=float, default=None,
                   help="global abs tolerance (optional)")
    return p.parse_args()

def default_tol(ref: float, override: float | None):
    """
    Return the tolerance to use for a given reference value.

    • If the user passed --tol, honour that.
    • Otherwise use 0.01  (1 % relative error).
    """
    if override is not None:
        return override
    return 0.01

success_re   = re.compile(r"NLP solver reports:\s*The problem has been solved!", re.I)
cost_line_re = re.compile(r"Optimal \(unscaled\) cost function value:\s+([-+0-9.eE]+)")

def run_example(exe_path: pathlib.Path, name: str):
    log = {"passed": False, "cost": None, "elapsed_s": None, "solver_ok": False}
    start = datetime.now()
    try:
        subprocess.run([str(exe_path)], check=True, timeout=600, cwd=exe_path.parent)
 # 10 min timeout
    except subprocess.CalledProcessError as e:
        log["error"] = f"exited with code {e.returncode}"
        return log
    except subprocess.TimeoutExpired:
        log["error"] = "timeout"
        return log

# ------------------------------------------------------------------
# locate the solution file (pattern is always psopt_solution_*.txt)
# ------------------------------------------------------------------
    candidates = list(exe_path.parent.glob("psopt_solution_*.txt"))
    if not candidates:                       # nothing next to the executable?
      candidates = list(pathlib.Path.cwd().glob("psopt_solution_*.txt"))

    if not candidates:
      log["error"] = "solution file not found"
      return log

# if several exist, take the most recently modified one
    sol_file = max(candidates, key=lambda p: p.stat().st_mtime)
    log["solution_file"] = str(sol_file)     # optional: keep for debugging

    with sol_file.open() as fh:
        for line in fh:
            if cost_line_re.search(line):
                log["cost"] = float(cost_line_re.search(line).group(1))
            if success_re.search(line):
                log["solver_ok"] = True
    log["elapsed_s"] = (datetime.now()-start).total_seconds()
    return log

def main():
    args   = parse_args()
    ex_dir = pathlib.Path(args.exe_dir).resolve()

    examples = [e.strip() for e in args.examples.split(",") if e.strip()]
    ref_vals = [float(x) for x in args.ref_costs.split(",")] if args.ref_costs else []
    if ref_vals and len(ref_vals)!=len(examples):
        print("::error ::Number of REF_COSTS entries does not match EXAMPLES list")
        sys.exit(1)
    ref_map  = dict(zip(examples, ref_vals))

    summary, all_ok = {}, True
    for ex in examples:
        exe = ex_dir / ex / ex          # ./build/examples/twoburn/twoburn
        result = run_example(exe, ex)
        ref    = ref_map.get(ex)
        tol    = default_tol(ref, args.tol) if ref is not None else None

# ───────── cost‑closeness check (relative unless ref≈0) ──────────
    if ref is not None and result["cost"] is not None:
      eps = 1e-12                     # threshold for “too small to divide by”
      if abs(ref) < eps:              # reference is zero → use absolute diff
          rel_err   = None
          abs_err   = abs(result["cost"] - ref)
          result["cost_pass"] = abs_err <= tol       # tol interpreted ABSOLUTELY
      else:                           # normal relative‑error check
          rel_err   = abs(result["cost"] - ref) / abs(ref)
          result["cost_pass"] = rel_err <= tol       # tol interpreted RELATIVELY
    result["rel_error"] = rel_err                  # may be None if abs test
    else:
        result["cost_pass"] = (ref is None)            # no reference ⇒ skip
# ─────────────────────────────────────────────────────────────────


        # overall pass flag
        result["passed"] = result["solver_ok"] and result["cost_pass"]
        if not result["passed"]:
            all_ok = False
            print(f"::error ::{ex} failed  "
                  f"(solver_ok={result['solver_ok']}, cost={result['cost']}, ref={ref})")

        summary[ex] = result

    # write JSON summary
    pathlib.Path(args.summary).write_text(json.dumps(summary, indent=2))
    sys.exit(0 if all_ok else 1)

if __name__ == "__main__":
    main()
