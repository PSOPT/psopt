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

Tolerance is absolute  (`abs(cost-ref) <= tol`) where the default
is 1e‑6 for |ref| < 1 and 1e‑3 otherwise; override with `--tol` CLI.
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

def default_tol(ref: float, override: float|None):
    if override is not None:
        return override
    return 1e-6 if abs(ref) < 1.0 else 1e-3

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

    # locate solution file
    sol_file = exe_path.parent / f"psopt_solution_{name}.txt"
    if not sol_file.exists():
        log["error"] = "solution file missing"
        return log

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

        # numerical check
        if ref is not None and result["cost"] is not None:
            result["cost_pass"] = abs(result["cost"]-ref) <= tol
        else:
            result["cost_pass"] = (ref is None)  # if no ref provided, ignore

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
