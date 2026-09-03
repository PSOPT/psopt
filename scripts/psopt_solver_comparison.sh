#!/usr/bin/env bash
#
# psopt_solver_comparison.sh -- run every shipped PSOPT example under IPOPT and under
#                               PSOPT's own SQP, and record what each does.
#
# The comparison this produces is the evidence for the claim that IPOPT solves more of
# the example set than the SQP does and solves it faster. Earlier figures were measured
# on a two-core cloud sandbox with a 120-second limit, which is not a defensible basis
# for a statement in a book: on that machine the limit was doing much of the deciding.
# This script exists to redo it on real hardware with a limit generous enough that a
# timeout means something.
#
#   ./psopt_solver_comparison.sh --build ~/GitHub/psopt/build
#   ./psopt_solver_comparison.sh --build ~/GitHub/psopt/build --limit 3600
#   ./psopt_solver_comparison.sh --help
#
# Output: results.csv (one row per example per solver) and results.md (a summary), in
# the output directory. Send both back; the CSV is what the numbers are drawn from.
#
# Requirements
# ------------
# A PSOPT build with the SQP and a QP backend, configured with the environment-override
# option so that each example can be run under either solver without editing its source:
#
#   cmake -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON \
#         -DWITH_SQP=ON -DWITH_GALAHAD=ON -DGALAHAD_DIR=<prefix> \
#         -DPSOPT_ALLOW_ENV_OVERRIDES=ON
#   cmake --build build -j
#
# PSOPT_ALLOW_ENV_OVERRIDES is off by default and must be on, or every run silently uses
# whichever solver the example's own source names and the comparison is meaningless.
# The script refuses to start if it cannot confirm the option took effect.
#
set -euo pipefail
 
BUILD=""
OUTDIR="psopt_comparison_$(date +%Y%m%d_%H%M%S)"
LIMIT=1800          # seconds per example per solver
ONLY=""
SKIP=""
QP="GALAHAD"        # which QP backend the SQP runs use; it matters more than expected
TRUST_REGION=""     # "box" or "l2"; empty leaves the build's default alone
SOLVERS="IPOPT SQP" # which of the two to run; SQP alone halves a sweep that is
                    # comparing two SQP configurations against each other
FIXED_MESH=0        # compare on each example's initial mesh, with refinement off
 
say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
info() { printf '    %s\n' "$*"; }
die()  { printf '\n\033[1;31mError: %s\033[0m\n\n' "$*" >&2; exit 1; }
warn() { printf '\n\033[1;33mWarning: %s\033[0m\n\n' "$*" >&2; }
 
usage() {
    awk 'NR>1 { if ($0 !~ /^#/) exit; sub(/^# ?/, ""); print }' "$0"
    cat <<'EOF'
 
Options:
  --build DIR    the PSOPT build directory (required)
  --out DIR      where to write results       (default: psopt_comparison_<timestamp>)
  --limit SECS   time limit per run           (default: 1800)
  --only LIST    comma-separated example names, for a quick trial run
  --skip LIST    comma-separated example names to leave out (e.g. unpublished work)
  --qp NAME      QP backend for the SQP runs (default GALAHAD; try PIQP, OSQP, Clarabel)
  --sqp-only     run only the SQP, not IPOPT. When the question is which of two SQP
                 configurations is better -- two backends, or a box against a Euclidean
                 trust region -- the IPOPT column is identical in both runs and costs
                 half the sweep to produce twice.
  --trust-region SHAPE
                 shape of the SQP's trust region: "box" (the default) or "l2". The
                 Euclidean region is a second-order cone, so it needs a conic backend --
                 --qp Clarabel -- and every other backend refuses it outright.
  --fixed-mesh   turn mesh refinement off, so that both solvers answer the same
                 discretized problem. Without this the objective column is not a
                 controlled comparison: each solver's answer steers its own next mesh.
  --help         this message
 
Runs are sequential, deliberately: wall-clock times are comparable between examples
only if they are not competing for cores. Please do not run anything else demanding
on the machine while it works.
EOF
    exit 0
}
 
while [ $# -gt 0 ]; do
    case "$1" in
        --build) BUILD="$2"; shift 2 ;;
        --out)   OUTDIR="$2"; shift 2 ;;
        --limit) LIMIT="$2"; shift 2 ;;
        --only)  ONLY="$2"; shift 2 ;;
        --skip)  SKIP="$2"; shift 2 ;;
        --qp)    QP="$2";   shift 2 ;;
        --trust-region) TRUST_REGION="$2"; shift 2 ;;
        --sqp-only) SOLVERS="SQP"; shift ;;
        --fixed-mesh) FIXED_MESH=1; shift ;;
        --help|-h) usage ;;
        *) die "unknown option '$1' (try --help)" ;;
    esac
done
 
[ -n "$BUILD" ] || die "give --build, the PSOPT build directory (try --help)"
[ -d "$BUILD/examples" ] || die "no examples directory under $BUILD -- was it configured with -DBUILD_EXAMPLES=ON?"
 
# Each example is run with the working directory set to its own build directory, because
# several of them open data files by a path relative to it. Anything this script writes
# while that is true -- the per-run logs above all -- must therefore be named absolutely,
# or the redirect lands nowhere and every example is recorded as having failed in zero
# seconds. Resolve both paths here, once, rather than at each use.
abspath() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *)  printf '%s\n' "$PWD/${1#./}" ;;
    esac
}
BUILD=$(abspath "$BUILD")
OUTDIR=$(abspath "$OUTDIR")
 
# ---------------------------------------------------------------------------------
# Preconditions, checked rather than assumed.
#
# Both of these have already cost a wasted sweep. Without the environment overrides
# every example runs its own solver and the SQP column silently becomes a second IPOPT
# column; without OMP_CANCELLATION GALAHAD cannot solve anything and the SQP appears to
# fail everywhere. Each is invisible in the output unless looked for.
# ---------------------------------------------------------------------------------
say "Checking the build"
 
command -v gtimeout >/dev/null 2>&1 && TIMEOUT=gtimeout || TIMEOUT=timeout
command -v "$TIMEOUT" >/dev/null 2>&1 || die "no timeout command. On macOS: brew install coreutils (for gtimeout), or MacPorts: sudo port install coreutils"
info "timeout command: $TIMEOUT"
 
# Every probe below depends on the environment overrides being compiled in: without them
# each example silently uses whatever its own source names, and the SQP column becomes a
# second IPOPT column. The build records whether the option was given, so ask the cache
# rather than inferring it from behaviour further down. Inferring it produces a confident
# diagnosis of the wrong thing -- a build reconfigured without the option failed the
# fixed-mesh check below, which then blamed the age of the tree.
if [ -f "$BUILD/CMakeCache.txt" ]; then
    case "$(grep -m1 '^PSOPT_ALLOW_ENV_OVERRIDES:' "$BUILD/CMakeCache.txt" 2>/dev/null)" in
        *=ON) info "environment overrides: enabled in this build" ;;
        "")   warn "this build's CMakeCache.txt does not mention PSOPT_ALLOW_ENV_OVERRIDES,
    which is expected only for a tree older than the option. The probes below
    will still catch it if the overrides turn out not to work." ;;
        *)    die "this build was configured with PSOPT_ALLOW_ENV_OVERRIDES off, so none of
    the environment variables this script relies on will do anything: every example
    would run whichever solver, backend and mesh mode its own source names, and the
    comparison would be meaningless. Reconfigure and rebuild:
      cmake -B $BUILD -DPSOPT_ALLOW_ENV_OVERRIDES=ON <your other options>
      cmake --build $BUILD -j
    A reconfigure that omits the option turns it back off, so it has to be passed
    every time, not only the first." ;;
    esac
else
    warn "no CMakeCache.txt under $BUILD, so PSOPT_ALLOW_ENV_OVERRIDES cannot be confirmed."
fi

PROBE=""
for c in brac1 hypersensitive interior_point; do
    [ -x "$BUILD/examples/$c/$c" ] && PROBE="$c" && break
done
[ -n "$PROBE" ] || die "could not find a small example to probe with under $BUILD/examples"
 
# --fixed-mesh is carried by PSOPT_MESH_REFINEMENT, and that override arrived later than
# the others, with the fixed-mesh mode itself. A build predating it, but built with
# -DPSOPT_ALLOW_ENV_OVERRIDES=ON, passes every other check here -- PSOPT_NLP_METHOD and
# PSOPT_QP_SOLVER work -- and then ignores PSOPT_MESH_REFINEMENT in silence. The sweep
# completes, looks entirely normal, and every row was produced under automatic refinement.
# That has happened: of four sweeps run with --fixed-mesh, two came back automatic, and
# nothing in the output said so. It was only caught afterwards, by noticing that IPOPT's
# own mesh counts differed between files in which they should have been identical.
#
# The test is two-way on purpose. Looking for the override's announcement is not enough:
# PSOPT announces an override only when it changes something, and several examples leave
# mesh_refinement at the library default, so a probe that is already manual announces
# nothing whether the override works or not. Forcing each setting in turn and reading back
# what PSOPT says it used answers the question directly, and needs no assumption about
# what the probe's source asks for.
if [ "$FIXED_MESH" = "1" ]; then
    # Runs the probe forcing one mesh setting; echoes the setting PSOPT reports using.
    # PSOPT writes that into the example's own output file, which is regenerated by every
    # run, so stale copies are cleared first.
    mesh_probe_setting() {
        ( cd "$BUILD/examples/$PROBE" && rm -f ./*.txt 2>/dev/null || true )
        ( cd "$BUILD/examples/$PROBE" && env OMP_CANCELLATION=TRUE OMP_PROC_BIND=TRUE \
            PSOPT_NLP_METHOD=IPOPT PSOPT_MESH_REFINEMENT="$1" \
            "$TIMEOUT" 300 "./$PROBE" > /dev/null 2>&1 ) || true
        grep -h "^MESH REFINEMENT:" "$BUILD/examples/$PROBE"/*.txt 2>/dev/null \
            | tail -1 | awk '{print $3}'
    }
    got_manual=$(mesh_probe_setting manual)
    got_auto=$(mesh_probe_setting automatic)

    if [ -z "$got_manual" ] || [ -z "$got_auto" ]; then
        warn "could not confirm that --fixed-mesh takes effect: the probe '$PROBE' wrote no
    MESH REFINEMENT line. Check one sweep log by hand before trusting the objective
    column -- if the run was really automatic, that column compares two pipelines
    rather than two solvers."
    elif [ "$got_manual" = "$got_auto" ]; then
        die "--fixed-mesh was asked for, but PSOPT_MESH_REFINEMENT does not change anything:
    the probe '$PROBE' reported '$got_manual' whichever value it was given.
    The overrides are enabled, so this build is older than the fixed-mesh mode itself.
    Rebuild from a tree whose src/NLP_interface.cxx has
    psopt_apply_mesh_environment_override, keeping -DPSOPT_ALLOW_ENV_OVERRIDES=ON.
    Going on would produce a sweep that looks entirely normal and silently answers a
    different question."
    else
        info "fixed mesh: working (probe reports '$got_manual' and '$got_auto' as asked)"
    fi
fi

probe_log=$(mktemp)
( cd "$BUILD/examples/$PROBE" && env OMP_CANCELLATION=TRUE OMP_PROC_BIND=TRUE \
    PSOPT_NLP_METHOD=SQP PSOPT_HESSIAN=exact PSOPT_QP_SOLVER="$QP" \
    ${TRUST_REGION:+PSOPT_TRUST_REGION="$TRUST_REGION"} \
    "$TIMEOUT" 300 "./$PROBE" > "$probe_log" 2>&1 ) || true
 
if ! grep -q "^SQP (" "$probe_log"; then
    tail -20 "$probe_log"
    die "the probe run of '$PROBE' did not use the SQP.
    Either the build lacks -DPSOPT_ALLOW_ENV_OVERRIDES=ON, or it lacks -DWITH_SQP=ON.
    Rebuild with both before running this script -- without them every SQP row below
    would silently be an IPOPT result."
fi
info "environment overrides: working (the probe ran under the SQP)"
 
# The banner above is printed before the first subproblem is attempted, so seeing it
# proves only that the SQP started. It does not prove that the QP backend exists. An
# earlier version of this script stopped there, and a whole sweep came back with every
# SQP row failed in under four seconds because the requested backend had never been
# built: the plugin is loaded by dlopen at the first subproblem, and a backend that is
# not there fails every one of them instantly.
if grep -q "could not load the QP backend plugin" "$probe_log"; then
    grep -m1 -A2 "could not load the QP backend plugin" "$probe_log"
    die "the QP backend '$QP' is not available in this build.
    Rebuild PSOPT with -DWITH_$(echo "$QP" | tr '[:lower:]' '[:upper:]')=ON, or pass
    --qp with a backend that is built. GALAHAD is the one the default build enables."
fi
 
# And a backend that loads still has to solve something. A trivial problem it cannot
# finish means the sweep below would measure the backend's installation, not the SQP.
if ! grep -q "Optimal solution found" "$probe_log"; then
    tail -25 "$probe_log"
    die "the QP backend '$QP' loaded but did not solve the probe problem '$PROBE'.
    Every SQP row in the sweep would be measuring that, so the run is stopped here."
fi
info "QP backend '$QP': working (the probe solved under it)"
 
if grep -qi "cannot run in this environment\|OMP_CANCELLATION" "$probe_log"; then
    die "GALAHAD reports that it cannot run. Export OMP_CANCELLATION=TRUE and
    OMP_PROC_BIND=TRUE before starting -- the OpenMP runtime reads them once, at
    process start, so this script setting them per-run is not sufficient if something
    in your environment overrides them."
fi
info "GALAHAD environment: working"
rm -f "$probe_log"
 
mkdir -p "$OUTDIR/logs"
CSV="$OUTDIR/results.csv"
echo "example,solver,status,meshes,meshes_solved,meshes_acceptable,sqp_iterations,objective,seconds,exit_code" > "$CSV"
 
# ---------------------------------------------------------------------------------
# The examples
# ---------------------------------------------------------------------------------
if [ -n "$ONLY" ]; then
    EXAMPLES=$(echo "$ONLY" | tr ',' ' ')
else
    EXAMPLES=""
    for d in "$BUILD"/examples/*/; do
        e=$(basename "$d")
        [ -x "$d/$e" ] && EXAMPLES="$EXAMPLES $e"
    done
fi
 
if [ -n "$SKIP" ]; then
    KEPT=""
    for e in $EXAMPLES; do
        drop=0
        for k in $(echo "$SKIP" | tr ',' ' '); do
            [ "$e" = "$k" ] && drop=1
        done
        [ "$drop" -eq 0 ] && KEPT="$KEPT $e"
    done
    info "skipping on request: $(echo "$SKIP" | tr ',' ' ')"
    EXAMPLES="$KEPT"
fi
 
N=$(echo "$EXAMPLES" | wc -w | tr -d ' ')
 
say "Running $N examples under two solvers, limit ${LIMIT}s each"
info "SQP QP backend: $QP"
info "SQP QP backend: $QP"
info "output: $OUTDIR"
info "this will take a while; results.csv is written as it goes"
 
run_one() {
    local e="$1" solver="$2"
    local d="$BUILD/examples/$e"
    local log="$OUTDIR/logs/${e}_${solver}.log"
    local env_args=""
 
    if [ "$solver" = "SQP" ]; then
        env_args="PSOPT_NLP_METHOD=SQP PSOPT_HESSIAN=exact PSOPT_QP_SOLVER=$QP"
        [ -n "$TRUST_REGION" ] && env_args="$env_args PSOPT_TRUST_REGION=$TRUST_REGION"
    else
        env_args="PSOPT_NLP_METHOD=IPOPT"
    fi
    [ "$FIXED_MESH" = "1" ] && env_args="$env_args PSOPT_MESH_REFINEMENT=manual"
 
    local start end secs rc
    start=$(date +%s)
    # env_args is deliberately unquoted: it carries several NAME=VALUE assignments
    # that env must see as separate arguments.
    # shellcheck disable=SC2086
    ( cd "$d" && env OMP_CANCELLATION=TRUE OMP_PROC_BIND=TRUE $env_args \
        "$TIMEOUT" -s KILL "$LIMIT" "./$e" > "$log" 2>&1 ) && rc=0 || rc=$?
    end=$(date +%s)
    secs=$((end - start))
 
    # Classification. An example counts as solved only if it returned an objective and
    # every mesh refinement reached optimality; anything else is reported as what it was
    # rather than folded into a single "failed".
    local meshes solved accept iters obj status
    meshes=$(grep -c "mesh refinement iteration:" "$log" || true)
    # Both solvers can stop at an acceptable level rather than at the requested tolerance,
    # and both mean the same thing by it: feasible, and dual-feasible to a looser bound.
    # They say it differently. The SQP prints "Optimal solution found to acceptable
    # tolerance", which matches the pattern below; IPOPT prints "EXIT: Solved To Acceptable
    # Level", which did not. So the identical outcome counted as solved for the SQP and as
    # partial for IPOPT -- a bias in the SQP's favour, in every sweep run so far. It is what
    # made examples/lqr_radau read as the one example the SQP solves and IPOPT does not,
    # when IPOPT in fact returns the same objective to seven figures with a dual
    # infeasibility of 8e-09. Count both, and report how often it happened.
    solved=$(grep -cE "Optimal solution found|EXIT: Optimal Solution Found|EXIT: Solved To Acceptable Level" "$log" || true)
    accept=$(grep -cE "Optimal solution found to acceptable tolerance|EXIT: Solved To Acceptable Level" "$log" || true)
    iters=$(grep -oE "SQP finished after [0-9]+" "$log" | grep -oE "[0-9]+" | paste -sd+ - | bc 2>/dev/null || true)
    obj=$(grep -E "^Optimal .unscaled. cost function value:" "$log" | tail -1 | grep -oE "[-0-9.e+]+$" || true)
 
    if   [ "$rc" = "137" ] || [ "$rc" = "124" ]; then status="timeout"
    # A self-checking study rather than a single solve: it ran to completion, printed no
    # objective, and has no mesh-refinement trace at all. Recognising these by the word
    # PASS, as this did, missed examples/fuller -- which compares collocation against an
    # integrated-residual formulation across several residual boxes and prints its own
    # table without that word -- so fuller was recorded as a failure for both solvers in
    # every sweep, and read as IPOPT being unable to solve it. Nothing that solves has
    # zero mesh-refinement iterations, and nothing that fails exits zero with no output,
    # so the three conditions together are exact on the shipped set.
    elif [ "$rc" = "0" ] && [ "$meshes" = "0" ] && [ -z "$obj" ]; then status="pass_selfcheck"
    elif [ -n "$obj" ] && [ "$meshes" -gt 0 ] && [ "$solved" = "$meshes" ]; then status="solved"
    elif [ -n "$obj" ]; then status="partial"
    else status="failed"
    fi
 
    printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
        "$e" "$solver" "$status" "$meshes" "$solved" "$accept" "${iters:-}" "${obj:-}" "$secs" "$rc" >> "$CSV"
    printf "    %-22s %-6s %-14s %6ss  %s\n" "$e" "$solver" "$status" "$secs" "${obj:-}"
}
 
for e in $EXAMPLES; do
    for sv in $SOLVERS; do run_one "$e" "$sv"; done
done
 
# ---------------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------------
say "Summary"
MESH_MODE="automatic refinement"
[ "$FIXED_MESH" = "1" ] && MESH_MODE="prescribed mesh (--fixed-mesh)"

python3 - "$CSV" "$OUTDIR/results.md" "$LIMIT" "$QP" "$MESH_MODE" <<'PYEOF'
import csv, sys, collections
csv_path, md_path, limit = sys.argv[1], sys.argv[2], sys.argv[3]
qp, mesh_mode = sys.argv[4], sys.argv[5]
rows = list(csv.DictReader(open(csv_path)))
by = collections.defaultdict(dict)
for r in rows:
    by[r['example']][r['solver']] = r
 
def ok(r): return r and r['status'] in ('solved', 'pass_selfcheck')
 
# Some examples are self-checking studies rather than a single optimal control solve:
# they run several configurations and print their own verdict, so they have no single
# objective and no mesh-refinement trace. They are named, but kept out of the counts and
# the statistics, where they would compare nothing against nothing.
def is_study(e):
    for sv in ('IPOPT', 'SQP'):
        r = by[e].get(sv)
        if not r: continue
        if r['status'] == 'pass_selfcheck': return True
        if r['status'] == 'failed' and not r['objective'] and r['meshes'] == '0': return True
    return False
 
ex = sorted(by)
studies = [e for e in ex if is_study(e)]
cmp_ex = [e for e in ex if e not in studies]
ni = sum(1 for e in cmp_ex if ok(by[e].get('IPOPT')))
ns = sum(1 for e in cmp_ex if ok(by[e].get('SQP')))
both = [e for e in cmp_ex if ok(by[e].get('IPOPT')) and ok(by[e].get('SQP'))]
 
out = []
out.append(f"# PSOPT: IPOPT against the built-in SQP\n")
out.append(f"{len(cmp_ex)} optimal control problems, time limit {limit} s per run.\n")
# The run's own conditions, recorded here because a sweep is meaningless without them and
# two sweeps that differ in them must never be put in the same table.
out.append(f"QP backend: **{qp}**. Mesh: **{mesh_mode}**.\n")

# How often each solver stopped at an acceptable level rather than at the requested
# tolerance. Both are counted as solved, because both mean the same thing by it, but the
# counts belong in the report: a column of solves that are mostly acceptable-level is a
# different claim from one that is not.
acc = {}
for sv in ('IPOPT', 'SQP'):
    acc[sv] = sum(1 for e in cmp_ex
                  if (by[e].get(sv) or {}).get('meshes_acceptable', '0') not in ('', '0'))
if acc['IPOPT'] or acc['SQP']:
    out.append(f"Runs that stopped at an acceptable level rather than the requested "
               f"tolerance, counted as solved for both: IPOPT {acc['IPOPT']}, "
               f"SQP {acc['SQP']}.\n")
if studies:
    out.append(f"({len(studies)} further examples are self-checking studies rather than a single "
               f"solve, and are excluded from everything below: " + ", ".join(studies) + ".)\n")
out.append(f"| | IPOPT | SQP |\n|---|---|---|")
out.append(f"| solved | **{ni}** | **{ns}** |")
for st in ('timeout', 'partial', 'failed'):
    out.append(f"| {st} | {sum(1 for e in cmp_ex if (by[e].get('IPOPT') or {}).get('status')==st)} "
               f"| {sum(1 for e in cmp_ex if (by[e].get('SQP') or {}).get('status')==st)} |")
out.append("")
 
only_i = [e for e in cmp_ex if ok(by[e].get('IPOPT')) and not ok(by[e].get('SQP'))]
only_s = [e for e in cmp_ex if ok(by[e].get('SQP')) and not ok(by[e].get('IPOPT'))]
out.append(f"**IPOPT solves, the SQP does not ({len(only_i)}):** " + (", ".join(only_i) or "none") + "\n")
out.append(f"**The SQP solves, IPOPT does not ({len(only_s)}):** " + (", ".join(only_s) or "none") + "\n")
 
# agreement and speed, over the examples both solve
dis, ratios = [], []
tot_i = tot_s = 0.0
for e in both:
    try:
        a, b = float(by[e]['IPOPT']['objective']), float(by[e]['SQP']['objective'])
        # A relative comparison of two objectives that are both numerically zero says
        # nothing: 2.4e-18 against 4.4e-13 is a relative difference of 1, and reporting
        # it as disagreement invites exactly the wrong conclusion. Require the larger of
        # the two to be meaningfully non-zero before comparing at all.
        if max(abs(a), abs(b)) > 1e-8:
            rel = abs(a-b)/max(abs(a), abs(b))
            if rel > 1e-4: dis.append((e, a, b, rel))
    except (ValueError, KeyError): pass
    try:
        ti, ts = float(by[e]['IPOPT']['seconds']), float(by[e]['SQP']['seconds'])
        # Per-example ratios are only meaningful where the denominator is well clear of
        # the one-second timing granularity; below that a ratio is mostly quantisation.
        if ti >= 10.0: ratios.append((ts/ti, e))
        tot_i += ti; tot_s += ts
    except (ValueError, KeyError): pass
 
out.append(f"Both solve {len(both)} examples.\n")
# Aggregate wall clock is the honest headline: most examples finish in a second or
# two, so a median of per-example ratios is dominated by 1s/1s = 1.0 and says almost
# nothing. Report the total, and name the individual gaps that make it up.
if tot_i > 0:
    out.append(f"Total wall clock over those {len(both)}: IPOPT **{tot_i:.0f} s**, "
               f"SQP **{tot_s:.0f} s** (**{tot_s/tot_i:.1f}x**).\n")
# No median of per-example ratios is reported. Most of these problems finish in a second
# or two, and at one-second granularity such a ratio is quantisation, not a measurement.
# Count instead the runs the SQP stretched by a wide margin.
wide = sorted(((ts/ti, e) for e in both
               for ti, ts in [(float(by[e]['IPOPT']['seconds']), float(by[e]['SQP']['seconds']))]
               if ti >= 2.0 and ts > 5.0*ti), reverse=True)
if wide:
    out.append(f"The SQP took more than five times as long on {len(wide)} of the runs where IPOPT "
               f"needed at least two seconds: " +
               ", ".join(f"{e} ({r:.0f}x)" for r, e in wide) + ".\n")
gaps = sorted(((float(by[e]['SQP']['seconds']) - float(by[e]['IPOPT']['seconds']), e) for e in both),
              reverse=True)[:5]
gaps = [(g, e) for g, e in gaps if g > 0]
if gaps:
    out.append("Largest absolute gaps: " +
               ", ".join(f"{e} (+{g:.0f} s)" for g, e in gaps) + ".\n")
if dis:
    out.append(f"Objectives differing by more than 1e-4 relative ({len(dis)}):\n")
    out.append("| example | IPOPT | SQP | relative |\n|---|---|---|---|")
    for e,a,b,r in dis: out.append(f"| {e} | {a:.6e} | {b:.6e} | {r:.1e} |")
elif both:
    out.append("Every objective agrees to better than 1e-4 relative.\n")
else:
    # Saying "every objective agrees" when nothing solved is worse than saying nothing:
    # it reads as a clean bill of health for a run that produced no results at all.
    out.append("**No example was solved by both solvers, so there is nothing to compare.** "
               "A whole column of failures usually means the run was misconfigured rather "
               "than that the solver cannot do these problems -- check the QP backend "
               "first.\n")
 
text = "\n".join(out) + "\n"
open(md_path,'w').write(text)
print(text)
PYEOF
 
say "Done"
info "send back: $OUTDIR/results.csv and $OUTDIR/results.md"
info "the logs are in $OUTDIR/logs if anything needs explaining"