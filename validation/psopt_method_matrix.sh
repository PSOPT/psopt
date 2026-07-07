#!/bin/bash
# Matrix: every PSOPT example x {native IPOPT, CASADI-ipopt, CASADI-sqpmethod}.
# For each example+method: sed the nlp_method line, rebuild that example target in
# the superbuild (GNU+CASADI) tree, run it, capture solved/objective/time.
set -uo pipefail
export PATH=/opt/local/bin:/usr/bin:/bin
export SDKROOT="$(xcrun --show-sdk-path 2>/dev/null)"
PSOPT=~/src/psopt
EPB=$PSOPT/build-sb/ep_psopt-prefix/src/ep_psopt-build
DEPS=$PSOPT/build-sb/_deps/install
export DYLD_LIBRARY_PATH=$DEPS/lib:$DEPS/lib/casadi:$DEPS/lib64
OUT=/opt/psopt/method-matrix; mkdir -p "$OUT"
RES=$OUT/RESULTS.tsv
printf "example\tmethod\tstatus\tsecs\tobjective\n" > "$RES"

# timeouts (s): sqpmethod is FD-slow, cap tighter
declare -A TMO=( [IPOPT]=120 [CASADI-ipopt]=90 [CASADI-sqpmethod]=40 )

sed_for () { # method -> replacement text for the nlp_method assignment
  case "$1" in
    IPOPT)            echo 'algorithm.nlp_method = "IPOPT";';;
    CASADI-ipopt)     echo 'algorithm.nlp_method = "CASADI"; algorithm.casadi_solver = "ipopt";';;
    CASADI-sqpmethod) echo 'algorithm.nlp_method = "CASADI"; algorithm.casadi_solver = "sqpmethod";';;
  esac
}

runcase () { # example method srcfile
  local ex=$1 m=$2 src=$3
  cp "$src" "$src.mbak"
  perl -0pi -e 's/algorithm\.nlp_method\s*=\s*"IPOPT";/'"$(sed_for "$m")"'/' "$src"
  if ! cmake --build "$EPB" --target "$ex" >/tmp/mm_build.log 2>&1; then
    printf "%s\t%s\tBUILD-FAIL\t\t\n" "$ex" "$m" | tee -a "$RES"; mv "$src.mbak" "$src"; return; fi
  mv "$src.mbak" "$src"
  local bin; bin=$(find "$EPB/examples/$ex" -name "$ex" -type f -perm +111 2>/dev/null | head -1)
  [ -x "$bin" ] || { printf "%s\t%s\tNO-BIN\t\t\n" "$ex" "$m" | tee -a "$RES"; return; }
  local d; d=$(dirname "$bin"); local t0 t1
  t0=$(date +%s)
  ( cd "$d" && /opt/local/bin/gtimeout ${TMO[$m]} "./$ex" >/tmp/mm_run.log 2>&1 ); local rc=$?
  t1=$(date +%s)
  local obj st
  obj=$(grep -m1 'cost function value' /tmp/mm_run.log | grep -oE '[-0-9.eE+]+' | tail -1)
  if grep -qiE 'problem has been solved|Optimal Solution Found' /tmp/mm_run.log; then st=SOLVED
  elif [ $rc -eq 124 ]; then st=TIMEOUT; else st="rc=$rc"; fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$ex" "$m" "$st" "$((t1-t0))" "${obj:-}" | tee -a "$RES"
}

for src in $PSOPT/examples/*/[a-z]*.cxx; do
  ex=$(basename "$(dirname "$src")")
  # only the file that actually sets nlp_method (skip helper .cxx)
  grep -q 'nlp_method' "$src" || continue
  for m in IPOPT CASADI-ipopt CASADI-sqpmethod; do runcase "$ex" "$m" "$src"; done
done
echo "### MATRIX DONE $(date)"
