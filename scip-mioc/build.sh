#!/bin/bash
# Build the SCIP mixed-integer optimal control example against /opt/claude/scip.
set -e
SCIP=/opt/claude/scip
cc -O2 -I$SCIP/include mioc_scip.c -o mioc_scip \
   -L$SCIP/lib -lscip -Wl,-rpath,$SCIP/lib -Wl,-rpath,/opt/local/lib
echo "built mioc_scip"
