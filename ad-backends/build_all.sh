#!/bin/bash
# Build all five AD-backend drivers for the PSOPT AD shootout and run the sweep.
# Operator-overloading backends use GNU g++-mp-15 (libstdc++, ABI-matched to the
# superbuilt ADOL-C 2.7.2); Enzyme uses the llvm-20 clang + ClangEnzyme plugin.
#
# Prereqs (adjust paths for your machine):
#   ADOLC   : superbuilt ADOL-C 2.7.2  (PSOPT superbuild _deps/install)
#   CPPAD   : coin-or/CppAD installed  (+ libcppad_lib)
#   CPPADCG : joaoleal/CppADCodeGen installed (generates *_c.hpp)
#   SACADO  : Trilinos/Sacado built (-DTrilinos_ENABLE_Sacado, C++20, Kokkos OFF)
#   ENZYME  : EnzymeAD/Enzyme built against llvm-20 -> ClangEnzyme-20.dylib
set -euo pipefail
export PATH=/opt/local/bin:/usr/bin:/bin
export SDKROOT="$(xcrun --show-sdk-path)"

GXX=/opt/local/bin/g++-mp-15
CLANGXX=/opt/local/libexec/llvm-20/bin/clang++
HERE=/opt/claude/ad-shootout
PFX=$HERE/install
DEPS=$HOME/src/psopt/build-sb/_deps/install     # superbuilt ADOL-C
EIG=/opt/local/include/eigen3
CPPADLIB=$HERE/src/CppAD/build/cppad_lib
ENZ=$HERE/src/Enzyme/build/Enzyme/ClangEnzyme-20.dylib

echo ">> ADOL-C";       $GXX -O2 -std=c++17 -I. -I"$DEPS/include" driver_adolc.cpp   -o bin_adolc   -L"$DEPS/lib64" -L"$DEPS/lib" -ladolc     -Wl,-rpath,"$DEPS/lib64" -Wl,-rpath,"$DEPS/lib"
echo ">> CppAD";        $GXX -O2 -std=c++17 -I. -I"$PFX/include" -I"$EIG" driver_cppad.cpp -o bin_cppad -L"$CPPADLIB" -lcppad_lib -Wl,-rpath,"$CPPADLIB"
echo ">> CppADCodeGen"; $GXX -O2 -std=c++17 -I. -I"$PFX/include" -I"$EIG" driver_cppadcg.cpp -o bin_cppadcg -L"$CPPADLIB" -lcppad_lib -Wl,-rpath,"$CPPADLIB" -ldl
echo ">> Sacado";       $GXX -O2 -std=c++20 -I. -I"$PFX/include" driver_sacado.cpp -o bin_sacado -L"$PFX/lib" -lsacado -Wl,-rpath,"$PFX/lib"
echo ">> Enzyme";       $CLANGXX -O2 -std=c++17 -I. -fpass-plugin="$ENZ" -Xclang -load -Xclang "$ENZ" driver_enzyme.cpp -o bin_enzyme

export DYLD_LIBRARY_PATH="$DEPS/lib64:$DEPS/lib:$CPPADLIB:$PFX/lib"
echo; echo "== sweep =="
for n in 10 20 50 100 200; do
  for b in adolc cppad cppadcg sacado enzyme; do ./bin_$b "$n" 1000; done
done
