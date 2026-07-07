# Building PSOPT in Visual Studio (Windows, MSVC)

> **Status: Supported on MSVC.** PSOPT has historically been built with
> GCC/Clang on Linux and macOS, but it now fully supports building with MSVC.
> The C++-only from-source **superbuild (`-DPSOPT_SUPERBUILD=ON`) is supported on
> MSVC** and automates building ADOL-C and CasADi using MSBuild/CMake and the vcpkg toolchain.
> Alternatively, you can use vcpkg + a manual ADOL-C build.

## 1. Prerequisites (Visual Studio 18)

In the VS Installer, under **Desktop development with C++**, ensure:
- MSVC v14x toolset, Windows SDK
- **C++ CMake tools for Windows** (gives CMake + Ninja; VS reads `CMakePresets.json`)

Install **vcpkg** (if not already):
```pwsh
git clone https://github.com/microsoft/vcpkg C:\vcpkg
C:\vcpkg\bootstrap-vcpkg.bat
setx VCPKG_ROOT C:\vcpkg     # the preset reads $env{VCPKG_ROOT}
```

## 2. Dependencies

| Dependency | How on Windows/MSVC |
|---|---|
| **Eigen3** (header) | `vcpkg install eigen3:x64-windows` (or CMake auto-fetches it) |
| **Boost** | `vcpkg install boost-system boost-filesystem:x64-windows` |
| **BLAS/LAPACK** | `vcpkg install openblas:x64-windows` (or use Intel oneMKL) |
| **IPOPT (+MUMPS)** | `vcpkg install ipopt:x64-windows` *(if the port is available)*, **or** download the official IPOPT Windows binaries from the COIN-OR releases and set `IPOPT_ROOT` |
| **ADOL-C** (required) | **Not in vcpkg.** Build it yourself: clone `coin-or/ADOL-C`, open `ADOL-C/MSVisualStudio/v14/adolc.sln`, build `x64/Release`, then set `ADOLC_ROOT` to the install/output dir |

PSOPT's `FindIPOPT`/`Findadolc` modules prefer pkg-config (Unix). On Windows, help
them by pointing at install roots on the CMake command line / in the preset's
`cacheVariables`, e.g. `-DIPOPT_ROOT=C:/ipopt -DADOLC_ROOT=C:/adolc`
(add those to `CMakePresets.json` once your paths are known).

## 3. Build in Visual Studio

1. **File → Open → Folder…** and select the PSOPT source folder.
2. VS detects `CMakePresets.json`. In the toolbar preset dropdown pick
   **"Windows x64 (MSVC, vcpkg)"** (`windows-msvc`).
3. VS runs CMake configure (using the vcpkg toolchain). Fix any missing-dependency
   errors by installing the vcpkg package or setting `*_ROOT`.
4. **Build → Build All** (or `Ctrl+Shift+B`). Binaries land in
   `out/build/windows-msvc/`.

Or from a *Developer PowerShell for VS*:
```pwsh
cmake --preset windows-msvc
cmake --build --preset windows-msvc
```

## 4. The Fortran problem (MSVC builds only C++)

Visual Studio's MSVC compiles **C/C++ only** — it has **no Fortran compiler**. But
IPOPT's default linear solver **MUMPS** (and HSL) is **Fortran**, so a plain MSVC
build cannot produce them. There are two ways forward.

### Option A — supply Fortran + MKL via Intel oneAPI

Install **Intel oneAPI Base Toolkit** (oneMKL) and **HPC Toolkit** (Intel Fortran
`ifx`). Run configure from an *Intel oneAPI command prompt* (it runs `setvars.bat`
so `MKLROOT` and `ifx` are on `PATH`), and use the **`windows-intel`** preset:
- **oneMKL** supplies BLAS/LAPACK and the **PARDISO** linear solver. Point IPOPT at
  MKL PARDISO with `algorithm.ipopt_linear_solver = "pardisomkl"` — this lets you
  build IPOPT **without MUMPS**, minimising the Fortran you must compile.
- **`ifx`** builds whatever Fortran remains (MUMPS if you still want it, or IPOPT's
  Fortran pieces). MSVC stays the C/C++ compiler; `ifx` is only the Fortran one.

### Option B — avoid Fortran entirely with the CasADi backend

PSOPT itself has **no Fortran** (`project(... LANGUAGES C CXX)`); Fortran enters
*only* through IPOPT/MUMPS/HSL. Drive the NLP through CasADi with a **C++** solver
and you need **no Fortran at all**:
- Build **CasADi** (CMake, MSVC-friendly) with a pure-C++ solver — **`fatrop`**, or
  **`sqpmethod` + qpOASES/qrqp** — and *without* the IPOPT/MUMPS plugins.
- Configure PSOPT with the **`windows-casadi-nofortran`** preset
  (`PSOPT_WITH_MUMPS=OFF`, `PSOPT_WITH_CASADI=ON`).
- At run time: `algorithm.nlp_method = "CASADI"; algorithm.casadi_solver = "fatrop";`
  (PSOPT passes the plugin straight to `casadi::nlpsol`).

**Caveats:** PSOPT **still requires ADOL-C** — but that's C++ (no Fortran), so it
doesn't reintroduce the problem. And `sqpmethod` presently runs with
finite-difference derivatives (slow; only practical with exact derivatives — see
`docs/PSOPT_vs_CasADi.md`), so **`fatrop` is the better default** for a Fortran-free
build.

**Summary:** MSVC + Intel MKL/PARDISO + `ifx` = the classic IPOPT path; MSVC +
CasADi(`fatrop`) + MKL/OpenBLAS = a **Fortran-free** path. Both still need ADOL-C.

## Building the dependencies from source (superbuild) on Windows

The `-DPSOPT_SUPERBUILD=ON` from-source dependency build uses **autotools + `make`
+ `gfortran`** (for ColPack/ADOL-C/METIS/MUMPS), so it **cannot run under MSVC**.
On Windows it requires the **MinGW-w64 GNU toolchain under MSYS2** — a separate path
from the MSVC presets above — and it *also* gives you the Fortran linear solvers
(MUMPS) without needing Intel Fortran.

1. Install **MSYS2** (https://www.msys2.org), open a **MINGW64** shell, install the
   toolchain:
   ```bash
   pacman -S --needed mingw-w64-x86_64-gcc mingw-w64-x86_64-gcc-fortran \
     mingw-w64-x86_64-cmake mingw-w64-x86_64-openblas \
     autoconf automake libtool make patch git gettext-devel pkgconf wget tar diffutils
   ```
2. Configure + build with the superbuild preset (from the MINGW64 shell):
   ```bash
   cmake --preset windows-superbuild
   cmake --build --preset windows-superbuild
   ```

**Caveats (experimental, untested):** `cmake/superbuild.cmake` is currently tuned
for macOS (MacPorts compiler names, `/opt/local` BLAS hints, libc++ notes). Expect
to adjust it for MSYS2 — e.g. the OpenBLAS `find_library` HINTS, and the COIN-OR
ThirdParty `get.Metis`/`get.Mumps` download scripts need a working `wget`/`curl`.
PSOPT source may also need minor MinGW portability fixes. If you get it working,
please upstream the superbuild changes.

### MSVC C++-only from-source superbuild

If you prefer to automate building Eigen, ADOL-C, and CasADi using MSVC, use the `windows-msvc-superbuild` preset:
1. Ensure **vcpkg** is installed and `VCPKG_ROOT` environment variable is set to the installation folder (e.g. `C:\vcpkg`).
2. Run from a **Developer PowerShell for VS 18**:
   ```powershell
   $env:VCPKG_ROOT = "C:\vcpkg"   # or your vcpkg installation path
   cmake --preset windows-msvc-superbuild
   cmake --build --preset windows-msvc-superbuild
   ```
This automatically compiles ADOL-C via MSBuild (disabling Boost pool allocator dependencies) and CasADi via CMake, linking CasADi with OpenBLAS provided by vcpkg.

## Notes / gotchas

- The preset sets **`PSOPT_AUTO_COMPILER=OFF`** so the build uses VS's MSVC rather
  than any MinGW/Clang on `PATH` (PSOPT's auto-selector otherwise prefers GNU).
- **`HEADLESS=ON`** makes gnuplot write PDFs instead of spawning windows (and the
  examples simply skip plotting if gnuplot isn't installed).
- `BUILD_TESTS=OFF` by default here — the GoogleTest suite is fine on MSVC, but get
  a plain build working first, then flip it on.
- Expect possible source edits for MSVC: replace POSIX-only includes, define
  `_USE_MATH_DEFINES` before `<cmath>` for `M_PI`, and address any
  `__attribute__`/VLA usage. Please upstream fixes you make.

If you get a clean MSVC build, that's genuinely new ground for PSOPT — worth a PR.

## Appendix: dependency recipes for the Fortran-free MSVC + CasADi build

Full dependency chain for the `windows-casadi-nofortran` preset (`PSOPT_WITH_IPOPT=OFF`,
`PSOPT_WITH_MUMPS=OFF`, `PSOPT_WITH_CASADI=ON`) — **no Fortran, no IPOPT, no MUMPS, no
ColPack**, only C++ libraries. Run from a **Developer PowerShell for VS 18** (x64).
Untested end-to-end; expect to adjust paths and a few source-level MSVC fixes.

### 0. vcpkg C++ libraries
```powershell
vcpkg install eigen3 boost-system boost-filesystem openblas:x64-windows   # VCPKG_ROOT set
```

### 1. ADOL-C (required — the `adouble` problem API). Not in vcpkg; no CMake — build the VS solution:
```powershell
git clone -b releases/2.7.2 https://github.com/coin-or/ADOL-C C:\deps\ADOL-C
# Open C:\deps\ADOL-C\ADOL-C\MSVisualStudio\v14\adolc.sln in VS 18 -> x64/Release -> Build.
# Output: adolc.lib (+ adolc.dll); headers under C:\deps\ADOL-C\ADOL-C\include
```
Point PSOPT at it (no ColPack — the CasADi path is finite-difference, so ADOL-C's sparse
drivers are unused). Either a clean `include\`+`lib\` prefix via `-DADOLC_ROOT=...`, or:
```powershell
-Dadolc_INCLUDE_DIR=C:\deps\ADOL-C\ADOL-C\include `
-Dadolc_LIBRARY=C:\deps\ADOL-C\ADOL-C\MSVisualStudio\v14\x64\Release\adolc.lib
```
Runtime: `adolc.dll` on `PATH`.

### 2. CasADi with a Fortran-free NLP solver (required). Not in vcpkg; build from source, no IPOPT plugin:
**Option A — `fatrop` (better solver; harder build — its blasfeo backend is assembly/autotools-oriented and may not build under MSVC):**
```powershell
git clone https://github.com/casadi/casadi C:\deps\casadi
cmake -S C:\deps\casadi -B C:\deps\casadi\build -G Ninja `
  -DCMAKE_INSTALL_PREFIX=C:\deps\casadi-install -DWITH_IPOPT=OFF -DWITH_PYTHON=OFF `
  -DBUILD_TESTING=OFF -DWITH_FATROP=ON -DWITH_BUILD_FATROP=ON -DWITH_BUILD_BLASFEO=ON
cmake --build C:\deps\casadi\build --target install
```
**Option B — `sqpmethod` + qpOASES (most reliable MSVC build; qpOASES builds inline, pure C++). Trade-off: sqpmethod runs finite-difference derivatives (slower):**
```powershell
cmake -S C:\deps\casadi -B C:\deps\casadi\build -G Ninja `
  -DCMAKE_INSTALL_PREFIX=C:\deps\casadi-install -DWITH_IPOPT=OFF -DWITH_PYTHON=OFF `
  -DBUILD_TESTING=OFF -DWITH_LAPACK=ON -DWITH_QPOASES=ON -DWITH_BUILD_QPOASES=ON
cmake --build C:\deps\casadi\build --target install
```

### 3. Configure + build PSOPT
```powershell
cmake --preset windows-casadi-nofortran `
  -Dadolc_INCLUDE_DIR=C:\deps\ADOL-C\ADOL-C\include `
  -Dadolc_LIBRARY=C:\deps\ADOL-C\ADOL-C\MSVisualStudio\v14\x64\Release\adolc.lib `
  -DCMAKE_PREFIX_PATH="C:\deps\casadi-install"
cmake --build --preset windows-casadi-nofortran
```

### 4. Select the backend at run time
```cpp
algorithm.nlp_method    = "CASADI";
algorithm.casadi_solver = "fatrop";   // or "sqpmethod" (Option B)
```
Put `adolc.dll` and CasADi's DLLs on `PATH`.

### Dependency summary
| dep | source | Fortran? |
|---|---|---|
| Eigen3, Boost, OpenBLAS | vcpkg | no |
| ADOL-C 2.7.2 | VS `.sln` build, `adolc_INCLUDE_DIR`/`adolc_LIBRARY` | no |
| CasADi (+ fatrop **or** qpOASES) | CMake build, `WITH_IPOPT=OFF` | no |
| ~~IPOPT / MUMPS / ColPack~~ | **not built** (`PSOPT_WITH_IPOPT=OFF`) | — |

### Could this be a superbuild? (partly)
CasADi is CMake, so it drops into the superbuild's `ExternalProject_Add` cleanly (like
the existing `ep_casadi`, with `WITH_IPOPT=OFF`). **ADOL-C 2.7.2 is the snag** — it has no
CMake, only autotools (needs MSYS2, not MSVC) or the VS `.sln`. A true MSVC superbuild would
`ExternalProject_Add` ADOL-C by invoking `msbuild adolc.sln /p:Configuration=Release
/p:Platform=x64` and then copying `adolc.lib`/headers into the install prefix. That's doable
but fiddly and untested — get one manual build working first, then codify it.
