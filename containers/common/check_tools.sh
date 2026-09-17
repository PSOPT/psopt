#!/usr/bin/env bash
#
# The command-line tools every later step assumes, checked before anything uses one.
#
# This exists because of what happened on both openSUSE images. They were missing
# findutils and diffutils, which nothing noticed, because of how those tools are called:
#
#     find "$PREFIX" -name osqp-config.cmake -print -quit 2>/dev/null || true
#
# The redirection is there so that a prefix with nothing in it is an empty answer rather
# than an error. It cannot tell "no such file" apart from "no such program", so an absent
# find returned "not installed" for every file it was asked about. The diagnosis that came
# back was self-contradictory -- the install log showed cmake writing
# /opt/qp/lib64/cmake/osqp/osqp-config.cmake, and the check below it reported that no such
# file existed anywhere under /opt/qp -- and two rounds of work went into the wrong
# question, which was why cmake had not installed something it plainly had.
#
# The same absence had already passed through IPOPT's configure, which reported "cmp:
# command not found" and "xargs: command not found" eleven times and carried on to produce
# an IPOPT that then built and passed its tests. That is the more alarming half: a source
# build silently taking a different path because a tool it probes for is missing.
#
# So the tools are named, once, and a missing one is reported as a missing tool by name
# before any of them is used. A dependency stanza that forgets one now fails in the
# Dockerfile with a sentence saying which, rather than three steps later as a phantom
# defect in something else.
#
set -euo pipefail

say() { printf '\n\033[1m== %s\033[0m\n' "$*"; }

# Each entry is "command:required|advisory:what needs it", so a failure says why the tool
# is wanted rather than leaving the reader to guess which step will break.
#
# A required tool stops the image. An advisory one is reported and allowed: file and which
# are probed by configure scripts that fall back when they are absent, and neither carries
# a step of its own, so a missing one is worth seeing in the log and is not worth failing
# an otherwise working distribution over -- particularly which, which Debian has been
# moving out of its essential set. If an advisory tool ever turns out to matter, this list
# is where to promote it, and the promotion will be a measurement rather than a guess.
TOOLS="
find:required:the QP backend installer scripts locate what they installed with it
xargs:required:configure scripts in the IPOPT source build probe for it
cmp:required:configure scripts in the IPOPT source build use it
diff:required:configure scripts in the IPOPT source build use it
sed:required:the shared scripts and several configure scripts
awk:required:the shared scripts read versions with it
grep:required:everything
tar:required:unpacking source archives
gzip:required:unpacking source archives
patch:required:coinbrew applies patches with it
curl:required:rustup, where Clarabel needs a Rust toolchain
git:required:cloning CppAD, OSQP, PIQP, Clarabel and GALAHAD
make:required:building IPOPT
cmake:required:building everything else
pkg-config:required:finding IPOPT and Eigen
python3:required:running the example driver
file:advisory:the shared scripts identify binaries with it
which:advisory:configure scripts in the IPOPT source build probe for it
"

say "Command-line tools"

missing=""
while IFS=: read -r tool need why; do
    [ -n "${tool}" ] || continue
    if command -v "${tool}" >/dev/null 2>&1; then
        printf '  %-12s %s\n' "${tool}" "$(command -v "${tool}")"
    elif [ "${need}" = "advisory" ]; then
        printf '  %-12s absent, and not required -- %s\n' "${tool}" "${why}"
    else
        printf '  %-12s MISSING -- %s\n' "${tool}" "${why}"
        missing="${missing} ${tool}"
    fi
done <<EOF
${TOOLS}
EOF

if [ -n "${missing}" ]; then
    echo
    echo "These tools are not installed:${missing}"
    echo
    echo "They are the Dockerfile's job, like every other distribution difference. Add the"
    echo "package that carries each one to this image's dependency stanza. find and xargs"
    echo "are findutils; cmp and diff are diffutils; both are absent from the openSUSE"
    echo "base images and present in every other image in this matrix, which is how this"
    echo "check came to exist."
    exit 1
fi
