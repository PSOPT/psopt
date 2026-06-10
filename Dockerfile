# Use latest ArchLinux as the base
FROM archlinux:latest

# 1. Update system packages
RUN pacman -Syu --noconfirm

# 2. Install official repo dependencies for PSOPT
#    (python is used by the CI example-runner script, .github/scripts/run_examples.py)
RUN pacman -S --noconfirm \
base-devel \
cmake \
gnuplot \
eigen \
boost \
blas \
lapack \
git \
python

# 3. Create a non-root user (required to build AUR packages cleanly)
RUN useradd -m builduser && echo "builduser ALL=(ALL) NOPASSWD: ALL" >> /etc/sudoers

# Switch to the builduser
USER builduser
WORKDIR /home/builduser

# 4. Install 'yay' (an AUR helper) from the AUR
#    - This step clones the yay repository and uses makepkg to build and install.
RUN git clone https://aur.archlinux.org/yay.git && \
cd yay && \
makepkg -si --noconfirm && \
cd .. && rm -rf yay

# 5. Use yay to install the IPOPT NLP solver from the AUR.
#    ADOL-C and ColPack are no longer needed: PSOPT now uses CppAD for automatic
#    differentiation (in-memory tapes). CppAD does its own sparsity colouring, so ColPack
#    -- which was only ever an ADOL-C dependency -- is dropped as well.
RUN yay -S --noconfirm coin-or-ipopt

# Switch back to root for system-wide tasks
USER root

# 6. Build and install CppAD from source under /usr.
#    PSOPT links BOTH the CppAD headers (cppad/cppad.hpp) and the compiled libcppad_lib,
#    so a full source install is the reliable way to provide both. Building from source
#    (rather than the AUR 'cppad' package) guarantees libcppad_lib is present and installed
#    to a standard location. NOTE: CppAD renamed its install-prefix variable to
#    'cppad_prefix' (passing CMAKE_INSTALL_PREFIX is a hard error). This installs
#    /usr/include/cppad/... and /usr/lib/libcppad_lib.so.
#    (To pin a specific release for reproducibility, add e.g. `--branch 20250000.2` to the clone.)
WORKDIR /opt/cppad
RUN git clone --depth 1 https://github.com/coin-or/CppAD.git src && \
    cd src && mkdir build && cd build && \
    cmake -D cppad_prefix=/usr .. && \
    make -j"$(nproc)" && \
    make install && \
    cd / && rm -rf /opt/cppad

# 7. Clone PSOPT (or COPY your local code instead).
#    This expects the CppAD-based PSOPT (ADOL-C removed; CppAD is the default backend).
WORKDIR /opt/psopt
RUN git clone --tags https://github.com/psopt/psopt.git psopt

# 8. Build & install PSOPT with the CppAD backend under /usr
WORKDIR /opt/psopt/psopt
RUN mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DBUILD_EXAMPLES=ON -DHEADLESS=ON \
             -DCMAKE_INSTALL_PREFIX=/usr \
             -DPSOPT_AD_BACKEND=CPPAD \
             -DCPPAD_INCLUDE_DIR=/usr/include \
             -DCPPAD_LIBRARY=/usr/lib/libcppad_lib.so && \
    make -j"$(nproc)" && \
    make install && \
    ldconfig                                     # refresh loader cache

# 9. Set default command
CMD ["/bin/bash"]
