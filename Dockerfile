# Use latest ArchLinux as the base
FROM archlinux:latest

# 1. Update system packages
RUN pacman -Syu --noconfirm

# 2. Install official repo dependencies for PSOPT
RUN pacman -S --noconfirm \
base-devel \
cmake \
gnuplot \
eigen \
boost \
blas \
lapack \
git

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

# 5. Use yay to install AUR packages: coin-or-ipopt, colpack, adol-c
#    - If these have dependencies also from the AUR, yay will build them automatically.
RUN yay -S --noconfirm coin-or-ipopt colpack adol-c

# Switch back to root for system-wide tasks
USER root
WORKDIR /opt/psopt

# 6. Clone PSOPT (or COPY your local code instead)
RUN git clone --tags https://github.com/PSOPT/psopt.git psopt

# 7. Build & install PSOPT under /usr
WORKDIR /opt/psopt/psopt
RUN mkdir build && cd build && \
    cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_EXAMPLES=ON \
      -DHEADLESS=ON \
      -DCMAKE_INSTALL_PREFIX=/usr && \      # <<< important
    make -j$(nproc) && \
    make install && \
    ldconfig                                     # refresh loader cache
# 8. Set default command
CMD ["/bin/bash"]
