PSOPT installation on Ubuntu 22.04
=====


With [Ubuntu 22.04](https://releases.ubuntu.com/22.04/), a runtime error related to the `adolc` library is currently reported when executing PSOPT code if you follow the instructions that work for Ubuntu 20.04. Here are the specific instructions that are needed to install PSOPT under Ubuntu 22.04.


First, you should run the following from a terminal window to install some packages that are required:

1. `sudo apt-get install git`
2. `sudo apt-get install cmake`
3. `sudo apt-get install gfortran`
4. `sudo apt-get install g++`
5. `sudo apt-get install libboost-dev`
6. `sudo apt-get install libboost-system-dev`
7. `sudo apt-get install coinor-libipopt-dev`
8. `sudo apt-get install gnuplot`
9. `sudo apt-get install libeigen3-dev`
10. `sudo apt-get install libblas-dev`
11. `sudo apt-get install liblapack-dev`

Second, you should run the following commands to download, compile and install adolc and ColPack.

1. `cd $HOME/Downloads`
2. `wget --continue www.coin-or.org/download/source/ADOL-C/ADOL-C-2.6.3.tgz`
3. `cd $HOME`
4. `tar zxvf ./Downloads/ADOL-C-2.6.3.tgz`
5. `cd $HOME/ADOL-C-2.6.3`
6. `mkdir ./ThirdParty`
7. `cd ./ThirdParty`
8. `wget --continue http://archive.ubuntu.com/ubuntu/pool/universe/c/colpack/colpack_1.0.10.orig.tar.gz`
9. `tar zxvf colpack_1.0.10.orig.tar.gz`
10. `mv ColPack-1.0.10 ColPack`
11. `cd ColPack`
12. `./autoconf.sh`
13. `make`
14. `sudo make install`
15. `sudo cp -P ./build/lib/libCol* /usr/lib`
16. `cd $HOME/ADOL-C-2.6.3`
17. `./configure --enable-sparse --with-colpack=$HOME/ADOL-C-2.6.3/ThirdParty/ColPack/build`
18. `make`
19. `make install`
20. `sudo cp -P $HOME/adolc_base/lib64/lib* /usr/lib`
21. `sudo cp -r $HOME/adolc_base/include/* /usr/include/`
22. `cd $HOME`

Then, you can run the following commands to download, compile and install PSOPT.

1. `git clone https://github.com/PSOPT/psopt.git`
2. `cd psopt; mkdir build; cd build`
3. `cmake -DBUILD_EXAMPLES=ON ..`
4. `make`
5. `sudo make install`

You should create a file called adolc.pc with the following content and copy it to the folder /usr/lib/pkgconfig

**/usr/lib/pkgconfig/adolc.pc**

prefix=/usr/
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include
Name: adolc
Version: 2.6.3
Description: Algorithmic Differentiation Library for C/C++
Requires:
Libs: -L${libdir} -ladolc -Wl,-rpath,${libdir}
-L${libdir} -lColPack -Wl,-rpath,${libdir}
Cflags: -I${includedir}
