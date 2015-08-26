# Script for Ubuntu 14.04, taken from PSOPT manual.
# Install miscellaneous libraries
sudo apt-get -y install g++ gfortran f2c libf2c2-dev libf2c2 libblas-dev libopenblas-base libopenblas-dev libblas3gf libatlas-base-dev liblapack-dev liblapack3gf
cd $HOME
# Download and install Ipopt, Metis and Mumps
wget --continue http://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.3.tgz
tar xzvf Ipopt-3.12.3.tgz
cd $HOME/Ipopt-3.12.3/ThirdParty/Metis
./get.Metis
cd $HOME/Ipopt-3.12.3/ThirdParty/Mumps
./get.Mumps
cd $HOME/Ipopt-3.12.3
./configure --enable-static coin_skip_warn_cxxflags=yes
make -j
make install
# Download and install ADOLC and ColPack
cd $HOME
wget --continue www.coin-or.org/download/source/ADOL-C/ADOL-C-2.5.2.tgz
tar zxvf ADOL-C-2.5.2.tgz
cd $HOME/ADOL-C-2.5.2
mkdir ./ThirdParty
cd ./ThirdParty
wget --continue http://cscapes.cs.purdue.edu/download/ColPack/ColPack-1.0.9.tar.gz
tar zxvf ColPack-1.0.9.tar.gz
mv ColPack-1.0.9 ColPack
cd ColPack
./configure
make
sudo make install
sudo cp /usr/local/lib/libCol* /usr/lib
cd $HOME/ADOL-C-2.5.2
./configure --enable-sparse --with-colpack=$HOME/ADOL-C-2.5.2/ThirdParty/ColPack
make
make install
sudo cp $HOME/adolc_base/lib64/*.a /usr/lib
sudo cp -r $HOME/adolc_base/include/* /usr/include/
# Download and install PDFlib
cd $HOME/Downloads
wget --continue http://www.pdflib.com/binaries/PDFlib/705/PDFlib-Lite-7.0.5p3.tar.gz
tar zxvf PDFlib-Lite-7.0.5p3.tar.gz
cd PDFlib-Lite-7.0.5p3 $ ./configure
make; sudo make install
sudo ldconfig
# Download and install GNUplot
cd $HOME/Downloads
wget --continue http://sourceforge.net/projects/gnuplot/files/gnuplot/4.2.2/gnuplot-4.2.2.tar.gz/download
mv download gnuplot-4.2.2.tar.gz
tar zxvf gnuplot-4.2.2.tar.gz
cd gnuplot-4.2.2
./configure -with-readline=gnu -without-tutorial
make;sudo make install
# Download and extract PSOPT
cd $HOME
wget --continue https://github.com/PSOPT/psopt/archive/master.zip
unzip master.zip
cd $HOME/psopt-master
# Download and extract SuiteSparse
wget --continue http://faculty.cse.tamu.edu/davis/SuiteSparse/SuiteSparse-4.4.3.tar.gz
tar zxvf SuiteSparse-4.4.3.tar.gz
cd $HOME/psopt-master
# Download and extract LUSOL
wget --continue http://www.stanford.edu/group/SOL/software/lusol/lusol.zip
unzip lusol.zip
cd $HOME/psopt-master
# Compile SuiteSparse, LUSOL, dmatrix and PSOPT
make all
