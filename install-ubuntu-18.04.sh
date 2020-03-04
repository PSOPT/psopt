# Script for installing PSOPT on Ubuntu 18.04
# Install miscellaneous libraries

export IPOPT_FILE=Ipopt-3.12.12.tgz
export ADOL_FILE=ADOL-C-2.6.3.tgz
export COLPACK_FILE=colpack_1.0.10.orig.tar.gz
export PDFLIB_FILE=PDFlib-Lite-7.0.5p3.tar.gz
export SUITESPARSE_FILE=SuiteSparse-4.4.3.tar.gz
export LUSOL_FILE=lusol.zip

export MYHOME=`pwd`
echo $MYHOME
mkdir Downloads

sudo apt-get install build-essential pkg-config
sudo apt-get install dh-autoreconf
sudo apt-get -y install g++ gfortran f2c libf2c2-dev libf2c2 libblas-dev libopenblas-base libopenblas-dev libblas3 libatlas-base-dev liblapack-dev liblapack3
echo -e '\n\n\e[33m \e[1mHit enter apt installation was succesfull. Otherwise press Ctrl-C and resolve arisen errors.\e[0m'
read x

####  Download and install Ipopt, Metis and Mumps #### 
cd $MYHOME
wget --continue http://www.coin-or.org/download/source/Ipopt/$IPOPT_FILE
mkdir Ipopt && tar zxvf $IPOPT_FILE -C Ipopt --strip-components 1
mv $IPOPT_FILE $MYHOME/Downloads
cd $MYHOME/Ipopt/ThirdParty/Metis
./get.Metis
cd $MYHOME/Ipopt/ThirdParty/Mumps
./get.Mumps
cd $MYHOME/Ipopt
./configure --enable-static coin_skip_warn_cxxflags=yes
make -j4
make install
echo -e '\n\n\e[33m \e[1mHit enter if compilation of Ipopt, Metis, and Mumps has been succesfull. Otherwise press Ctrl-C and resolve arisen errors.\e[0m'
read x


#### Download and install ADOLC and ColPack ####
### First ColPack ### 
cd $MYHOME
wget --continue www.coin-or.org/download/source/ADOL-C/$ADOL_FILE
mkdir Adol-c && tar zxvf $ADOL_FILE -C Adol-c --strip-components 1
mv $ADOL_FILE $MYHOME/Downloads
cd $MYHOME/Adol-c

sed -i 's|export HOME=`pwd`||g' configure
sed '2iexport HOME=`pwd`' configure > tmp.txt
mv tmp.txt configure

mkdir ./ThirdParty
cd ./ThirdParty
wget --continue http://archive.ubuntu.com/ubuntu/pool/universe/c/colpack/$COLPACK_FILE
mkdir ColPack && tar zxvf $COLPACK_FILE -C ColPack --strip-components 1
cd ColPack
./autoconf.sh
make
sudo make install
echo -e '\n\n\e[33m \e[1mHit enter if compilation of Colpack has been succesfull. Otherwise press Ctrl-C and resolve arisen errors.\e[0m'
read x


### Second Adol ###
sudo cp -P ./build/lib/libCol* /usr/lib
cd $MYHOME/Adol-c
chmod +x configure
./configure --enable-sparse --with-colpack=$MYHOME/Adol-c/ThirdParty/ColPack/build
make
make install
sudo cp -P $MYHOME/Adol-c/adolc_base/lib64/lib* /usr/lib
sudo cp -r $MYHOME/Adol-c/adolc_base/include/* /usr/include/
echo -e '\n\n\e[33m \e[1mHit enter if compilation of Adol has been succesfull. Otherwise press Ctrl-C and resolve arisen errors.\e[0m'
read x


#### Download and install PDFlib ####
cd $MYHOME
wget --continue https://fossies.org/linux/misc/old/$PDFLIB_FILE
mkdir PDFLib && tar zxvf $PDFLIB_FILE -C PDFLib --strip-components 1
mv $PDFLIB_FILE $MYHOME/Downloads
cd PDFLib
./configure
make; sudo make install
sudo ldconfig
echo -e '\n\n\e[33m \e[1mHit enter if compilation of PDFLib has been succesfull. Otherwise press Ctrl-C and resolve arisen errors.\e[0m'
read x


#### Download and install GNUplot ####
# cd $MYHOME/Downloads
# wget --continue https://sourceforge.net/projects/gnuplot/files/gnuplot/4.4.0/gnuplot-4.4.0.tar.gz/download
# mv download gnuplot-4.4.0.tar.gz
# tar zxvf gnuplot-4.4.0.tar.gz
# sudo apt-get -y install libx11-dev libxt-dev libgd2-xpm-dev libreadline6-dev
sudo apt-get -y install libx11-dev libxt-dev libreadline6-dev libgd-dev gnuplot-qt
# cd gnuplot-4.4.0
# ./configure -with-readline=gnu -without-tutorial
# make;sudo make install


#### Download and extract PSOPT #### 
# cd $MYHOME
# wget --continue https://github.com/PSOPT/psopt/archive/master.zip
# unzip master.zip
# mv master.zip $MYHOME/Downloads
cd $MYHOME


#### Download and extract SuiteSparse #### 
wget --continue http://faculty.cse.tamu.edu/davis/SuiteSparse/$SUITESPARSE_FILE
mkdir SuiteSparse && tar zxvf $SUITESPARSE_FILE -C SuiteSparse --strip-components 1
mv $SUITESPARSE_FILE $MYHOME/Downloads
cd $MYHOME


#### Download and extract LUSOL #### 
wget --continue http://www.stanford.edu/group/SOL/software/lusol/$LUSOL_FILE
unzip $LUSOL_FILE
mv $LUSOL_FILE $MYHOME/Downloads
cd $MYHOME


####  Compile SuiteSparse, LUSOL, dmatrix and PSOPT #### 
make all
echo -e '\n\n\e[33m \e[1mPSOPT installation script completed\e[0m'
echo        '\e[33m \e[1mPSOPT installed in \e[0m' $MYHOME 

