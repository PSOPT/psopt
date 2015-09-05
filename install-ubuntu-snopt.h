sudo apt-get install autoconf libtool
cd snopt-interface
sh autogen.sh
./configure
make
make install
cd ..
sh use_ipopt_and_snopt.sh
