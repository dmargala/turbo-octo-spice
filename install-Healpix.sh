#!/bin/sh
set -ex
curl -k -L http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz | tar zx
cd Minuit2-5.34.14 && ./configure --disable-openmp --prefix=/usr && make && sudo make install


# curl -k -L http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz | tar zx
# cd Minuit2-5.34.14
# ./configure --disable-openmp --prefix=/usr/local
# make
# sudo make install