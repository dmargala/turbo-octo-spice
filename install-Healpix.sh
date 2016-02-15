#!/bin/sh
set -ex
wget http://downloads.sourceforge.net/project/healpix/Healpix_3.30/Healpix_3.30_2015Oct08.tar.gz
tar -xzvf Healpix_3.30_2015Oct08.tar.gz
cd Healpix_3.30/src/cxx && autoconf &&./configure --prefix=/usr && make && sudo make install
