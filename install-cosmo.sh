#!/bin/sh
set -ex
wget https://github.com/deepzot/cosmo/archive/master.tar.gz 
tar -xzvf master.tar.gz
cd cosmo-master/build && ../configure --prefix=/usr && make && sudo make install