#!/bin/sh
set -ex
wget https://github.com/deepzot/likely/archive/master.tar.gz 
tar -xzvf master.tar.gz
cd likely-master/build && ../configure --prefix=/usr && make && sudo make install