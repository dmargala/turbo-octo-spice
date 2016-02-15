#!/bin/sh
set -ex
curl -k -L https://github.com/deepzot/cosmo/archive/master.tar.gz | tar xz
cd cosmo-master/build && ../configure --prefix=/usr && make && sudo make install