turbo-octo-spice 
================

[![Build Status](https://travis-ci.org/dmargala/turbo-octo-spice.svg?branch=master)](https://travis-ci.org/dmargala/turbo-octo-spice)

Two-point correlation function estimator for SDSS LyA observations.

 * [github repo](https://github.com/dmargala/turbo-octo-spice)
 * [documentation](http://dmargala.github.io/turbo-octo-spice)

Dependencies
------------

 * boost
 * likely
  * gsl
  * Minuit2
 * cosmo
  * fftw
 * libhdf5 (`sudo apt-get install libhdf5-dev`)
 * cfitso 
 * ccfits (`sudo apt-get install libccfits-dev`)
 * healpix


Instructions for installing `likely`, `cosmo`, and their dependencies are available [here](http://darkmatter.ps.uci.edu/wiki/DeepZot). Some tips for installing boost on OSX can be found [here](http://darkmatter.ps.uci.edu/wiki/BoostMavericks). The `.travis.yml` file and various `install-*.sh` shell scripts might also be helpful.

Motivation
----------

This project started out as an attempt at efficient mock generation and analysis of Lyman-alpha forest surveys. Mostly still in the analysis phase of things...

The BOSS quasar sample contains ~100k quasar spectra. Individual spectra sample the Lyman-alpha forest along the line of sight to the quasar. There are typically about ~100 pixels in that contain information about the Lyman-alpha forest density field in each spectrum.

The plot below shows a comparison for generating and analyzing a gaussian random field on a cubic grid. Each line color corresponds to a different system the benchmark was performed on. Blue corresponds to my laptop, red corresponds to UCI HPC public GPU server, green corresponds to our UCI group server 'darkmatter', and gray corresponds to Table 1 of "Cosmological calculations on the GPU" (Bard et al., 2012). The solid lines indicate a "brute force" GPU implementation xi estimator while the dashed lines indicate a "brute force" CPU implementation xi estimator. The dot-dashed line indicates the gaussian random field generation performance of 'cosmogrf'.

![](images/xiperfplot.png "Xi Performance Comparison")
