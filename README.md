# sgd_grid - generate grids of 1D self-gravitating disc models
================================================================

This repository computes physical properties of self-gravitating protostellar discs.
The code returns a grid of disc models as a function of input accretion rate and outer disc radius (for a given set of input disc parameters).

This code has been used in the following publications:

*Forgan & Rice (2011), MNRAS, 417, pp 1928-1937*

*Forgan & Rice (2013a), MNRAS, 430, pp 2082-2089*

*Forgan & Rice (2013b), MNRAS, 432, pp 1796-1801*

*Forgan et al (2016), MNRAS, 463, pp 957-964*

with modified versions also appearing in

*Hall, Forgan et al (2016), MNRAS, 458, pp 306-318*

If you plan to use this code for your own research, or if you would like to contribute to this repository then please get in touch with a brief description of what you would like to do.  I will seek co-authorship for any subsequent publications.


Features
--------

This code produces three Fortran programs:

`sgd_grid`
`calc_observables`
`calc_pebble_accretion`

`sgd_grid` calculates a grid of steady state self-gravitating disc models,
given a range of values for accretion rate and disc outer radii. 

The models are computed assuming a fixed Toomre Q parameter, and also provide fragmentation criteria using the Jeans mass in 
spiral arm formalism of (full details in Forgan and Rice (2011), MNRAS, 417, pp 1928-1937).  This also allows an estimate of the initial disc fragment mass.

`calc_observables` is run on the output from `sgd_grid`
to produce observed fluxes at different wavelengths, as well as estimates of the observed disc mass from dust emission and the true disc mass (see Forgan & Rice 2013b, Forgan et al 2016)

`calc_pebble_accretion` is also run on output from `sgd_grid` to determine whether streaming instabilities are likely for grains of a given size/Stokes number,
as well as the expected pebble accretion of disc fragments


Compiling and Running
---------------------

The code is written in FORTRAN 90 throughout. The supplied Makefile compiles the main 
program using gfortran through the command

`> make`

And run using the command

`> ./sgd_grid`

The `calc_observables` and `calc_pebble_accretion` codes are compiled using the command
`> make calc_observables`
`> make calc_pebble_accretion`

And run with

`> ./calc_observables`
`> ./calc_pebble_accretion`

Alternatively, all three can be compiled with

`> make all`

The input parameters are specified in `sgd_grid.params`, `calc_observables.params' and `calc_pebble_accretion.params` respectively

Plotting
--------

The output files can be plotted using Python scripts found in the
`plot/` directory

The scripts were developed in Python 2.7, and depend on numpy and matplotlib

License
-------

This code is licensed under the MIT license - see LICENSE.md
