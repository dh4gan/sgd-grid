The main program calculates a grid of steady state self-gravitating disc models,
given a range of values for accretion rate and disc outer radii.

The program also calculates fragmentation criteria using the Jeans mass in 
spiral arm formalism of Forgan and Rice (2011), MNRAS, 417, pp 1928-1937.

A second program, calc_observables, can also be run on the output from the main program
to produce observed fluxes at different wavelengths, as well as estimates of the observed
disc mass from dust emission, and the true disc mass.

The code is written in FORTRAN 90 throughout. The supplied Makefile compiles the main 
program using gfortran through the command

`> make`

And run using the command

`> ./sgd_grid`

The calc_observables code is compiled using the command
`> make calc_observables`

And run with

`> ./calc_observables`


The input parameters are specified in `sgd_grid.params`
and `calc_observables.params' respectively

The output files can be plotted using Python scripts found in the
`plot_sgd_grid` repository

