#####################################################
###                                               ###
###  Makefile for disc survey routine   	  ###
###  (updated throughout to Fortran 90)           ###
###                                               ###
###         Duncan H. Forgan 18/02/2011           ###
###       				          ###
###                                               ###
#####################################################

# Compiler variables for WKMR/DHF:
FC     = gfortran

# For big endian files generated from stacpolly use these flags
#FFLAGS = -O3 -fPIC -frecord-marker=4 -fconvert=swap -fdefault-real-8

# For files generated from stacpolly use these flags
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8 -Wunused 

# For files generated on seaforth use these flags
#FFLAGS = -O3 -frecord-marker=4  -Wall

# Create object files:
%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = eosmodule.f90 main.f90 eosread.f90 \
	eos_cs.f90 get_kappa.f90
OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: selfgravdisc_modelgrid

selfgravdisc_modelgrid:  $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)

calc_observables: calc_observables.f90
	 $(FC) $(FFLAGS) calc_observables.f90 -o calc_observables
	rm -f calc_observables.o

# Clean statements:
clean: 
	\rm *.o *.mod selfgravdisc_modelgrid calc_observables

# End Makefile
