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
FFLAGS = -O3 -frecord-marker=4 -fdefault-real-8  

# For files generated on seaforth use these flags
#FFLAGS = -O3 -frecord-marker=4  -Wall

# Create object files:
%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

SOURCESAF90 = eosmodule.f90 main.f90 eos_T.f90 eosread.f90 \
	eos_cs.f90 get_kappa.f90
OBJECTSA    = $(SOURCESAF90:.f90=.o)

# Create executable files:
build: survey_jeans

survey_jeans:  $(OBJECTSA)
	$(FC) $(FFLAGS) -o $@ $(OBJECTSA)

# Clean statements:
clean: 
	\rm *.o *.mod survey_jeans

# End Makefile
