PROGRAM calc_observables
  
  ! Program reads data from sgd_grid
  ! calculates astronomical observables for the model grid
  ! Assumes a dust temperature and opacity properties to do so

  use eosdata
  implicit none

  integer :: irad,imdot,irrchoice
  integer :: nrad, nmdot

  real :: beta_k,distance, lambda0,kappa0,lambda,Tdust
  real :: Mstar, metallicity, gamma_sigma,gamma_omega
  real :: nu0, nu, kappa_nu,flux_nu,fluxtot,mflux
  real :: rmin, rmax, mdotmin,mdotmax, mdotvisc
  real :: Q_irr, T_irr, r, qratio,mtot
  real :: sigma, omega, T, cs,betac,alpha,dr

  character(100) ::inputfile,inputprefix, outputfile

  !*******************************************
  ! 1. Get input parameters and set up header:
  !*******************************************

  print*, " "
  print*, "-----------------------------------------------"  
  print*, " OBSERVABLES FROM SELF-GRAVITATING DISC MODELS " 
  print*, "     Created by D. Forgan, Jan 2013            "
  print*, "     Current Version: Jan 2018"
  print*, " (Relies on input from sgd_grid) "
  print*, "-----------------------------------------------"
  print*, " "
  print*, " "
  print*, " input files: ./calc_observables.params" 
  print*, " "
  print*, "-----------------------------------------------"
  print*, " "

  OPEN(10,file='calc_observables.params', status='unknown')
  read(10,*) inputprefix
  READ(10,*) outputfile
  READ(10,*) distance
  READ(10,*) lambda  ! Input in nanometres
  READ(10,*) kappa0
  read(10,*) lambda0 ! input in nanometres
  read(10,*) beta_k 
  read(10,*) Tdust
  close(10)

  ! Convert wavelengths into cm
  lambda = lambda*1e-4
  lambda0 = lambda0*1e-4
 
  inputfile = trim(inputprefix)//'.sgdmodel'
  print*, 'Reading file ', inputfile
  outputfile = trim(inputprefix)//'.observe'
  print*, 'Outputting to file ', outputfile

  OPEN(10,file=inputfile ,status='unknown')
  read(10,*) nrad, nmdot, rmin,rmax,mdotmin,mdotmax,Mstar, metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

  dr = (rmax-rmin)/REAL(nrad)
  dr = dr*udist
  distance = distance*pc
  nu = c/lambda
  nu0 = c/lambda0
  kappa_nu = kappa0*(nu/nu0)**(beta_k)

  OPEN(20,file=outputfile,status='unknown')
  write(20,*) nrad,nmdot,nu,lambda*1e4, rmin, rmax,mdotmin,mdotmax,Mstar,metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

  !**************************
  ! 2. Compute observables
  !**************************

  ! Loop over parameters

  DO imdot = 1, nmdot

     mtot = 0.0
     fluxtot = 0.0

     ! Loop over radii
     DO irad = 1,nrad

        read(10,*) r, mdotvisc, qratio, sigma, cs,omega,T,betac,alpha

        r = r*udist
        
        ! Calculate flux generated in this annulus

        if (sigma*kappa_nu < 1.0) THEN
           flux_nu = (2.0*Boltzmannk/(c*c))*nu*nu*kappa_nu*sigma*T*twopi*r*dr/(distance*distance)
        else
           flux_nu = (2.0*Boltzmannk/(c*c))*nu*nu*(T/(sigma*kappa_nu)**0.25)*twopi*r*dr/(distance*distance)

        ENDIF

        ! Calculate total flux emitted up to and including this radius

        mtot = qratio*Mstar
        fluxtot = fluxtot + flux_nu

        ! Calculate predicted mass from the observed flux
        ! (given estimated dust temperature and assuming optically thin)

        mflux = distance*distance*fluxtot*c*c/(2.0*kappa_nu*nu*nu*Boltzmannk*Tdust)

        !********************
        ! 3.	Write to file 
        !********************

        ! flux output in mJy

        WRITE(20,*) r/udist, mdotvisc, qratio,mtot, flux_nu/1e-26, fluxtot/1e-26,mflux/umass, mflux/(umass*mtot)

     ENDDO
  ENDDO

close(10)
close(20)

END PROGRAM calc_observables
