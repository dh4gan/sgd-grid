PROGRAM calc_observables
  ! Program calculates observables for a grid of pre-generated
  ! self-gravitating disc models
  ! Assumes a dust temperature and opacity properties to do so

  use eosdata, only: umass, udist, Boltzmannk
  implicit none

  integer :: irad,imdot,irrchoice
  integer :: nrad, nmdot

  real, parameter :: pi = 3.14159265285
  real, parameter :: twopi = 2.0*pi
  real, parameter :: sigma_SB = 5.67e-5
  real, parameter :: c = 3.0e10
  real, parameter :: tolerance = 1.0e-5
  real, parameter :: Qfrag = 2.0
  real, parameter :: mjup = 1.8986e30
  real, parameter :: pc = 206265.0*1.496e13  ! Parsec in cm


  real :: beta_k,distance, lambda0,kappa0,lambda,Tdust
  real :: Mstar, metallicity, gamma_sigma,gamma_omega
  real :: nu0, nu, kappa_nu,flux_nu,fluxtot,mflux
  real :: rmin, rmax, mdotmin,mdotmax, mdotvisc
  real :: Q_irr, T_irr, r, qratio,mtot
  real :: sigma, omega, T, cs,betac,alpha,dr

  character(100) ::inputfile, outputfile


  ! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"  
  print*, " OBSERVABLES FROM SELF-GRAVITATING DISC MODELS " 
  print*, "     Created by D. Forgan, Jan 2013            "
  print*, " (Relies on input from selfgravdisc_modelgrid) "
  print*, "-----------------------------------------------"
  print*, " "
  print*, " "
  print*, " input files: ./calc_observables.params ./myeos.dat" 
  print*, " "
  print*, "-----------------------------------------------"
  print*, " "

  OPEN(10,file='calc_observables.params', status='unknown')
  read(10,*) inputfile
  READ(10,*) outputfile
  READ(10,*) distance
  READ(10,*) lambda
  READ(10,*) kappa0
  read(10,*) lambda0
  read(10,*) beta_k
  read(10,*) Tdust
  close(10)


 

  print*, 'Reading file ', inputfile
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
        print*, flux_nu, Boltzmannk, c, nu, T,sigma,kappa_nu, distance, twopi*r*dr

        !	Calculate total flux emitted up to and including this radius

        mtot = qratio*Mstar
        fluxtot = fluxtot + flux_nu

        ! Calculate predicted mass from two different fluxes (given estimated dust temperature)

        mflux = distance*distance*fluxtot*c*c/(2.0*kappa_nu*nu*nu*Boltzmannk*Tdust)

        !	Write to file

        WRITE(20,*) r/udist, mdotvisc, qratio,mtot, flux_nu/1e-23, fluxtot/1e-23,mflux/umass

     ENDDO
  ENDDO


END PROGRAM calc_observables
