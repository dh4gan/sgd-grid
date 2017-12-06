PROGRAM calc_pebble_accretion
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
  print*, " PEBBLE ACCRETION IN SELF-GRAVITATING DISC MODELS " 
  print*, "     Created by D. Forgan, Dec 2017            "
  print*, " (Relies on input from sgd_grid) "
  print*, "-----------------------------------------------"
  print*, " "
  print*, " "
  print*, " input files: ./calc_observables.params" 
  print*, " "
  print*, "-----------------------------------------------"
  print*, " "

  OPEN(10,file='calc_pebble_accretion.params', status='unknown')
  read(10,*) inputfile
  READ(10,*) outputfile
  READ(10,*) tstop ! Stopping time (dimensionless)
  read(10,*) rmax_peb ! Maximum radius to consider pebble accretion

  !TODO - finish read in of parameters
  close(10)

  print*, 'Reading file ', inputfile
  print*, 'Outputting to file ', outputfile

  OPEN(10,file=inputfile ,status='unknown')

  ! Read input file header

  read(10,*) nrad, nmdot, rmin,rmax,mdotmin,mdotmax,Mstar, metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr


  ! TODO - find index of maximum radius for pebble accretion (npebrad)
  dr = (rmax-rmin)/REAL(nrad)
  dr = dr*udist
  distance = distance*pc
  nu = c/lambda
  nu0 = c/lambda0
  kappa_nu = kappa0*(nu/nu0)**(beta_k)

  OPEN(20,file=outputfile,status='unknown')
  write(20,*) nrad,nmdot,nu,lambda*1e4, rmin, rmax,mdotmin,mdotmax,Mstar,metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr
  ! Loop over accretion rates

  DO imdot = 1, nmdot

     mtot = 0.0
     fluxtot = 0.0

     ! Read entire disc over pebble growth radius
     DO irad = 1,npebrad

        read(10,*) r, mdotvisc, qratio, sigma, cs,omega,T,betac,alpha

        r = r*udist
        
     ! Calculate pressure gradient and eta - sub Keplerian parameter

     enddo

     ! Now loop over radius, where radius now refers to rpeb

     do irad =1,npebrad

        rpeb = 
        ! Compute mdotpebble


        do jrad =1,irad
        ! Compute rhop/rhog interior to pebble radius
           r=

        ! If rhop/rhog >1, mark this radius as streaming unstable region
           
        enddo

        ! Compute crossing mass at this rpeb

        ! Compute Pebble Accretion Rate for Mjup
        ! compute pebble accretion efficiency



        ! Write to file for this rpeb
        ! mdot rpeb --> mdotpebble, r1, r2, Mcross, Mpl
        
     enddo

     ! Now write file for this mdot
     ! mdot: maximum mdotp, rg(max mdotp), Mcross(max mdotp), efficiency
  enddo
        

!END of program

end PROGRAM calc_pebble_accretion
