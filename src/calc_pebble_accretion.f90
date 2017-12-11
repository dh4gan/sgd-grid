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

  integer, allocatable,dimension(:) :: stream_unstable
  real, allocatable, dimension(:) :: r, sigma,cs,omega,alpha
  real, allocatable, dimension(:) :: H, eta, vrpeb, rhogas

  integer :: irad, ipebrad,irmin_unstable,ir_max_unstable, jrad
  real :: tstop, zpeb, beta_peb, rmax_peb, mplanet
  


  character(100) ::inputfile, outputprefix, outputlog, output_mdotfile

  ! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"  
  print*, " PEBBLE ACCRETION IN SELF-GRAVITATING DISC MODELS " 
  print*, "     Created by D. Forgan, Dec 2017            "
  print*, " (Relies on input from sgd_grid) "
  print*, "-----------------------------------------------"
  print*, " "
  print*, " "
  print*, " input files: ./calc_pebble_accretion.params" 
  print*, " "
  print*, "-----------------------------------------------"
  print*, " "

  OPEN(10,file='calc_pebble_accretion.params', status='unknown')
  read(10,*) inputfile
  READ(10,*) outputprefix ! Prefix for output files 
  READ(10,*) tstop ! Stopping time (dimensionless)
  read(10,*) zpeb ! Metallicity of disc
  read(10,*) beta_peb ! growth rate of pebbles t_peb = beta_peb *(zpeb*omega(irad)^-1
  read(10,*) rmax_peb ! Maximum radius to consider pebble accretion  (AU)
  read(10,*) mplanet ! Typical planet mass for pebble accretion (Jupiter masses)
  close(10)

  print*, 'Reading file ', inputfile
  print*, 'Outputting pebble accretion data to files of prefix ', outputprefix


  outputlog = trim(outputprefix)//"_log.pebble"
  print*, "Log data for all gas accretion rates to be outputted to ", trim(outputlog)
  print*, "Data for each accretion rate to be outputted to files of form", trim(outputlog)//"mdot_+++.pebble"

  !******************************************
  ! 1. Read in disc and recompute essential gas  properties
  !******************************************

  OPEN(10,file=inputfile ,status='unknown')

  ! Read input file header

  read(10,*) nrad, nmdot, rmin,rmax,mdotmin,mdotmax,Mstar, metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

  allocate(r(nrad))
  allocate(sigma(nrad))
  allocate(cs(nrad))
  allocate(H(nrad))
  allocate(rhogas(nrad))
  allocate(omega(nrad))
  allocate(eta(nrad))
  allocate(vrpeb(nrad))
  allocate(stream_unstable(nrad))



  ! Find index of maximum radius for pebble accretion (npebrad)
  dr = (rmax-rmin)/REAL(nrad)
  !dr = dr*udist
  imax_peb = int(rmax_peb/dr) + 1

  ! Write header for log output file
  OPEN(20,file=outputlog,status='unknown')
  write(20,*) nrad,nmdot,rmin,rmax,mdotmin,mdotmax,Mstar,metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr, tstop,zpeb,beta_peb

  
  ! Loop over accretion rates
  DO imdot = 1, nmdot


     r(:) = 0.0
     sigma(:) = 0.0
     cs(:) = 0.0
     H(:) = 0.0
     rhogas(:) = 0.0
     omega(:) = 0.0
     eta(:) = 0.0
     vrpeb(:) = 0.0
     stream_unstable(:) = 0

     ! Read entire disc at this accretion rate
     DO irad = 1,nrad

        read(10,*) r(irad), mdotvisc, qratio, sigma(irad), cs(irad),omega(irad),&
             T,betac,alpha(irad)

        ! TODO - check units
        H(irad) = cs(irad)/omega(irad)
        rhogas(irad) = sigma(irad)/(2.0*H(irad))
        !r(irad) = r(irad)*udist
        
        ! Calculate pressure gradient and eta - sub Keplerian parameter

        if(irad>1) then
           dlogrhodr = (log(rhogas(irad)*cs(irad)*cs(irad)) - log(rhogas(irad)*cs(irad)*cs(irad)))/(log(r(irad) - log(r(irad-1))))
        else
           dlogrhodr = 0.0
        endif

        eta(irad) = -0.5(H(irad)*H(irad))/(r(irad)*r(irad)) * dlogrhodr

        ! Radial velocity of pebbles
        vrpeb(irad) = -2.0*eta(irad)*omega(irad)*r(irad)*tstop/(1.0 + tstop*tstop)

     enddo


     !******************************************************************
     ! 2. Now compute pebble accretion properties as a function of rpeb
     !******************************************************************

    

     ! Specify outputfile for this accretion rate

     write(mdot_char,'(e5.2)') mdotvisc
     output_mdotfile = trim(outputprefix)//"mdot_"//trim(mdot_char)//".pebble"

     open(30, file=output_mdotfile, status='unknown')

     ! Variables store maximum value at this accretion rate
     planet_accretemax = 0.0
     rpeb_accretemax = 0.0
     mcross_accretemax = 0.0
     eff_accretemax = 0.0

     mdotvisc = mdotvisc*umass/3.15e7

     ! Now loop again over radius, where radius now refers to rpeb

     do ipebrad =1,npebrad

        rpeb = r(ipebrad)

        !**************************************
        ! 2a Calculate growth rate of pebbles
        !**************************************

        ! Pebble growth timescale at this radius
        tpeb = beta_peb/(zpeb*omega(ipebrad))

        ! Pebble front growth rate
        rdotpeb = 0.6666*(G*Mstar*zpeb*zpeb/(tpeb*beta_peb*beta_peb))**(0.333)

        ! Compute mdotpebble
        mdotpebble = 2.0*pi*rpeb*rdotpeb*zpeb*sigma(ipebrad)



        !************************************************
        ! 2b Find regions where streaming instability active
        !************************************************

        stream_unstable(:) = 0
        inner_radius = .false.
        outer_radius = .false.

        do jrad =1,ipebrad

           Hp_to_Hg = sqrt(tstop/alpha(jrad))
           ! Compute rhop/rhog interior to pebble radius
           rhop_rhog = mdotpebble/(Hp_to_Hg*sigma(jrad)*2.0*pi*r(jrad)*vrpeb(jrad))

           ! If rhop/rhog >1, mark this radius as streaming unstable region
           
           if(rhop_rhog>1.0) then
              stream_unstable(jrad) = 1

              ! Record minimum and maximum values of streaming regions (assuming single region only
              if(inner_radius.eqv..false.) then
                 rmin_unstable = r(jrad)
                 irmin_unstable = jrad
                 inner_radius = .true.
              endif
           else
              if(stream_unstable(jrad-1)==1) then
                 rmax_unstable = r(jrad)
                 irmax_unstable = jrad
                 outer_radius = .true.
              endif
           endif

        enddo

        !****************************************************
        ! 3 Calculate expected core accretion rates
        !****************************************************

        ! Find width of unstable region 

        width_unstable = rmax_unstable - rmin_unstable

        ! Compute surface density of pebbles in here
        sigma_p = mdotpebble/(2.0*pi*r(irmin_unstable)*vr(irmin_unstable))

        ! Compute Pebble Accretion Rate for M= 1 Mjup
        ! Mdot = 2 R_H^2 omega tstop *sigma p

        rhill = r(irmin_unstable)*(mplanet/3.0*mstar)**0.333

        planet_pebaccrete = 2.0* rhill*rhill*omega(irmin_unstable)*(tstop**0.66666)*sigma_p

        ! compute pebble accretion efficiency
        eff_pebble = planet_pebaccrete/mdotpebble

        ! Compute crossing mass at this rpeb
        mcross = sqrt(3.0*pi*width_unstable*alpha(irmin_unstable)*mdotpebble*eff_pebble/&
             (mdotvisc*gamma1))*(H(irmin_unstable)/r(irmin_unstable))**2 * Mstar


        ! Check if accretion rate is at a maximum
        if(planet_pebaccrete > planet_accretemax) then
           planet_accretemax = planet_pebaccrete
           rpeb_accretemax = rpeb
           tpeb_accretemax = tpeb
           mcross_accretemax = mcross
           eff_accretemax = eff_pebble
        endif


        !****************************
        ! 4. Write output data to files
        !****************************

        ! Write to file for this rpeb
        ! mdot rpeb --> mdotpebble, r1, r2, Mcross, Mpl
        
        write(30,*) mdotvisc, rpeb,tpeb, rdotpeb,  mdotpebble, rmin_unstable, rmax_unstable, mcross, planet_pebaccrete, eff_pebble

     enddo
     close(30)

     ! Now write log file for this mdot
     ! mdot: maximum mdotp, rg(max mdotp), Mcross(max mdotp), efficiency
     write(20,*) mdotvisc, rpeb_accretemax, tpeb_accretemax, planet_accretemax, mcross_accretemax, eff_accretemax

  enddo
        
close(20)

!END of program

end PROGRAM calc_pebble_accretion
