PROGRAM calc_pebble_accretion
  ! Program calculates observables for a grid of pre-generated
  ! self-gravitating disc models
  ! Assumes a dust temperature and opacity properties to do so

  use eosdata, only: umass, udist
  implicit none

  integer :: irad,imdot,irrchoice
  integer :: nrad, nmdot

  real :: mdotmax,mdotmin, mdotvisc, metallicity, mstar
  real :: qratio, T, gamma_sigma,gamma_omega
  real :: Q_irr, T_irr,dr, rmax,rmin, betac

  real, parameter :: G = 6.67e-8
  real, parameter :: pi = 3.14159265285
  real, parameter :: twopi = 2.0*pi
  real, parameter :: sigma_SB = 5.67e-5
  real, parameter :: c = 3.0e10
  real, parameter :: tolerance = 1.0e-5
  real, parameter :: Qfrag = 2.0
  real, parameter :: mjup = 1.8986e30
  real, parameter :: mearth = 5.972e27
  real, parameter :: pc = 206265.0*1.496e13  ! Parsec in cm
  real, parameter :: gamma1 = 4.0

  integer, allocatable,dimension(:) :: stream_unstable
  real, allocatable, dimension(:) :: r, sigma,cs,omega,alpha
  real, allocatable, dimension(:) :: H, eta, vrpeb, rhogas

  integer :: ipebrad,npebrad, irmin_unstable,irmax_unstable, jrad, imax_peb
  real :: tstop, zpeb, beta_peb, rmax_peb, mplanet, Hp_to_Hg, dlogrhodr
  
  real :: rpeb, rdotpeb, tpeb, mdotpebble, rhill, rhop_rhog
  real :: rmin_unstable, rmax_unstable, width_unstable, h_unstable
  real :: planet_pebaccrete, mcross, eff_pebble, sigma_p
  real :: rpeb_accretemax, tpeb_accretemax
  real :: planet_accretemax, mcross_accretemax, eff_accretemax

  logical :: inner_radius
  logical :: outer_radius

  character(3) :: imdot_char
  character(14) :: mdot_char
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


  outputlog = trim(outputprefix)//".log"
  print*, "Log data for all gas accretion rates to be outputted to ", trim(outputlog)
  print*, "Data for each accretion rate to be outputted to files of form ", trim(outputlog)//".mdot.<num>"

  !******************************************
  ! 1. Read in disc and recompute essential gas  properties
  !******************************************

  OPEN(10,file=inputfile ,status='unknown')

  ! Read input file header

  read(10,*) nrad, nmdot, rmin,rmax,mdotmin,mdotmax,Mstar, metallicity,&
       gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

  mstar = mstar*umass
  mplanet = mplanet*mearth

  print*, 'There are ', nrad, ' radii in the input file'

  allocate(r(nrad))
  allocate(sigma(nrad))
  allocate(cs(nrad))
  allocate(H(nrad))
  allocate(rhogas(nrad))
  allocate(omega(nrad))
  allocate(alpha(nrad))
  allocate(eta(nrad))
  allocate(vrpeb(nrad))
  allocate(stream_unstable(nrad))

  ! Find index of maximum radius for pebble accretion (npebrad)
  dr = (rmax-rmin)/REAL(nrad)
  !dr = dr*udist
  npebrad = int(rmax_peb/dr) + 1

  ! Write header for log output file
  OPEN(20,file=outputlog,status='unknown')
  write(20,*) nrad,nmdot,rmin,rmax,mdotmin,mdotmax,Mstar,metallicity,&
       gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr, tstop,zpeb,beta_peb

  
  ! Loop over accretion rates
  DO imdot = 1, nmdot

     r(:) = 0.0
     sigma(:) = 0.0
     cs(:) = 0.0
     H(:) = 0.0
     rhogas(:) = 0.0
     omega(:) = 0.0
     alpha(:) = 0.0
     eta(:) = 0.0
     vrpeb(:) = 0.0
     stream_unstable(:) = 0

     ! Read entire disc at this accretion rate
     DO irad = 1,nrad

        read(10,*) r(irad), mdotvisc, qratio, sigma(irad), cs(irad),omega(irad),&
             T,betac,alpha(irad)

        !print*, r(irad), mdotvisc, qratio, sigma(irad), cs(irad),omega(irad),&
        !     T,betac,alpha(irad)


        ! Be careful with units here
        H(irad) = (cs(irad)/omega(irad)) ! in units of cm
        rhogas(irad) = sigma(irad)/(2.0*H(irad))
        r(irad) = r(irad)*udist
                
        ! Calculate pressure gradient and eta - sub Keplerian parameter

        if(irad>1) then
           dlogrhodr = (log(rhogas(irad)*cs(irad)*cs(irad)) - log(rhogas(irad-1)*cs(irad-1)*cs(irad-1)))/&
                (log(r(irad) - log(r(irad-1))))
        else
           dlogrhodr = 0.0
        endif

        eta(irad) = -0.5*(H(irad)*H(irad))/(r(irad)*r(irad)) * dlogrhodr

        ! Radial velocity of pebbles
        vrpeb(irad) = 2.0*eta(irad)*omega(irad)*r(irad)*tstop/(1.0 + tstop*tstop)

     enddo


     !******************************************************************
     ! 2. Now compute pebble accretion properties as a function of rpeb
     !******************************************************************

     ! Specify outputfile for this accretion rate

     write(imdot_char, '(I3.3)') imdot
     output_mdotfile = trim(outputprefix)//".mdot."//trim(imdot_char)

     open(30, file=output_mdotfile, status='unknown')
     write(30,*) mdotvisc

     ! Variables store maximum value at this accretion rate
     planet_accretemax = 0.0
     rpeb_accretemax = 0.0
     mcross_accretemax = 0.0
     eff_accretemax = 0.0

     ! convert mdot into cgs
     mdotvisc = mdotvisc*umass/3.15e7

     ! Now loop again over radius, where radius now refers to rpeb

     print*, 'Number of pebble radii' , npebrad
     do ipebrad =1,npebrad

        rpeb = r(ipebrad)

        !**************************************
        ! 2a Calculate growth rate of pebbles
        !**************************************

        ! Pebble growth timescale at this radius
        tpeb = beta_peb/(zpeb*omega(ipebrad))

        ! Pebble front growth rate
        rdotpeb = 0.6666*(G*Mstar*zpeb*zpeb/(tpeb*beta_peb*beta_peb))**(0.333)

        ! Compute mdotpebble (g s^-1)
        mdotpebble = 2.0*pi*rpeb*rdotpeb*zpeb*sigma(ipebrad)

!        print*, mdotvisc*3.15e7/umass, mstar/umass, rpeb/udist, tpeb/3.15e7, rdotpeb*3.15e7/udist, mdotpebble*3.15e7/umass
        !************************************************
        ! 2b Find regions where streaming instability active
        !************************************************

        stream_unstable(:) = 0
        inner_radius = .false.
        outer_radius = .false.
        rmin_unstable = 0.0
        rmax_unstable = 0.0
        irmin_unstable = 0
        irmax_unstable = 0

        do jrad =2,ipebrad

           Hp_to_Hg = sqrt(tstop/alpha(jrad))
           ! Compute rhop/rhog interior to pebble radius
           rhop_rhog = mdotpebble/(Hp_to_Hg*sigma(jrad)*2.0*pi*r(jrad)*vrpeb(jrad))

           !print*, jrad, rhop_rhog, Hp_to_Hg, sigma(jrad), vrpeb(jrad)
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

        !print*, 'out', irmin_unstable, irmax_unstable
        ! If streaming zone goes right to pebble front, make sure outer radius is found
       if (inner_radius.eqv. .true. .and. outer_radius.eqv..false.) then           
          rmax_unstable = r(ipebrad)
           irmax_unstable = ipebrad
           outer_radius = .true.
        endif
        
        !****************************************************
        ! 3 Calculate expected core accretion rates
        !****************************************************

        ! Find width of unstable region 

        width_unstable = rmax_unstable - rmin_unstable
       
        ! Only do the calculations for non-zero widths
        if(width_unstable > 0.0) then

           h_unstable = H(irmin_unstable)/rmin_unstable
        ! Compute surface density of pebbles in here
        sigma_p = mdotpebble/(2.0*pi*r(irmin_unstable)*vrpeb(irmin_unstable))

       ! print*, rmin_unstable/udist, rmax_unstable/udist, sigma(irmin_unstable), sigma_p
        ! Compute Pebble Accretion Rate for M= 1 Mjup
        ! Mdot = 2 R_H^2 omega tstop *sigma p

        rhill = rmin_unstable*(mplanet/(3.0*mstar))**0.333

        !print*, rmin_unstable/udist, rhill/udist, mplanet, mstar

        planet_pebaccrete = 2.0* rhill*rhill*omega(irmin_unstable)*(tstop**0.66666)*sigma_p

        if(planet_pebaccrete>mdotpebble) planet_pebaccrete = mdotpebble

        ! compute pebble accretion efficiency
        eff_pebble = planet_pebaccrete/mdotpebble

        ! Compute crossing mass at this rpeb
        mcross = sqrt(3.0*pi*width_unstable*alpha(irmin_unstable)*mdotpebble*eff_pebble/&
             (rmin_unstable*mdotvisc*gamma1))*h_unstable * h_unstable* Mstar

        write(*,'(10(1P,e8.1,1X))'), mcross/mjup, eta(irmin_unstable)/(h_unstable*h_unstable), &
             width_unstable/udist, H(irmin_unstable)/r(irmin_unstable),mdotpebble/mdotvisc, &
             alpha(irmin_unstable), mdotpebble, eff_pebble, mdotvisc, Mstar/umass
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
        
        write(30,*) mdotvisc*3.15e7/umass, rpeb/udist,tpeb/3.15e7, rdotpeb*3.15e7/udist, &
             mdotpebble*3.15e7/mjup, rmin_unstable/udist, rmax_unstable/udist, &
             mcross/mjup, planet_pebaccrete*3.15e7/mjup, eff_pebble
        endif
     enddo
     close(30)

     ! Now write log file for this mdot
     ! mdot: maximum mdotp, rg(max mdotp), Mcross(max mdotp), efficiency
     write(20,*) mdotvisc*3.15e7/umass, rpeb_accretemax/udist, tpeb_accretemax/3.15e7, &
          planet_accretemax*3.15e7/mjup, &
          mcross_accretemax/mjup, eff_accretemax

  enddo
        
close(20)

!END of program

end PROGRAM calc_pebble_accretion
