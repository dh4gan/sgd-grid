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
  real, parameter :: year = 3.15e7
  real, parameter :: pc = 206265.0*1.496e13  ! Parsec in cm
  real, parameter :: gamma1 = 4.0

  integer, allocatable,dimension(:) :: stream_unstable
  real, allocatable, dimension(:) :: r, sigma,cs,omega,alpha,mjeans
  real, allocatable, dimension(:) :: H, eta, vrpeb, rhogas,etadash

  integer :: ipebrad,npebrad, irmin_unstable,irmax_unstable, jrad
  real :: tstop, zpeb, beta_peb, rmax_peb, mplanet, Hp_to_Hg, Hp, dlogrhodr
  
  real :: rpeb, rdotpeb, tpeb, mdotpebble, rhill, rhop_rhog
  real :: rmin_unstable, rmax_unstable, width_unstable, h_unstable
real :: bcross_reduce, bcross, rhill_reduce, zeta, Chi
  real :: planet_pebaccrete, mcross, eff_pebble, effpeb
  real :: rpeb_accretemax, tpeb_accretemax, mdotpebble_accretemax, mjeans_accretemax
  real :: planet_accretemax, mcross_accretemax, eff_accretemax

  logical :: inner_radius
  logical :: outer_radius

  character(3) :: imdot_char
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


  outputlog = trim(inputfile)//".pebble.log"
  print*, "Log data for all gas accretion rates to be outputted to ", trim(outputlog)
  print*, "Data for each accretion rate to be outputted to files of form ", trim(inputfile)//"<num>.pebble"

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
  allocate(mjeans(nrad))
  allocate(eta(nrad))
  allocate(etadash(nrad))
  allocate(vrpeb(nrad))
  allocate(stream_unstable(nrad))

  ! Find index of maximum radius for pebble accretion (npebrad)
  dr = (rmax-rmin)/REAL(nrad)
  !dr = dr*udist
  npebrad = int(rmax_peb/dr) + 1

  if(npebrad > nrad) npebrad = nrad


  ! Convenience functions for calculating pebble accretion (Ida et al 2016)

  zeta = 1/(1+tstop*tstop)
  Chi = sqrt((1.0 + 4.0*tstop*tstop))/(1.0+tstop*tstop)


  ! Write header for log output file
  OPEN(20,file=outputlog,status='unknown')
  write(20,*) nrad,nmdot,rmin,rmax,mdotmin,mdotmax,Mstar,&
       tstop,zpeb,beta_peb,&
       metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

  
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
     mjeans(:) = 0.0
     stream_unstable(:) = 0

     ! Read entire disc at this accretion rate
     DO irad = 1,nrad

        read(10,*) r(irad), mdotvisc, qratio, sigma(irad), cs(irad),omega(irad),&
             T,betac,alpha(irad), mjeans(irad)

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

        ! Compute etadash function (for pebble accretion rates)
        etadash(irad) = Chi*eta(irad)


        ! Radial velocity of pebbles
        vrpeb(irad) = 2.0*eta(irad)*omega(irad)*r(irad)*tstop/(1.0 + tstop*tstop)
!        write(76,*) r(irad)/udist, rhogas(irad)*cs(irad)*cs(irad), &
!             H(irad)/udist, rhogas(irad), dlogrhodr, eta(irad), (r(irad)/vrpeb(irad))/3.15e7
     enddo

     !******************************************************************
     ! 2. Now compute pebble accretion properties as a function of rpeb
     !******************************************************************

     ! Specify outputfile for this accretion rate

     write(imdot_char, '(I3.3)') imdot
     output_mdotfile = trim(inputfile)//"."//trim(imdot_char)//".pebble"

     open(30, file=output_mdotfile, status='unknown')
     write(30,*) nrad,nmdot,rmin,rmax,mdotvisc,Mstar,tstop,zpeb,beta_peb,&
    metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

     ! Variables store maximum values at this gas accretion rate
     planet_accretemax = 0.0
     mjeans_accretemax = 0.0
     mdotpebble_accretemax = 0.0
     rpeb_accretemax = 0.0
     tpeb_accretemax = 0.0
     mcross_accretemax = 0.0
     eff_accretemax = 0.0

     ! convert mdot into cgs
     mdotvisc = mdotvisc*umass/3.15e7

     ! Now loop again over radius, where radius now refers to rpeb

     print*, 'Number of pebble radii' , npebrad
     do ipebrad =1,npebrad

        rpeb = r(ipebrad)

        !**************************************
        ! 2 Calculate growth rate of pebble front
        !**************************************

        ! Pebble growth timescale at this radius
        tpeb = beta_peb/(zpeb*omega(ipebrad))

        ! Pebble front growth rate
        rdotpeb = 0.6666*(G*Mstar*zpeb*zpeb/(tpeb*beta_peb*beta_peb))**(0.333)

        ! Compute mdotpebble (g s^-1)
        mdotpebble = 2.0*pi*rpeb*rdotpeb*zpeb*sigma(ipebrad)

!        print*, mdotvisc*3.15e7/umass, mstar/umass, rpeb/udist, tpeb/3.15e7, rdotpeb*3.15e7/udist, mdotpebble*3.15e7/umass


        !*************************************************************
        ! 3 Compute the pebble accretion rate of any fragment present
        ! (according to equations of Ida et al (2016)
        !*************************************************************

        planet_pebaccrete = 0.0
        eff_pebble = 0.0
        if (mjeans(ipebrad)>0.0) then

            Hp = sqrt(tstop/alpha(ipebrad))*(1 + tstop/alpha(ipebrad))**(-0.5) *H(ipebrad)

            ! Compute Hill Radius
            rhill_reduce = (mjeans(ipebrad)*mjup/(3.0*mstar))**0.3333
            rhill = rhill_reduce*r(ipebrad)

            ! Compute cross section of pebble flow for accretion

            bcross_reduce = 3*tstop**0.333*rhill_reduce/etadash(irad)

            if(bcross_reduce>1.0) bcross_reduce = 1.0

            bcross_reduce = bcross_reduce*2.0*tstop**0.333*rhill_reduce
            bcross = bcross_reduce*r(ipebrad)

            ! Now compute fragment pebble accretion rate

            planet_pebaccrete = sqrt(8.0/pi)*Hp/bcross

            if(planet_pebaccrete>1.0) planet_pebaccrete = 1.0

            planet_pebaccrete = sqrt(pi/2)*(bcross*bcross/Hp)* &
            mdotpebble/(4.0*pi*r(ipebrad)*tstop*zeta)

            planet_pebaccrete = planet_pebaccrete*Chi*(1.0+ 3*bcross_reduce/(2*Chi*eta(irad)))

            if(planet_pebaccrete>mdotpebble) planet_pebaccrete = mdotpebble

            ! compute pebble accretion efficiency
            eff_pebble = planet_pebaccrete/mdotpebble

        endif

        !*************************************************************************
        ! 4 Find regions where streaming instability permits planetesimal formation
        !**************************************************************************

        !*****************************************************************
        ! 4a Find regions where streaming instability active (rhod/rhog~1)
        !*****************************************************************

        stream_unstable(:) = 0
        inner_radius = .false.
        outer_radius = .false.
        rmin_unstable = 0.0
        rmax_unstable = 0.0
        irmin_unstable = 0
        irmax_unstable = 0

        do jrad =2,ipebrad

           Hp_to_Hg = sqrt(tstop/alpha(ipebrad))*(1 + tstop/alpha(ipebrad))**(-0.5)
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

        ! If streaming zone goes right to pebble front, make sure outer radius is found
       if (inner_radius.eqv. .true. .and. outer_radius.eqv..false.) then           
          rmax_unstable = r(ipebrad)
           irmax_unstable = ipebrad
           outer_radius = .true.
        endif
        
        !****************************************************
        ! 4b compute crossing radius for planetesimal growth
        !****************************************************

        ! Find width of unstable region 

        width_unstable = rmax_unstable - rmin_unstable

        h_unstable = 0.0
        mcross = 0.0

        ! Only do the calculations for non-zero widths
        if(width_unstable > 0.0) then

           h_unstable = H(irmin_unstable)/rmin_unstable

           ! Compute crossing mass at this rpeb (Ormel et al 2017)

           if(eff_pebble>0.0) then
              effpeb = eff_pebble
           else
              effpeb = 0.1
           endif

           mcross = sqrt(3.0*pi*width_unstable*alpha(irmin_unstable)*mdotpebble*effpeb/&
                (rmin_unstable*mdotvisc*gamma1))*h_unstable * h_unstable* Mstar

           !write(*,'(10(1P,e8.1,1X))'), mcross/mjup, &
           !     width_unstable/udist, h_unstable,mdotpebble/mdotvisc, &
           !     alpha(irmin_unstable), mdotpebble, effpeb, mdotvisc, Mstar/umass


        endif


        !****************************
        ! 4. Write output data to files
        !****************************


        ! Check if accretion rate is at a maximum
        if(planet_pebaccrete > planet_accretemax) then
              planet_accretemax = planet_pebaccrete
              mjeans_accretemax = mjeans(ipebrad)
              mdotpebble_accretemax = mdotpebble
              rpeb_accretemax = rpeb
              tpeb_accretemax = tpeb
              mcross_accretemax = mcross
              eff_accretemax = eff_pebble
        endif



           ! Write to file for this rpeb
           ! mdot rpeb --> mdotpebble, r1, r2, Mcross, Mpl
        
           write(30,*) rpeb/udist,tpeb/year, rdotpeb*year/udist, &
                mdotpebble*year/mjup, rmin_unstable/udist, rmax_unstable/udist, &
                mcross/mjup, mjeans(ipebrad),planet_pebaccrete*year/mjup, eff_pebble

     enddo
     close(30)

     ! Now write log file for this mdot
     ! mdot: maximum mdotp, rg(max mdotp), Mcross(max mdotp), efficiency
     write(20,*) mdotvisc*year/umass, rpeb_accretemax/udist, tpeb_accretemax/3.15e7, &
          mdotpebble_accretemax*year/mjup, mcross_accretemax/mjup, &
          mjeans_accretemax, planet_accretemax*year/mjup, &
          eff_accretemax

  enddo
        
close(20)

!END of program

end PROGRAM calc_pebble_accretion
