PROGRAM calc_pebble_accretion
  ! Program calculates observables for a grid of pre-generated
  ! self-gravitating disc models
  ! Assumes a dust temperature and opacity properties to do so

  use eosdata
  implicit none

  integer :: irad,imdot,irrchoice
  integer :: nrad, nmdot

  real :: mdotmax,mdotmin, mdotvisc, metallicity, mstar
  real :: qratio, T, gamma_sigma,gamma_omega
  real :: Q_irr, T_irr,dr, rmax,rmin, betac

  real, parameter :: mu = 2.4
 
  real, parameter :: tolerance = 1.0e-5
  real, parameter :: Qfrag = 2.0
 
  real, parameter :: gamma1 = 4.0
  real, parameter :: vfrag = 1.0e3 ! Empirically determined fragmentation velocity (cm/s)

  integer, allocatable,dimension(:) :: stream_unstable
  real, allocatable, dimension(:) :: r, sigma,cs,omega,alpha,mjeans
  real, allocatable, dimension(:) :: H, eta, vrpeb, rhogas,etadash,tstop,grainsize
  real,allocatable,dimension(:) :: tstop_frag,maxgrainsize

  integer :: ipebrad,npebrad, irmin_unstable,irmax_unstable, jrad
  real :: zpeb, beta_peb, rmax_peb, mplanet, Hp_to_Hg, Hp, dlogrhodr
  real :: grainsize_in,tstop_in, rhosolid,mfp
  
  real :: rpeb, rdotpeb, tpeb, mdotpebble, rhill, rhop_rhog
  real :: rmin_unstable, rmax_unstable, width_unstable, h_unstable
real :: bcross_reduce, bcross, rhill_reduce, zeta, Chi
  real :: planet_pebaccrete, mcross, eff_pebble, effpeb
  real :: rpeb_accretemax, tpeb_accretemax, mdotpebble_accretemax, mjeans_accretemax
  real :: planet_accretemax, mcross_accretemax, eff_accretemax

  logical :: inner_radius
  logical :: outer_radius

  character(1) :: fixparam
  character(100) ::inputfile, outputfile, outputlog

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
  READ(10,*) outputfile ! Output file name
  read(10,*) fixparam ! Fix either tstop or grain size
  READ(10,*) tstop_in ! Input Stopping time (dimensionless)
  read(10,*) grainsize_in ! Input Grain size (cm)
  read(10,*) rhosolid ! Density of grains (g cm ^-3)
  read(10,*) zpeb ! Metallicity of disc
  read(10,*) beta_peb ! growth rate of pebbles t_peb = beta_peb *(zpeb*omega(irad)^-1
  read(10,*) rmax_peb ! Maximum radius to consider pebble accretion  (AU)
  read(10,*) mplanet ! Typical planet mass for pebble accretion (Jupiter masses)
  close(10)

  print*, 'Reading file ', inputfile
  print*, 'Outputting pebble accretion data to file ', outputfile

  outputlog = trim(outputfile)//".max"
  print*, "Maximal pebble accretion data per gas accretion rates to be written to ", trim(outputlog)

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
  allocate(grainsize(nrad))
  allocate(tstop(nrad))
  allocate(tstop_frag(nrad))
  allocate(maxgrainsize(nrad))

  ! Find index of maximum radius for pebble accretion (npebrad)
  dr = (rmax-rmin)/REAL(nrad)
  !dr = dr*udist
  npebrad = int(rmax_peb/dr) + 1

  if(npebrad > nrad) npebrad = nrad


  ! Write header for main output file
  open(30, file=outputfile, status='unknown')
  write(30,*) npebrad,nmdot,rmin,rmax,Mstar,fixparam,zpeb,beta_peb,&
       metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

  ! Write header for global maxima output file
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

     ! Fill tstop, grainsize arrays with default values (changed later)
     tstop(:) = tstop_in 
     grainsize(:) = grainsize_in
     tstop_frag(:) = 0.0
     maxgrainsize(:) = 0.0

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
                
        !************************************************************************
        ! 1a  Calculate pressure gradient and eta - sub Keplerian parameter
        !************************************************************************

        if(irad>1) then
           dlogrhodr = (log(rhogas(irad)*cs(irad)*cs(irad)) - log(rhogas(irad-1)*cs(irad-1)*cs(irad-1)))/&
                (log(r(irad) - log(r(irad-1))))
        else
           dlogrhodr = 0.0
        endif

        eta(irad) = -0.5*(H(irad)*H(irad))/(r(irad)*r(irad)) * dlogrhodr

        ! Compute etadash function (for pebble accretion rates)
        etadash(irad) = Chi*eta(irad)


     
!        write(76,*) r(irad)/udist, rhogas(irad)*cs(irad)*cs(irad), &
!             H(irad)/udist, rhogas(irad), dlogrhodr, eta(irad), (r(irad)/vrpeb(irad))/3.15e7

        !****************************************************************************
        ! 1b Compute either grain size or stopping time depending on inputs
        !****************************************************************************

        ! Compute mean free path
        mfp = sqrt(twopi)*mu*mH*H(irad)/(coll_H*sigma(irad))

        ! Compute maximum stopping time (from grain fragmentation)       
        ! (see Birnstiel et al 2009, A&A, 503, L5-L8)

        tstop_frag(irad) = vfrag*vfrag/(alpha(irad)*cs(irad)*cs(irad))
        call calc_grainsize(maxgrainsize(irad),mfp,rhosolid,sigma(irad),tstop_frag(irad))

        ! If stopping time fixed
        if(fixparam=='t') then

           if(tstop(irad)>tstop_frag(irad)) tstop(irad) = tstop_frag(irad)
           call calc_grainsize(grainsize(irad),mfp,rhosolid,sigma(irad),tstop(irad))

           ! Otherwise if grain size fixed
        else

           call calc_tstop(tstop(irad),mfp,rhosolid,sigma(irad),grainsize(irad))

           ! if tstop > tstop_frag, must recompute both tstop and grainsize
           if(tstop(irad) > tstop_frag(irad)) then
              tstop(irad) = tstop_frag(irad)
              call calc_grainsize(grainsize(irad),mfp,rhosolid,sigma(irad),tstop(irad))
           endif

        endif

        ! Now we have tstop, compute radial velocity of pebbles
        vrpeb(irad) = 2.0*eta(irad)*omega(irad)*r(irad)*tstop(ipebrad)/(1.0 + tstop(irad)*tstop(irad))
     enddo

     !******************************************************************
     ! 2. Now compute pebble accretion properties as a function of rpeb
     !******************************************************************
   

     ! Variables store maximum values at this gas accretion rate
     planet_accretemax = 0.0
     mjeans_accretemax = 0.0
     mdotpebble_accretemax = 0.0
     rpeb_accretemax = 0.0
     tpeb_accretemax = 0.0
     mcross_accretemax = 0.0
     eff_accretemax = 0.0

     ! convert mdot into cgs
     mdotvisc = mdotvisc*umass/year

     ! Now loop again over radius, where radius now refers to rpeb

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

       
        ! Initialise parameters for calculation
        planet_pebaccrete = 0.0
        eff_pebble = 0.0
        mcross = 0.0
        h_unstable = 0.0
        stream_unstable(:) = 0
        inner_radius = .false.
        outer_radius = .false.
        rmin_unstable = 0.0
        rmax_unstable = 0.0
        irmin_unstable = 0
        irmax_unstable = 0
        width_unstable = 0.0

        !*************************************************************
        ! 3 Compute the pebble accretion rate of any fragment present
        ! (according to equations of Ida et al (2016)
        !*************************************************************
           
        if (mjeans(ipebrad)>0.0) then

           ! Convenience functions for calculating pebble accretion 
           ! (Ida et al 2016)
              
           zeta = 1/(1+tstop(ipebrad)*tstop(ipebrad))
           Chi = sqrt((1.0 + 4.0*tstop(ipebrad)*tstop(ipebrad)))/(1.0+tstop(ipebrad)*tstop(ipebrad))
              
           Hp = sqrt(tstop(ipebrad)/alpha(ipebrad))*(1 + tstop(ipebrad)/alpha(ipebrad))**(-0.5) *H(ipebrad)
              
           ! Compute Hill Radius
           rhill_reduce = (mjeans(ipebrad)*mjup/(3.0*mstar))**0.3333
           rhill = rhill_reduce*r(ipebrad)
              
           ! Compute cross section of pebble flow for accretion
              
           bcross_reduce = 3*tstop(ipebrad)**0.333*rhill_reduce/etadash(irad)
           
           if(bcross_reduce>1.0) bcross_reduce = 1.0
           
           bcross_reduce = bcross_reduce*2.0*tstop(ipebrad)**0.333*rhill_reduce
           bcross = bcross_reduce*r(ipebrad)
           
           ! Now compute disc fragment's pebble accretion rate
           
           planet_pebaccrete = sqrt(8.0/pi)*Hp/bcross
           
           if(planet_pebaccrete>1.0) planet_pebaccrete = 1.0
           
           planet_pebaccrete = sqrt(pi/2)*(bcross*bcross/Hp)* &
                mdotpebble/(4.0*pi*r(ipebrad)*tstop(ipebrad)*zeta)
           
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
        
        
        do jrad =2,ipebrad
           
           Hp_to_Hg = sqrt(tstop(ipebrad)/alpha(ipebrad))*(1 + tstop(ipebrad)/alpha(ipebrad))**(-0.5)
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
        
        
        
        ! Write to main output file
        write(30,*) mdotvisc*year/umass, rpeb/udist,grainsize(ipebrad),tstop(ipebrad),&
             tstop(ipebrad)/tstop_frag(ipebrad),maxgrainsize(ipebrad),tpeb/year, &
             rdotpeb*year/udist, mdotpebble*year/mjup, &
             rmin_unstable/udist, rmax_unstable/udist, &
             mcross/mjup, mjeans(ipebrad),planet_pebaccrete*year/mjup, eff_pebble
        
     enddo
     
     ! Now write to global maxima file for this mdot
     
     write(20,*) mdotvisc*year/umass, rpeb_accretemax/udist, tpeb_accretemax/3.15e7, &
          mdotpebble_accretemax*year/mjup, mcross_accretemax/mjup, &
          mjeans_accretemax, planet_accretemax*year/mjup, &
          eff_accretemax
     
  enddo
  
  close(20)
  close(30)
  
end PROGRAM calc_pebble_accretion



subroutine calc_tstop(tstop,mfp,rhosolid,sigmagas,grainsize)
  ! Calculates the stopping time for Epstein and Stokes drag regimes
  ! (Ida et al 2016)
  use eosdata, only: roottwopi
  implicit none

  real, intent(in) :: mfp,rhosolid,sigmagas,grainsize
  real, intent(inout) :: tstop

  if(grainsize<9.0*mfp/4.0) then
     tstop = roottwopi*rhosolid*grainsize/sigmagas
  else
     tstop = 4.0*roottwopi*rhosolid*grainsize*grainsize &
          /(9.0*mfp*sigmagas)
  endif

end subroutine calc_tstop


subroutine calc_grainsize(grainsize,mfp,rhosolid,sigmagas,tstop)
  ! Calculates the grain size for a given stopping time
  ! (Ida et al 2016)

  use eosdata, only: roottwopi
  implicit none

  real, intent(in) :: mfp,rhosolid,sigmagas,tstop
  real, intent(inout) :: grainsize

   ! Calculate assuming Epstein drag, then check
  grainsize = tstop*sigmagas/(roottwopi*rhosolid)

  ! If in Stokes regime, recompute
  if(grainsize >= 9.0*mfp/4.0) then
     grainsize = 9.0*mfp*sigmagas*tstop/(4.0*roottwopi*rhosolid)
     grainsize = sqrt(grainsize)
  endif

end subroutine calc_grainsize
