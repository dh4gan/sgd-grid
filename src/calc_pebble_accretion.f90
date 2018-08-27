PROGRAM calc_pebble_accretion
  ! Program reads in output from sgd_grid
  ! Computes properties of the disc related to pebble accretion:
  ! i) Where streaming instability may occur, and
  ! ii) Pebble accretion rates of fragments
  ! Assumes either a fixed grain size or fixed dimensionless stopping time

  use eosdata
  implicit none

  integer :: irad,imdot,irrchoice
  integer :: nrad, nmdot

  real :: mdotmax,mdotmin, mdotvisc, metallicity, mstar
  real :: qratio, T, gamma_sigma,gamma_omega
  real :: Q_irr, T_irr,dr, rmax,rmin, betac, vrgas, mpebble
  real :: gapcriterion, tcross,tgap,tmig1,aspectratio,dustaspectratio,mratio
  real :: tstop_turb, mturb,maxgrow, miso_peb, effpeb
  real :: percentcount, displaypercent,increment

  real, parameter :: mu = 2.4

  real, parameter :: gamma1 = 1.0
  real, parameter :: vfrag = 1.0e3 ! Empirically determined fragmentation velocity (cm/s)

  real, parameter :: alpha_max = 1.0e20 ! Maximum turbulence for streaming instability

  integer, allocatable,dimension(:) :: stream_unstable
  real, allocatable, dimension(:) :: r, sigma,cs,omega,alpha,mjeans,Hp_to_Hg,rhop_rhog
  real, allocatable, dimension(:) :: H, eta, vrpeb, rhogas,etadash,tstop,grainsize
  real,allocatable,dimension(:) :: tstop_frag,maxgrainsize, sigma_peb, sigma_peb_max,mpebtot

  integer :: ipebrad,npebrad, irmin_unstable,irmax_unstable
  real :: zpeb, fpeb,beta_peb, rmax_peb, mplanet, Hp, dlogrhodr
  real :: grainsize_in,tstop_in, rhosolid,mfp, actual_width

  real :: rpeb, rdotpeb, tpeb, mdotpebble, rhill
  real :: rmin_unstable, rmax_unstable, width_unstable, h_unstable
  real :: bcross_reduce, bcross, rhill_reduce, zeta, Chi
  real :: planet_pebaccrete, mcross, eff_pebble
  real :: rpeb_accretemax, tpeb_accretemax, mdotpebble_accretemax, mjeans_accretemax
  real :: planet_accretemax, mcross_accretemax, eff_accretemax

  logical :: inner_radius
  logical :: outer_radius

  character(1) :: fixparam
  character(100) ::inputfile, inputprefix,outputfile, outputlog

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
  read(10,*) inputprefix
  read(10,*) fixparam ! Fix either tstop or grain size
  READ(10,*) tstop_in ! Input Stopping time (dimensionless)
  read(10,*) grainsize_in ! Input Grain size (cm)
  read(10,*) rhosolid ! Density of grains (g cm ^-3)
  read(10,*) zpeb ! Metallicity of disc
  read(10,*) fpeb ! Fraction of disc in pebbles
  read(10,*) effpeb ! Pebble accretion efficiency
  read(10,*) beta_peb ! growth rate of pebbles t_peb = beta_peb *(zpeb*omega(irad)^-1
  read(10,*) rmax_peb ! Maximum radius to consider pebble accretion  (AU)
  close(10)

  inputfile = trim(inputprefix)//'.sgdmodel'
  print*, 'Reading file ', inputfile
  outputfile = trim(inputprefix)//'.pebble'
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
  allocate(Hp_to_Hg(nrad))
  allocate(rhop_rhog(nrad))
  allocate(sigma_peb(nrad))
  allocate(mpebtot(nrad))
  allocate(sigma_peb_max(nrad))
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
  dr = dr*udist
  npebrad = nrad

  if(npebrad > nrad) npebrad = nrad


  ! Write header for main output file
  open(30, file=outputfile, status='unknown')
  write(30,*) npebrad,nmdot,rmin,rmax,Mstar,fixparam,zpeb,beta_peb,&
       metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

  ! Write header for global maxima output file
  OPEN(20,file=outputlog,status='unknown')
  write(20,*) nrad,nmdot,rmin,rmax,mdotmin,mdotmax,Mstar,&
       tstop,zpeb,fpeb,beta_peb,&
       metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr


  !*******************************************************************
  ! 1 Compute disc properties needed for pebble accretion calculations
  !********************************************************************

  increment = 10.0
  displaypercent = 10.0

  ! Loop over accretion rates
  DO imdot = 1, nmdot

     percentcount = real(imdot)*100.0/real(nmdot)
     if(percentcount > displaypercent) then
        print('(F5.0,A)'), displaypercent, '% complete'
        displaypercent = displaypercent +increment
     endif


     r(:) = 0.0
     sigma(:) = 0.0
     cs(:) = 0.0
     H(:) = 0.0
     Hp_to_Hg(:) = 0.0
     rhop_rhog(:) = 0.0
     rhogas(:) = 0.0
     omega(:) = 0.0
     alpha(:) = 0.0
     eta(:) = 0.0
     vrpeb(:) = 0.0
     mjeans(:) = 0.0
     stream_unstable(:) = 0
     
     inner_radius = .false.
     outer_radius = .false.
     rmin_unstable = 0.0
     rmax_unstable = 0.0
     irmin_unstable = 0
     irmax_unstable = 0
     width_unstable = 0.0


     ! Fill tstop, grainsize arrays with default values (changed later)
     tstop(:) = tstop_in 
     grainsize(:) = grainsize_in
     tstop_frag(:) = 0.0
     maxgrainsize(:) = 0.0

     ! Read entire disc at this accretion rate
     DO irad = 1,nrad

        read(10,*) r(irad), mdotvisc, qratio, sigma(irad), cs(irad),omega(irad),&
             T,betac,alpha(irad), mjeans(irad)

        ! Be careful with units here
        H(irad) = (cs(irad)/omega(irad)) ! in units of cm
        rhogas(irad) = sigma(irad)/(2.0*H(irad))
        r(irad) = r(irad)*udist

        

        mpebtot(irad) = fpeb*zpeb*mstar*qratio

        if(irad>1) then
           sigma_peb_max(irad) = (mpebtot(irad)-mpebtot(irad-1))/(twopi*r(irad)*dr)
        else
           sigma_peb_max(irad) = 0.0
        endif

     enddo

     !print*, 'Total pebble mass: ',mpebtot(nrad)/mearth ,' Earth masses'

     mpebble = 0.0
     do irad = 1,nrad

        !**********************************************************************
        ! 1a  Calculate pressure gradient and eta - sub Keplerian parameter
        !**********************************************************************

        if(irad>1) then
           dlogrhodr = (log(rhogas(irad)*cs(irad)*cs(irad)) - log(rhogas(irad-1)*cs(irad-1)*cs(irad-1)))/&
                (log(r(irad) - log(r(irad-1))))
        else
           dlogrhodr = 0.0
        endif

        eta(irad) = abs(0.5*(H(irad)*H(irad))/(r(irad)*r(irad)) * dlogrhodr)

        !*********************************************************************
        ! 1b Compute either grain size or stopping time depending on inputs
        !**********************************************************************

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

        Chi = sqrt((1.0 + 4.0*tstop(irad)*tstop(irad)))/(1.0+tstop(irad)*tstop(irad))

        ! Compute etadash function (for pebble accretion rates)
        etadash(irad) = Chi*eta(irad)

        ! Pebble radial velocity = drift velocity + viscous radial velocity 
        ! (Takeuchi & Lin 2002)

        ! Drift velocity first
        vrpeb(irad) = 2.0*eta(irad)*omega(irad)*r(irad)*tstop(irad)/(1.0 + tstop(irad)*tstop(irad))

        ! now gas viscous radial velocity

	if(irad>1) then
        vrgas = (sqrt(r(irad)*alpha(irad)*cs(irad)*H(irad)*sigma(irad)-  &
                sqrt(r(irad-1)*alpha(irad-1)*cs(irad-1)*H(irad-1)*sigma(irad-1))))/&
                (r(irad)-r(irad-1))
	else
	vrgas = 0.0
	endif

        vrgas = -3.0*vrgas/(sqrt(r(irad))*sigma(irad))

        vrpeb(irad) = abs(vrgas/(1.0+tstop(irad)*tstop(irad)) + vrpeb(irad))

        rpeb = r(irad)

        !**************************************************
        ! 1c Compute susceptibility to streaming instability
        !***************************************************

        ! Need to compute mdotpebble first to do this here

        ! Pebble growth timescale at this radius
        tpeb = beta_peb/(zpeb*omega(irad))

        ! Pebble front growth rate
        rdotpeb = 0.6666*(G*Mstar*zpeb*zpeb/(tpeb*beta_peb*beta_peb))**(0.333)

        ! Compute mdotpebble (g s^-1)
        mdotpebble = 2.0*pi*rpeb*rdotpeb*zpeb*fpeb*sigma(irad)

        ! Scale height ratio of dust to gas (Dubrulle et al 1995)
        Hp_to_Hg(irad) = sqrt(alpha(irad)/(alpha(irad)+tstop(irad)))

        ! Surface density of pebbles
        if(vrpeb(irad)>1.0e-30) then
           sigma_peb(irad) = (rdotpeb/vrpeb(irad))*zpeb*fpeb*sigma(irad)
        else
           sigma_peb(irad) = 0.0
        endif

        ! Ensure surface density of pebbles can't exceed local maximum
        if(sigma_peb(irad)>sigma_peb_max(irad))then
           sigma_peb(irad) = sigma_peb_max(irad)
        endif
        
        if(Hp_to_Hg(irad)*sigma(irad) >1.0e-30) then
           ! Ratio of dust to gas density at this radius
           rhop_rhog(irad) = sigma_peb(irad)/(sigma(irad)*Hp_to_Hg(irad))
        else
           rhop_rhog(irad) = 0.0
        endif
        
        ! Double check if total pebble supply used up - compute pebble mass generated so far
        mpebble = mpebble + twopi*r(irad)*sigma_peb(irad)*dr
        
        ! If pebble supply exhausted, then set pebble density to zero
        if(mpebble>(mpebtot(nrad)/fpeb)) then
           rhop_rhog(irad) = 0.0
           Hp_to_Hg(irad) = 0.0
        endif

        ! If rhop/rhog >1 and not in fragmentation zone, 
        ! mark this radius as streaming unstable region           
        if(rhop_rhog(irad)>1.0 .and.mjeans(irad)<1.0e-30 .and.alpha(irad)<alpha_max) then
           stream_unstable(irad) = 1

           ! Record minimum and maximum values of streaming regions 
           ! (assuming single region only)
           ! If in a fragmentation zone, curtail streaming region

           if(inner_radius.eqv..false.) then
              rmin_unstable = r(irad)
              irmin_unstable = irad
              inner_radius = .true.
           endif
        else
           if(irad>1) then
              if(stream_unstable(irad-1)==1) then
                 rmax_unstable = r(irad)
                 irmax_unstable = irad
                 outer_radius = .true.
                 width_unstable = rmax_unstable-rmin_unstable
              endif
           endif
           endif

        enddo
  
       ! print*, rmin_unstable/udist, rmax_unstable/udist, outer_radius,sum(stream_unstable)>0,&
       !      (outer_radius .eqv. .false.) .and. sum(stream_unstable)>0


        if((outer_radius .eqv. .false.) .and. sum(stream_unstable)>0) then
           rmax_unstable = r(nrad)
           irmax_unstable = nrad
           width_unstable = rmax_unstable - rmin_unstable
        endif

        !print*, rmin_unstable/udist, rmax_unstable/udist
        !if(sum(stream_unstable)<1) STOP

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
        mdotpebble = 2.0*pi*rpeb*rdotpeb*zpeb*fpeb*sigma(ipebrad)

        ! Initialise parameters for calculation
        planet_pebaccrete = 0.0
        eff_pebble = 0.0
        mcross = 0.0
        mturb = 0.0
        h_unstable = 0.0
       

        !*************************************************************
        ! 3 Compute the pebble accretion rate of any fragment present
        ! (according to equations of Ida et al (2016)
        !*************************************************************

        ! Compute some useful expressions out front
        Hp = Hp_to_Hg(ipebrad)*H(ipebrad)
        aspectratio = H(ipebrad)/r(ipebrad)
        dustaspectratio = Hp/r(ipebrad)
        mratio = mjeans(ipebrad)*mjup/mstar


        ! Pebble isolation mass
        miso_peb = (dustaspectratio*dustaspectratio*dustaspectratio)*mstar


        if (mjeans(ipebrad)>0.0) then

           ! Convenience functions for calculating pebble accretion 
           ! (Ida et al 2016)

           zeta = 1/(1+tstop(ipebrad)*tstop(ipebrad))
           Chi = sqrt((1.0 + 4.0*tstop(ipebrad)*tstop(ipebrad)))/(1.0+tstop(ipebrad)*tstop(ipebrad))

           ! Compute Hill Radius
           rhill_reduce = (mjeans(ipebrad)*mjup/(3.0*mstar))**0.3333
           rhill = rhill_reduce*r(ipebrad)

           !***********************************************************
           ! Does fragment drive a gap in pebbles? If so, no accretion
           !***********************************************************

           gapcriterion = 0.75*Hp/rhill + &
                50.0*alpha(ipebrad)*dustaspectratio**2/mratio

           ! Also check that migration slow enough to allow gap formation

           ! Type I migration timescale (Ormel et al 2017)
           tmig1 = mstar*aspectratio*aspectratio/(gamma1*mratio*sigma(ipebrad)*r(ipebrad)*r(ipebrad)*omega(ipebrad))

           ! Gap crossing timescale
           tcross = 2.5*rhill*tmig1/r(ipebrad)
           
           ! Gap formation timescale
           tgap = 1.0e2*(dustaspectratio**5)/(omega(ipebrad)*mratio*mratio)          

           
           ! If a gap is opened in the pebbles, then no pebble accretion
           if(gapcriterion<1.0 .and. tgap< tcross) then
              planet_pebaccrete = 0.0
           
           else

              ! Compute cross section of pebble flow for accretion

              bcross_reduce = 3*tstop(ipebrad)**0.333*rhill_reduce/etadash(ipebrad)

              if(bcross_reduce>1.0) bcross_reduce = 1.0

              bcross_reduce = bcross_reduce*2.0*tstop(ipebrad)**0.333*rhill_reduce
              bcross = bcross_reduce*r(ipebrad)

              ! Now compute disc fragment's pebble accretion rate

              planet_pebaccrete = sqrt(8.0/pi)*Hp/bcross

              if(planet_pebaccrete>1.0) planet_pebaccrete = 1.0

              planet_pebaccrete = sqrt(pi/2)*(bcross*bcross/Hp)* &
                   mdotpebble/(4.0*pi*r(ipebrad)*tstop(ipebrad)*zeta)

              planet_pebaccrete = planet_pebaccrete*Chi*(1.0+ 3*bcross_reduce/(2*Chi*eta(ipebrad)))

              if(planet_pebaccrete>mdotpebble) planet_pebaccrete = mdotpebble

              if(abs(planet_pebaccrete*year/mjup)>1.0e10) then
                 print*, planet_pebaccrete, Hp, bcross, bcross_reduce,eta(ipebrad)
              endif

           endif
           ! measure pebble accretion efficiency
           eff_pebble = planet_pebaccrete/mdotpebble

        endif

        !****************************************************
        ! 4 compute crossing radius for planetesimal growth
        !****************************************************
    
        mcross = 0.0
        actual_width = 0.0
        

        ! Only do the calculations for non-zero widths
        ! Do calculation assuming ipebrad = inner distance of streaming
        if(width_unstable > 1.0e-40) then
           if(r(ipebrad)>rmin_unstable .and. r(ipebrad)<rmax_unstable) then
           
           actual_width = rmax_unstable-r(ipebrad)
           h_unstable = H(ipebrad)/r(ipebrad)                    

           ! Compute crossing mass at this rpeb (Ormel et al 2017)
           ! with our imposed pebble accretion efficiency (effpeb parameter)
           
           mcross = sqrt(3.0*pi*actual_width*alpha(ipebrad)*mdotpebble*effpeb/&
                (rmin_unstable*mdotvisc*gamma1))*h_unstable * h_unstable* Mstar
           
          
           ! If crossing mass exceeds local isolation mass, restrict it here
          
           if(mcross > miso_peb) mcross = miso_peb

           ! If entire pebble disc required to create embryo, 
           ! ensure mass cannot exceed total pebble mass (at given efficiency)

           if(mcross > mpebtot(nrad)*effpeb) then
              !print*, 'Pebble disc mass exceeded: ',mpebtot(nrad)*effpeb/mearth
              mcross = mpebtot(nrad)*effpeb
           endif
           
           ! Compute maximum growable mass in turbulent disc

           ! Find maximum allowed stopping time for bidisperse grain population
           ! (Booth & Clarke 2016, MNRAS, 458, 2676)

           tstop_turb = vfrag/(cs(ipebrad)*4.47*sqrt(alpha(ipebrad))) + tstop(ipebrad)

           ! Compute mean free path
           mfp = sqrt(twopi)*mu*mH*H(ipebrad)/(coll_H*sigma(ipebrad))
           
           call calc_grainsize(maxgrow,mfp,rhosolid,sigma(ipebrad),tstop_turb)

           ! maximum turbulent mass
           mturb = 4.0*pi*maxgrow*maxgrow*maxgrow*rhosolid

           !print*, mcross/mjup, miso_peb/mjup, dustaspectratio, aspectratio

           if(mturb>1.0e7) then
              print*, ipebrad, mturb/1.0e3, maxgrow, tstop_turb, tstop(ipebrad), tstop_turb/tstop(ipebrad)
           endif
           !write(*,'(11(1P,e8.1,1X))'), mcross/mjup, &
           !     actual_width/udist, gamma1,h_unstable,mdotpebble/mdotvisc, &
           !     alpha(irmin_unstable), mdotpebble, effpeb, mdotvisc, Mstar/umass

           endif
        endif

        !****************************
        ! 4. Write output data to files
        !****************************

        !print*, r(ipebrad)/udist, Hp_to_Hg(ipebrad)*H(ipebrad)/r(ipebrad)
        ! Check if fragment pebble accretion rate is at a maximum
        ! Record maximal values

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
             rdotpeb*year/udist, mdotpebble*year/mjup,  &
             Hp_to_Hg(ipebrad), rhop_rhog(ipebrad),vrpeb(ipebrad),&
             actual_width/udist, rmin_unstable/udist, &
             mcross/mearth, maxgrow, miso_peb/mearth, mjeans(ipebrad),planet_pebaccrete*year/mjup, eff_pebble

     enddo

     ! Now write to global maxima file for this mdot

     write(20,*) mdotvisc*year/umass, rpeb_accretemax/udist, tpeb_accretemax/3.15e7, &
          mdotpebble_accretemax*year/mjup, mcross_accretemax/mearth, &
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
