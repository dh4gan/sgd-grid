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
  read(10,*) zpeb ! Metallicity of disc
  read(10,*) growth_parameter ! growth rate of pebbles t_peb = growth_parameter *(zpeb*omega(irad)^-1
  read(10,*) rmax_peb ! Maximum radius to consider pebble accretion
  read(10,*) mplanet ! Typical planet mass for pebble accretion
  close(10)

  print*, 'Reading file ', inputfile
  print*, 'Outputting pebble accretion data to file ', outputfile

  OPEN(10,file=inputfile ,status='unknown')

  ! Read input file header

  read(10,*) nrad, nmdot, rmin,rmax,mdotmin,mdotmax,Mstar, metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr


  ! Find index of maximum radius for pebble accretion (npebrad)
  dr = (rmax-rmin)/REAL(nrad)
  !dr = dr*udist
  imax_peb = int(rmax_peb/dr) + 1

  ! Write header for log output file
  OPEN(20,file=outputfile,status='unknown')
  write(20,*) nrad,nmdot,nu,rmin,rmax,mdotmin,mdotmax,Mstar,metallicity,gamma_sigma,gamma_omega,irrchoice,Q_irr,T_irr

  ! Loop over accretion rates

  DO imdot = 1, nmdot

     ! Read entire disc over pebble growth radius
     DO irad = 1,imax_peb

        read(10,*) r(irad), mdotvisc(irad), qratio(irad), sigma(irad), cs(irad),omega(irad),T(irad),betac(irad),alpha(irad)

        ! TODO - check units
        H(irad) = cs(irad)/omega(irad)
        rho(irad) = sigma(irad)/(2.0*H(irad))
        !r(irad) = r(irad)*udist
        
     ! Calculate pressure gradient and eta - sub Keplerian parameter

        if(irad>1) then

           dlogrhodr = (log(rho(irad)*cs(irad)*cs(irad)) - log(rho(irad)*cs(irad)*cs(irad)))/(log(r(irad) - log(r(irad-1)))

        else
           dlogrhodr = 0.0
        endif

        eta(irad) = -0.5(H(irad)/(r(irad))**2 * dlogrhodr

         ! Radial velocity of pebbles
        vrpeb(irad) = -2.0*eta(irad)*omega(irad)*r(irad)*tstop/(1.0 + tstop*tstop)

     enddo

     ! Now loop over radius, where radius now refers to rpeb

     do ipebrad =1,npebrad

        rpeb = r(ipebrad)

        ! Calculate growth rate of pebbles

        ! Pebble growth timescale at this radius
        tpeb = growth_param/(zpeb*omega(ipebrad))

        ! Pebble front growth rate
        rdotpeb = 0.6666*(G*Mstar*zpeb*zpeb/(tpeb*growth_param*growth_param))**(0.333)

       
        ! Compute mdotpebble
        mdotpebble = 2.0*pi*rpeb*rdotpeb*zpeb*sigma(ipebrad)

        ! Find regions where streaming instability active
        stream_unstable(:) = 0
        inner_radius = .false.
        outer_radius = .false.

        do jrad =1,ipebrad

           Hp_to_Hg = sqrt(tstop/alpha(jrad))
           ! Compute rhop/rhog interior to pebble radius
           rhop_rhog(jrad) = mdotpebble/(Hp_to_Hg*sigma(jrad)*2.0*pi*r(jrad)*vrpeb(jrad))

           ! If rhop/rhog >1, mark this radius as streaming unstable region
           
           if(rhop_rhog(jrad)>1.0) then
              stream_unstable(jrad) = 1

              ! Record minimum and maximum values of streaming regions (assuming single region only
              if(inner_radius.eqv..false.) then
                 rmin_unstable = r(jrad)
                 inner_radius = .true.
              endif
           else
              if(stream_unstable(jrad-1)==1) then
                 rmax_unstable = r(jrad)
                 outer_radius = .true.
              endif

        enddo

        ! Find width of unstable region 

        width_unstable = rmax_unstable - rmin_unstable

        ! Compute surface density of pebbles in here
        sigma_p = mdotpebble/(2.0*pi*r(irmin_unstable)*vr(irmin_unstable))

        ! Compute Pebble Accretion Rate for M= 1 Mjup
        ! Mdot = 2 R_H^2 omega tstop *sigma p

        rhill = r(irmin_unstable)*(mplanet/3.0*mstar)**0.333

        planet_pebaccrete = 2.0* rhill*rhill*omega*tstop*sigma_p

        ! compute pebble accretion efficiency
        eff_pebble = planet_pebaccrete/mdotpebble

        ! Compute crossing mass at this rpeb

        mcross = sqrt(3.0*pi*width_unstable*alpha(irmin_unstable)*mdotpebble*eff_pebble/(mdotvisc*gamma1)*(H(irmin_unstable)/r(irmin_unstable))**2 * Mstar

        ! Write to file for this rpeb
        ! mdot rpeb --> mdotpebble, r1, r2, Mcross, Mpl
        
     enddo

     ! Now write file for this mdot
     ! mdot: maximum mdotp, rg(max mdotp), Mcross(max mdotp), efficiency
  enddo
        

!END of program

end PROGRAM calc_pebble_accretion
