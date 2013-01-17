PROGRAM selfgravdisc_modelgrid
  !	Program does a parameter survey of self-gravitating disc models
  !	Writes disc properties and fragmentation outcomes...
  !	******** Assumes fixed Q
  !	******** Surveys the mdot - r parameter space only
  !	Finds local Jeans mass ==> potential fragment masses

  use eosdata
  implicit none

  integer :: frag,selfgrav,irad,imdot,irrchoice
  integer :: nrad, nmdot

  real, parameter :: pi = 3.14159265285
  real, parameter :: twopi = 2.0*pi
  real, parameter :: sigma_SB = 5.67e-5
  real,parameter :: c = 2.99e10
  real, parameter :: tolerance = 1.0e-5
  real, parameter :: Qfrag = 2.0
  real, parameter :: mjup = 1.8986e30
  real,parameter :: Q_irrcrit = 1.99
  real,parameter :: alpha_sat = 0.1 ! Value of alpha at which self-gravity torque saturates
  real, parameter :: gamma_Jcrit = -5 ! Minimum GammaJ value to be exceeded for fragmentation
  real, parameter :: gamma_Qcrit = 10 ! abs(gamma_Q) must be greater than this for self-gravitating disc

  real :: Mstar, rmax,rmin,dr,r
  real :: mdotmin, mdotmax, mdotvisc, dmdot
  real :: sigma, omega, cs,rhomid,alpha, gamma, mu, T, kappa
  real :: mtot,qratio
  real :: dT,sigma_old,ntries,mdot_try,fine
  real :: deltasigma, metallicity

  real :: gamma_J, gamma_omega, gamma_sigma, gamma_Q

  real :: mdot_eq, alpha_eq
  real :: tau,H,betac,dsig
  real :: mjeans, ljeans, rhill
  real :: cs_irr, T_irr, Q_irr, Msol, rAU
  character(100) :: outputfile

  ! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"  
  print*, "  SURVEYING 1D SELF-GRAVITATING DISC MODELS  " 
  print*, "     Created by D. Forgan, May 2011 "
  print*, "     Majorly Revised, Jan 2013  "
  print*, "	                                        	  "
  print*, "-----------------------------------------------"
  print*, " "
  print*, " "
  print*, " input files: ./survey_jeans.params ./myeos.dat" 
  print*, " "
  print*, "-----------------------------------------------"
  print*, " "

  OPEN(10,file='survey_jeans.params', status='unknown')
  READ(10,*) outputfile
  READ(10,*) nrad
  READ(10,*) nmdot
  READ(10,*) Mstar
  READ(10,*) rmin
  READ(10,*) rmax
  read(10,*) mdotmin
  read(10,*) mdotmax
  read(10,*) metallicity
  read(10,*) gamma_sigma
  read(10,*) gamma_omega
  read(10,*) irrchoice
  read(10,*) Q_irr
  read(10,*) T_irr

  close(10)

  print*, 'Outputting to file ', outputfile
  print*, 'Mstar (msol) ', Mstar
  print*, 'Radial Limits: ', rmin, rmax
  print*, 'log(mdotmin), log(mdotmax): ', mdotmin,mdotmax
  print*, 'Metallicity (relative to solar): ', metallicity

  if(irrchoice==-1) then
     print*, 'No Irradiation'
  else if(irrchoice==0) then
     print*, 'Irradiation: fixed Qirr of ', Q_irr
  else if(irrchoice==1) then
     print*, 'Irradiation: fixed Tirr of ', T_irr
  else if (irrchoice==2) then
     print*, 'Irradiation: Ida Lin prescription'
  endif

  rmin = rmin*udist
  rmax = rmax*udist

  dr = (rmax-rmin)/REAL(nrad)			
  dmdot = (mdotmax-mdotmin)/REAL(nmdot)							
  Mstar = Mstar*umass	! star mass in cgs		

  ! gammamuT columns:
  ! 1) gamma
  ! 2) mu
  ! 3) T
  ! 4) Opacity
  ! 5) Mass Weighted Opacity (Stamatellos et al 2007)

  allocate(gammamuT(5)	)	
  !	Read in EoS

  CALL eosread

  OPEN(10,file=outputfile ,status='unknown')
  write(10,*) nrad, nmdot, Mstar/umass, metallicity, gamma_sigma,gamma_omega, irrchoice, Q_irr, T_irr

  ! Loop over parameters

  DO imdot = 1, nmdot

     mdotvisc = mdotmin + (imdot-1)*dmdot				

     print*, 'Testing mdot ', imdot, 10.0**mdotvisc

     mdotvisc = (10.0**mdotvisc)*umass/3.15e7 ! convert to cgs

     mtot = 0.0

     ! Loop over radii
     DO irad = 1,nrad

	r = rmin +(irad-1)*dr 
        sigma = 0.0
        omega = 0.0
        T = 0.0
        cs = 0.0
        mjeans = 0.0
	betac = 0.0
	alpha = 0.0
	alpha_eq = 0.0
	mdot_eq = 0.0				


 ! mjeans always has the same prefactor - do this first
        mjeans = 4.0*1.412*pi*pi*pi/(3.0*G)

        ! Set flags for selfgravitating, fragmenting 
	frag =0
        selfgrav = 0

        ! Calculate Omega

	omega = sqrt(G*Mstar/r**3)

 ! Begin iteration process

 ! Guess Sigma (start from previous answer)

	IF(irad==1) sigma_old = 50000.0

	sigma = 2.0*sigma_old

	dT = 1.0e30
	ntries = 0
        fine = 0.01
	DO WHILE(ABS(dT)> tolerance)

    !	Calculate sound speed assuming fixed Q

           cs = Qfrag*pi*G*sigma/omega

           !	Calculate scale height

           H = cs/omega

           rhomid = sigma/(2.0*H)

           !	Use EoS to calculate tau, gamma, betac

           CALL eos_cs(rhomid, cs)

           gamma = gammamuT(1)
           mu = gammamuT(2)
           T = gammamuT(3)
           kappa = gammamuT(4)*metallicity

           tau = sigma*kappa

           ! If no irradiation, then set all irradiation parameters to zero
           IF(irrchoice==-1) THEN
              T_irr = 0.0
              Q_irr = 0.0
              cs_irr = 0.0
           ENDIF

           ! If irradiation at fixed Q_irr, calculate other irradiation parameters
           IF(irrchoice==0) THEN
              cs_irr = Q_irr*pi*G*sigma/omega
              T_irr = 0.0
              IF(cs_irr/=0.0) THEN
                 CALL eos_cs(rhomid,cs_irr)                 
                 T_irr = gammamuT(3)
              ENDIF

              ! If irradiation uses T_irr, calculate this and other parameters
           ELSE IF(irrchoice>0) THEN
              ! Ida and Lin prescription
              ! Modify T_irr
              IF(irrchoice==2) THEN
                 rAU = r/udist
                 Msol = Mstar/umass
                 T_irr = 280.0*Msol/sqrt(rAU)
                 ! Must account for optical depth
                 IF(tau/=0.0)T_irr = T_irr/(tau+1.0/tau)**0.25
              ENDIF

              cs_irr = SQRT(gamma*Boltzmannk*T_irr/(mu*mH))
              Q_irr = cs_irr*omega/(pi*G*sigma)
              ! If Q_irr > crit then fix it to crit
              IF(Q_irr>Q_irrcrit) THEN

                 Q_irr = Q_irrcrit
                 cs_irr = Q_irr*pi*G*sigma/omega
                 T_irr = 0.0
                 IF(cs_irr/=0.0) THEN
                    CALL eos_cs(rhomid,cs_irr)                 
                    T_irr = gammamuT(3)
                 ENDIF
              ENDIF
           ENDIF

           ! Calculate cooling timescale for these parameters

           betac = (tau+1.0/tau)*(cs*cs)*omega/&
                (sigma_SB*(T**4.0-T_irr**4.0)*gamma*(gamma-1.0))

           !	Calculate alpha from this value --> accretion rate

           alpha = 4.0/(9.0*gamma*(gamma-1.0)*betac)

           !	Compare calculated mdot with imposed mdot

           mdot_try = 3.0*pi*alpha*cs*cs*sigma/omega
           dT = (mdotvisc-mdot_try)/mdotvisc

           IF(ntries> 1000) THEN 
              fine = fine/10.0
              ntries = 0.0
           ENDIF

           dsig = (sigma-sigma_old)/sigma_old         

           IF(ABS(dsig)<1.0e-5.and.ntries>500) exit  ! Exit if percentage change in sigma very small

           sigma_old = sigma
           sigma = sigma*(1.0 +dT/(abs(dT))*fine)

           ntries = ntries+1
	ENDDO


 ! Check for MRI activation
 ! If so, then set alpha=0.01 and readjust sigma to maintain mdotvisc

        IF(T>1000.0) THEN
           !	print*, 'MRI active ',r/udist, T, alpha, sigma
           alpha = 0.01
           sigma = mdotvisc*omega/(3.0*pi*alpha*cs*cs)
           !	print*, 'MRI active ', T, alpha, sigma
        ENDIF

        !	Calculate enclosed Mass and mass ratio

        mtot = mtot + twopi*r*sigma*dr
        qratio = mtot/Mstar

        
        !	Calculate Local Jeans Mass and jeans length

        ! First strength of density perturbations
        deltasigma = 1.0+4.47*sqrt(alpha)

        mjeans = mjeans*cs*cs*sqrt(Qfrag)*H/deltasigma
        ljeans = sqrt(2.0*pi*pi*Qfrag/deltasigma)*H
        ljeans = ljeans/udist

        !	Calculate fragment Hill Radius

        rhill = r*(mjeans/Mstar)**0.333
        rhill = rhill/udist

        ! Convert Jeans mass to Jupiter masses
        mjeans = mjeans/mjup

        ! Calculate if fragmentation is likely: requires two conditions to be met:
        ! 1. gamma_J small and negative (Jeans mass decreases rapidly)
        ! 2. gamma_Q large (Q relatively unchanging)

        ! First, calculate gamma_Q from gamma_sigma, gamma_omega

        IF(alpha<alpha_sat) THEN
           gamma_J = 1.25*(9.0*alpha*gamma*(gamma-1)/4.0 - 1/betac) - 1.5/gamma_sigma + 0.5/gamma_omega
           gamma_Q = 0.5*(9.0*alpha*gamma*(gamma-1)/4.0 - 1/betac) - 1.0/gamma_sigma + 1.0/gamma_omega
        ELSE
           gamma_J = 1.5*(9.0*alpha_sat*gamma*(gamma-1)/4.0 - 1/betac) - 1.5/gamma_sigma - 0.5/gamma_omega
           gamma_Q = 1.25*(9.0*alpha*gamma*(gamma-1)/4.0 - 1/betac) - 1.5/gamma_sigma -0.5/gamma_omega
        ENDIF

        gamma_J = 1.0/gamma_J
        gamma_Q = 1.0/gamma_Q

        ! Test for fragmentation

        frag = 1  ! Set fragmentation flag to 'true' initially
        selfgrav = 1 ! Set self-gravitating flag to 'true' initially

        ! If gamma_Q not large enough, disc not self-gravitating

        IF(abs(gamma_Q) < gamma_Qcrit) selfgrav = 0

        ! If gamma_J not small and negative, disc non-fragmenting
        IF(gamma_J < gamma_Jcrit .or. gamma_J > 0.0) frag = 0

        ! If irradiation completely suppresses self-gravity, not self-gravitating
        IF(Q_irr==Q_irrcrit) selfgrav = 0

        ! If disc will be completely accreted within 5 orbital periods, not of interest
        tevolve = mtot*omega/mdotvisc
        if(tevolve < 5) selfgrav = 0

        ! Now check fragmenting disc is self-gravitating
        frag = frag*selfgrav


        IF(frag==0) THEN
           mjeans = 0.0
           rhill=0.0
           ljeans = 0.0
        ENDIF

        IF(selfgrav==0) THEN
           sigma = 0.0
           cs = 0.0
           T = 0.0
           alpha = 0.0
           betac=0.0
           T_irr = 0.0
           cs_irr = 0.0
           Q_irr = 0.0
        ENDIF

        !	Write disc model to file

        WRITE(10,*) r/udist, mdotvisc*3.15e7/umass, qratio,sigma,cs,omega, T,betac,alpha,mjeans, &
             ljeans, rhill,H, T_irr,cs_irr,Q_irr


     ENDDO
  ENDDO


END PROGRAM selfgravdisc_modelgrid
