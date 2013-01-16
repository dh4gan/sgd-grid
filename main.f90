PROGRAM survey_jeans
  !	Program does a parameter survey of self-gravitating disc models
  !	Writes disc properties and fragmentation outcomes...
  !	********Assumes fixed Q***********
  !	****** *Surveys the mdot - r parameter space only****
  !	Finds local Jeans mass ==> potential fragment masses

  use eosdata
  implicit none

  integer :: rflag, frag,i,iparam1,iparam2,ichoose1,ichoose2,ifix,k,irrchoice

  integer, parameter :: nsurvey = 1000

  real, parameter :: pi = 3.14159265285
  real, parameter :: twopi = 2.0*pi
  real, parameter :: sigma_SB = 5.67e-5
  real,parameter :: c = 3.0e10
  real, parameter :: tolerance = 1.0e-5
  real, parameter :: Qfrag = 2.0
  real, parameter :: mjup = 1.8986e30
  real,parameter :: Q_irrcrit = 1.99
  real, parameter :: distance = 140.0*206265.0*1.496e13  ! Distance to object in cm
  real,parameter :: kappa0 = 0.035 ! Frequency dependent opacity constant
  real,parameter :: nu0 = c/(8.5e-2) ! Frequency constant
  real,parameter :: beta_k = 0.6 ! Power law index for frequency dependent opacity
real,parameter :: Tdust = 30.0 ! Assumed dust temperature for inferring disc mass from flux
real,parameter :: metallicity = 1.0 ! Metallicity relative to solar
  !  real,parameter :: Q_irrcrit = 2.0
  real :: q, Mstar, rmax,rmin,dr,r, qratio,kappa
  real :: mdotmax, mdotmin,sigma_old, p_sig, dmdot, mu
  real :: mdotvisc, Ttry,dT,ntries,fixvalue, r_c,mtot,gamma
  real :: dro,dq,romin,romax,qmin,qmax,jeansmin,mdot_try,fine

  real :: fluxtot1,fluxtot2, nu1,nu2,kappa1,kappa2, flux_nu1,flux_nu2
  real :: mflux1,mflux2
  real :: sigma, omega, T, cs,rhomid,alpha
  real :: mdot_eq, alpha_eq,betacrit
  real :: Tfrag,tau,H,betac,mjeans,njeans,dsig
  real :: betasigma, beta_J,rhill,rstall,l_jeans,R_HP
  real :: cs_irr, T_irr, Q_irr, Msol, rAU
  real,dimension(3) :: param, maxparam, minparam, dparam
  character(100) :: outputfile
  character(5) :: char1,char2

  ! Get input parameters and set up header:
  print*, " "
  print*, "-----------------------------------------------"  
  print*, "           JEANS MASS IN SIMPLE DISC MODELS    " 
  print*, "     Created by D. Forgan, May 2011       "
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
  READ(10,*) Mstar
  READ(10,*) rmin
  READ(10,*) rmax
  read(10,*) mdotmin
  read(10,*) mdotmax
  read(10,*) betasigma
  read(10,*) irrchoice
  read(10,*) Q_irr
  read(10,*) T_irr
  close(10)

  print*, 'Outputting to file ', outputfile

  rmin = rmin*udist
  rmax = rmax*udist


  dr = (rmax-rmin)/REAL(nsurvey)			
  dmdot = (mdotmax-mdotmin)/REAL(nsurvey)							
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
  write(10,*) nsurvey*nsurvey		

  ! Loop over parameters

  DO iparam1 = 1, nsurvey

     mdotvisc = mdotmin + (iparam1-1)*dmdot				

     print*, 'Testing mdot ', iparam1, 10.0**mdotvisc

     mdotvisc = (10.0**mdotvisc)*umass/3.15e7

     mtot = 0.0
     fluxtot1 = 0.0
     fluxtot2 = 0.0

     ! Loop over radii
     DO iparam2 = 1,nsurvey

	r = rmin +(iparam2-1)*dr 
        sigma = 0.0
        omega = 0.0
        T = 0.0
        cs = 0.0
        Tfrag = 0.0
        mjeans = 0.0
        njeans = 0.0
	betac = 0.0
	alpha = 0.0
	alpha_eq = 0.0
	mdot_eq = 0.0				

 ! begin with simple factors first

        mjeans = 4.0*1.412*pi*pi*pi/(3.0*G)

        ! print*, 'Calculating for ',param(ichoose1),param(ichoose2),param(ifix)

	frag =0
	jeansmin = 1.0e30	

 ! Calculate Omega

	omega = sqrt(G*Mstar/r**3)

 ! Begin iteration process

 ! Guess Sigma (start from previous answer)

	IF(iparam2==1) sigma_old = 50000.0

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
           IF(irrchoice==-1) THEN
              T_irr = 0.0
              Q_irr = 0.0
              cs_irr = 0.0
           ENDIF

           IF(irrchoice==0) THEN
              cs_irr = Q_irr*pi*G*sigma/omega
              T_irr = 0.0
              IF(cs_irr/=0.0) THEN
                 CALL eos_cs(rhomid,cs_irr)                 
                 T_irr = gammamuT(3)
              ENDIF
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


        ! Calculate flux generated in this annulus

        ! First, 850 micron flux

        nu1 = c/(8.7e-2)
        kappa1 = kappa0*(nu1/nu0)**(beta_k)

        if (sigma*kappa1 < 1.0) THEN
           flux_nu1 = (2.0*Boltzmannk/(c*c))*nu1*nu1*kappa1*sigma*T*twopi*r*dr/(distance*distance)
        else
           flux_nu1 = (2.0*Boltzmannk/(c*c))*nu1*nu1*(T/(sigma*kappa1)**0.25)*twopi*r*dr/(distance*distance)

        ENDIF

        ! Now 3.4 mm flux
        nu2 = c/(0.34)
        kappa2 = kappa0*(nu2/nu0)**(beta_k)

        if (sigma*kappa2 < 1.0) THEN
           flux_nu2 = (2.0*Boltzmannk/(c*c))*nu2*nu2*kappa2*sigma*T*twopi*r*dr/(distance*distance)
        else
           flux_nu2 = (2.0*Boltzmannk/(c*c))*nu2*nu2*(T/(sigma*kappa2)**0.25)*twopi*r*dr/(distance*distance)

        ENDIF


        !	Calculate enclosed Mass and mass ratio

        mtot = mtot + twopi*r*sigma*dr
        fluxtot1 = fluxtot1 + flux_nu1
        fluxtot2 = fluxtot2 + flux_nu2
        qratio = mtot/Mstar

        ! Calculate predicted mass from two different fluxes (given estimated dust temperature)

        mflux1 = distance*distance*fluxtot1*c*c/(2.0*kappa1*nu1*nu1*Boltzmannk*Tdust)
        mflux2 = distance*distance*fluxtot2*c*c/(2.0*kappa2*nu2*nu2*Boltzmannk*Tdust)

        !	Calculate Local Jeans Mass and jeans length
        !      In irradiated discs case, delta sigma/sigma taken from Rice et al (2011)
        mjeans = mjeans*cs*cs*sqrt(Qfrag)*H/(1.0+4.47*sqrt(alpha))
        l_jeans = (1.0/(1.0+1.0/sqrt(betac)))**(0.6666)*sqrt(2.0)*pi*sqrt(Qfrag)*H		
        njeans = 2.0*pi*H*H*sigma*(1.0+1.0/sqrt(betac))/(3.0*mjeans)			

        !	Calculate fragment survival radii

        rstall = (gamma*G*mjeans)/(cs*cs)
        rhill = r*(mjeans/Mstar)**0.333

        R_HP = ((3.0*gamma-4.0)/(3.0*gamma-1.0))**0.333*rhill

        !	Convert lengths to AU, and Masses to Jupiter Masses
        rstall = rstall/udist
        rhill = rhill/udist
        l_jeans = l_jeans/udist
        R_HP = R_HP/udist

        mjeans = mjeans/mjup		! convert to Jupiter masses

        !	Calculate if fragmentation is likely (is beta_J small and negative?)

        IF(alpha<0.1) THEN
           beta_J = 1.5*(9.0*alpha*gamma*(gamma-1)/4.0 - 1/betac) - 1.0/betasigma
        ELSE
           beta_J = 1.5*(9.0*0.1*gamma*(gamma-1)/4.0 - 1/betac) - 1.0/betasigma
        ENDIF

        beta_J = 1.0/beta_J

        IF(Q_irr/=Q_irrcrit.and.beta_J<0.0.and.ABS(beta_J)<3.and.frag==0) THEN
           frag=1
           !	print*, 'Fragmentation in progress'
           !	print*, mdotvisc*3.15e7/umass, r/udist,qratio, alpha, betac,beta_J
        ENDIF

        !	Write disc profiles to separate file

        !	betacrit = sqrt(2.0*pi)/(3.0*gamma-4.0)

        !	IF(betac<betacrit) THEN
        !	rstall = 0.0
        !	r_HP = 0.0
        !	ENDIF

        IF(frag==0) THEN
           mjeans = 0.0
           rhill=0.0
           rstall = 0.0
           R_HP = 0.0
           l_jeans = 0.0
        ENDIF
        !	print*, 'Writing to file: ',iparam1,iparam2, r/udist, mdotvisc*3.15e7/umass		
        !	Write fragmentation stats to separate file

        WRITE(10,*) iparam1,iparam2, r/udist, mdotvisc*3.15e7/umass, qratio,frag,sigma,betac,cs,T,alpha,mjeans, &
             l_jeans, rstall, rhill, R_HP, T_irr,cs_irr,Q_irr, flux_nu1/1e-23,fluxtot1/1e-23,mflux1/umass, &
             flux_nu2/1e-23,fluxtot2/1e-23,mflux2/umass

       ! print*, iparam1,iparam2, r/udist, mdotvisc*3.15e7/umass, T,tau,kappa1,nu1,sigma,flux_nu1/1e-23,fluxtot1/1e-23,mflux1/umass

     ENDDO
  ENDDO


END PROGRAM survey_jeans
