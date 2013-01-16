      subroutine eos_T(rho,T)
!-----------------------------------------------------------------------
! Reads and interpolates the equation of state tables to give
! temperatures, opacities etc
! Updated to fortran 90 PJC 22/05/2008
! Extended to opacities DHF 01/09/2009
!-----------------------------------------------------------------------
      
      use eosdata

      implicit none

      integer :: i,j,p
      real :: mT,mT1,mT2,cT,cT1,cT2,T1,T2
      real :: mmu,mmu1,mmu2,cmu,cmu1,cmu2,mu1,mu2
      real :: mkap,mkap1,mkap2,ckap,ckap1,ckap2,kap1,kap2 
      real :: mkbar,mkbar1,mkbar2,ckbar,ckbar1,ckbar2,kbar1,kbar2
      real :: cv, rho, T, gamma
      	  		      
! gammamuT columns:
! 1) gamma
! 2) mu
! 3) T
! 4) Opacity
! 5) Mass Weighted Opacity (Stamatellos et al 2007)
	gammamuT = 0.0
		      
! Find the relevant records in the table...
! ... for rho
         if (rho < 1.0e-24) rho = 1.0e-24
         i = 1
         do 
            if ((eostable(i,1,1) >= rho).or.(i == nrhopoints)) exit
            i = i + 1
         enddo
			
! ... and for internal energy
         if (T < eostable(i,1,2)) T = eostable(i,1,2)
         j = 1
         do 
            if ((eostable(i-1,j,2) >= T).or.(j==nUpoints)) exit
            j = j + 1
         enddo

	
! Interpolate over the j value at i-1
         mT1 = (eostable(i-1,j-1,3) - eostable(i-1,j,3))/ &
              (eostable(i-1,j-1,2) - eostable(i-1,j,2))
         cT1 = eostable(i-1,j,3) - mT1*eostable(i-1,j,2)
         T1 = mT1*T + cT1

         mmu1 = (eostable(i-1,j-1,4) - eostable(i-1,j,4))/ &
              (eostable(i-1,j-1,2) - eostable(i-1,j,2))
         cmu1 = eostable(i-1,j,4) - mmu1*eostable(i-1,j,2)
         mu1 = mmu1*T + cmu1
		 
	 mkap1 = (eostable(i-1,j-1,6) - eostable(i-1,j,6))/ &
              (eostable(i-1,j-1,2) - eostable(i-1,j,2))
	 ckap1 = eostable(i-1,j,6) - mkap1*eostable(i-1,j,2)
	 kap1 = mkap1*T + ckap1
		 
	 mkbar1 = (eostable(i-1,j-1,5) - eostable(i-1,j,5))/ &
              (eostable(i-1,j-1,2) - eostable(i-1,j,2))
	ckbar1 = eostable(i-1,j,5) - mkbar1*eostable(i-1,j,2)
	kbar1 = mkbar1*T + ckbar1
				 
! Then interpolate over the j value at i
! Update j value as necessary
         j = 1
         do 
            if ((eostable(i,j,2) >= T).or.(j==nUpoints)) exit
            j = j + 1
         enddo

         mT2 = (eostable(i,j-1,3) - eostable(i,j,3))/ &
              (eostable(i,j-1,2) - eostable(i,j,2))
         cT2 = eostable(i,j,3) - mT2*eostable(i,j,2)
         T2 = mT2*T + cT2
         
         mmu2 = (eostable(i,j-1,4) - eostable(i,j,4))/ &
              (eostable(i,j-1,2) - eostable(i,j,2))
         cmu2 = eostable(i,j,4) - mmu2*eostable(i,j,2)
         mu2 = mmu2*T + cmu2
		 
         mkap2 = (eostable(i,j-1,6) - eostable(i,j,6))/ &
              (eostable(i,j-1,2) - eostable(i,j,2))
         ckap2 = eostable(i,j,6) - mkap2*eostable(i,j,2)
         kap2 = mkap2*T + ckap2
		 
	 mkbar2 = (eostable(i,j-1,5) - eostable(i,j,5))/ &
              (eostable(i,j-1,2) - eostable(i,j,2))
         ckbar2 = eostable(i,j,5) - mkbar2*eostable(i,j,2)
         kbar2 = mkbar2*T + ckbar2

! Finally interpolate over i at the fractional j value
         mT = (T2 - T1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cT = T2 - mT*eostable(i,1,1)
         gammamuT(3) = mT*rho + cT

         mmu = (mu2 - mu1) / (eostable(i,1,1)-eostable(i-1,1,1))
         cmu = mu2 - mmu*eostable(i,1,1)
         gammamuT(2) = mmu*rho + cmu
		 
	mkap = (kap2 - kap1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ckap = kap2 - mkap*eostable(i,1,1)
         gammamuT(4) = mkap*rho + ckap
		 
	mkbar = (kbar2 - kbar1) / (eostable(i,1,1)-eostable(i-1,1,1))
         ckbar = kbar2 - mkbar*eostable(i,1,1)
         gammamuT(5) = mkbar*rho + ckbar

! Evaluate the gamma value from gamma = 1 + k/(mH*mu*cv), with cv = U/T
           cv = gammamuT(3)/T
         gammamuT(1) = 1.0d0 + (Boltzmannk / &
              (mH*gammamuT(2)*cv))

             gammamuT(3) = T

      end subroutine eos_T
