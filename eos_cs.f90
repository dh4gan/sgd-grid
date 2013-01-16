      subroutine eos_cs(rho,cs)
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
      real :: cv, rho, cs, gamma
      	  		      
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
         if (cs < cstab(1,1)) cs = cstab(1,1)
         j = 1
         do 
            if ((cstab(i-1,j) >= cs).or.(j==nUpoints)) exit
            j = j + 1
         enddo

	IF(j==1) j = j+1
			
! Interpolate over the j value at i-1
         mT1 = (eostable(i-1,j-1,2) - eostable(i-1,j,2))/ &
              (cstab(i-1,j-1) - cstab(i-1,j))
         cT1 = eostable(i-1,j,2) - mT1*cstab(i-1,j)
         T1 = mT1*cs + cT1

         mmu1 = (eostable(i-1,j-1,4) - eostable(i-1,j,4))/ &
              (cstab(i-1,j-1) - cstab(i-1,j))
         cmu1 = eostable(i-1,j,4) - mmu1*cstab(i-1,j)
         mu1 = mmu1*cs + cmu1
		 
	 mkap1 = (eostable(i-1,j-1,6) - eostable(i-1,j,6))/ &
              (cstab(i-1,j-1) - cstab(i-1,j))
	 ckap1 = eostable(i-1,j,6) - mkap1*cstab(i-1,j)
	 kap1 = mkap1*cs + ckap1
		 
	 mkbar1 = (eostable(i-1,j-1,5) - eostable(i-1,j,5))/ &
              (cstab(i-1,j-1) - cstab(i-1,j))
	ckbar1 = eostable(i-1,j,5) - mkbar1*cstab(i-1,j)
	kbar1 = mkbar1*cs + ckbar1
				 
! Then interpolate over the j value at i
! Update j value as necessary
         j = 1
         do 
            if ((cstab(i,j) >= cs).or.(j==nUpoints)) exit
            j = j + 1
         enddo

	IF(j==1) j = j+1
			
         mT2 = (eostable(i,j-1,2) - eostable(i,j,2))/ &
              (cstab(i,j-1) - cstab(i,j))
         cT2 = eostable(i,j,2) - mT2*cstab(i,j)
         T2 = mT2*cs + cT2 
         
         mmu2 = (eostable(i,j-1,4) - eostable(i,j,4))/ &
              (cstab(i,j-1) - cstab(i,j))
         cmu2 = eostable(i,j,4) - mmu2*cstab(i,j)
         mu2 = mmu2*cs + cmu2
		 
         mkap2 = (eostable(i,j-1,6) - eostable(i,j,6))/ &
              (cstab(i,j-1) - cstab(i,j))
         ckap2 = eostable(i,j,6) - mkap2*cstab(i,j)
         kap2 = mkap2*cs + ckap2
		 
	 mkbar2 = (eostable(i,j-1,5) - eostable(i,j,5))/ &
              (cstab(i,j-1) - cstab(i,j))
         ckbar2 = eostable(i,j,5) - mkbar2*cstab(i,j)
         kbar2 = mkbar2*cs + ckbar2

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
         
         gammamuT(1) = cs*cs*gammamuT(2)*mH/(Boltzmannk*gammamuT(3))
              

      end subroutine eos_cs
