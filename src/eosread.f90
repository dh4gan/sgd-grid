      subroutine eosread
!-----------------------------------------------------------------------
! Reads the equation of state tables
! Updated to fortran 90 PJC 22/05/2008
! Extended to opacities DHF 01/09/2009
!
! Note the columns (== k values) are as follows:
!   1) density
!   2) temperature
!   3) internal energy
!   4) mean molecular weight
!	5) Mass Weighted Opacity
!	6) Mass-weighted Opacity
!
!-----------------------------------------------------------------------

      use eosdata

      implicit none
      integer :: i,j,k,check
	real :: cv, gamma
! Open eos data file and read in values

      open(50,file="/Users/dhf/Programs/selfgravdisc_modelgrid/params/myeos.dat", &
           status='old',iostat=check,action='read')
      if (check /= 0) then
         print*, 'Input file myeos.dat not found'
         stop
      endif
      print*, "Reading equation of state tables"
      
      read(50,*) nrhopoints, nUpoints

! Allocate input data array
      allocate(eostable(nrhopoints,nUpoints,6))
	allocate(cstab(nrhopoints,nUpoints))
      do i = 1,nrhopoints
         do j = 1,nUpoints
            read(50,*) (eostable(i,j,k),k=1,6)
	  
!	Now calculate sound speed values for interpolation
	cv = eostable(i,j,3)/eostable(i,j,2)
	gamma = 1.0d0 + (Boltzmannk / &
              (mH*eostable(i,j,4)*cv))
	cstab(i,j) = sqrt(gamma*(gamma-1.0)*eostable(i,j,3))			    
         enddo
      enddo
      close(50)
	
      print*, "Equation of state tables read in successfully"
      print*, "-----------------------------------------------"
      
      return
      end subroutine eosread
        
