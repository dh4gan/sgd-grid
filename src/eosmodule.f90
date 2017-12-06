      module eosdata

!-----------------------------------------------------------------------
! Data module for saving equation of state values
! PJC 22/05/2008
! DHF 09/09/2009
!-----------------------------------------------------------------------

      implicit none
      save

! Units
      integer :: nrhopoints,nUpoints
      real(kind=8),parameter :: Boltzmannk = 1.3807d-16
      real(kind=8),parameter :: mH = 1.6726d-24
      real(kind=8),parameter :: G = 6.672041d-8
      real(kind=8),parameter :: udist = 1.50d13
      real(kind=8),parameter :: umass = 1.99d33
      real(kind=8),parameter :: utime = sqrt((udist**3)/(G*umass))
      real(kind=8),parameter :: uergg = udist*udist/(utime*utime)

! Arrays
     real,allocatable,dimension(:) :: gammamuT
      real,allocatable,dimension(:,:,:) :: eostable
	real,allocatable,dimension(:,:) :: cstab

      end module eosdata
