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
      real, parameter :: pi = 3.14159265285
      real, parameter :: twopi = 2.0*pi
      real, parameter :: roottwopi = sqrt(twopi)

      real, parameter :: c = 3.0e10

      real(kind=8), parameter :: Boltzmannk = 1.3807d-16
      real(kind=8), parameter :: mH = 1.6726d-24
      real(kind=8), parameter :: G = 6.672041d-8
      real(kind=8), parameter :: sigma_SB = 5.67e-5
      real(kind=8), parameter :: coll_H = 2e-15 ! Collisional cross section of H2
      real(kind=8), parameter :: udist = 1.50d13
      real(kind=8), parameter :: umass = 1.99d33
      real(kind=8), parameter :: mjup = 1.8986e30
      real(kind=8), parameter :: mearth = 5.972e27
      real(kind=8), parameter :: year = 3.15e7
      real(kind=8), parameter :: pc = 206265.0*1.496e13  ! Parsec in cm
      real(kind=8), parameter :: utime = sqrt((udist**3)/(G*umass))
      real(kind=8), parameter :: uergg = udist*udist/(utime*utime)

! Arrays
      real,allocatable,dimension(:) :: gammamuT
      real,allocatable,dimension(:,:,:) :: eostable
      real,allocatable,dimension(:,:) :: cstab

      end module eosdata
