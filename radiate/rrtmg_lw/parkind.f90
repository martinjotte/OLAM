      module parkind

      implicit none

!------------------------------------------------------------------
! rrtmg kinds
! Define integer and real kinds for various types.
!
! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!     integer kinds
!     -------------
!
!     integer, parameter :: kind_ib = selected_int_kind(13)  ! 8 byte integer
!     integer, parameter :: kind_im = selected_int_kind(6)   ! 4 byte integer
      integer, parameter :: kind_in = kind(1)                ! native integer
      integer, parameter :: kind_im = kind(1)                ! native integer
      integer, parameter :: kind_ib = kind(1)                ! native integer

!
!     real kinds
!     ----------
!
!     integer, parameter :: kind_rb = selected_real_kind(12) ! 8 byte real
!     integer, parameter :: kind_rm = selected_real_kind(6)  ! 4 byte real
      integer, parameter :: kind_rn = kind(1.0)              ! native real
      integer, parameter :: kind_rm = kind(1.0)              ! native real
      integer, parameter :: kind_rb = kind(1.0)              ! native real

      ! min cloud fraction below which clouds are not considered
      real, parameter :: cldmin = 0.025

      ! max cloud fraction above which layer is considers completely clouds
      real, parameter :: cldmax = 0.975

      end module parkind
