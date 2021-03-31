!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University;
   ! Colorado State University Research Foundation ; ATMET, LLC

   ! This software is free software; you can redistribute it and/or modify it
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version.

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.

   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA
   ! (http://www.gnu.org/licenses/gpl.html)
   !----------------------------------------------------------------------------

!===============================================================================

Module mem_average_vars

  implicit none

  ! Number of time points to average over in a day

  integer :: npoints_davg

  ! Number of levels on which to average the 3-d variables

  integer, parameter :: nz_avg = 7

  ! Daily values - atmosphere

  real, allocatable ::  press_davg(:)
  real, allocatable ::    vxe_davg(:)
  real, allocatable ::    vye_davg(:)
  real, allocatable ::    vze_davg(:)
  real, allocatable :: rshort_davg(:)

  real, allocatable :: tempk_davg(:)
  real, allocatable :: tempk_dmin(:)
  real, allocatable :: tempk_dmax(:)

  real, allocatable :: accpmic_dtot(:)
  real, allocatable :: accpcon_dtot(:)

  real, allocatable :: press_ul_davg(:)
  real, allocatable ::   vxe_ul_davg(:)
  real, allocatable ::   vye_ul_davg(:)
  real, allocatable ::   vze_ul_davg(:)

  ! Daily values - surface grid

  real, allocatable :: airtempk_davg(:)
  real, allocatable :: airtempk_dmin(:)
  real, allocatable :: airtempk_dmax(:)

  real, allocatable :: cantempk_davg(:)
  real, allocatable :: cantempk_dmin(:)
  real, allocatable :: cantempk_dmax(:)

  real, allocatable :: sfluxt_davg(:)
  real, allocatable :: sfluxr_davg(:)

  ! Daily values - land

  real, allocatable ::  vegtempk_davg(:)
  real, allocatable ::  vegtempk_dmin(:)
  real, allocatable ::  vegtempk_dmax(:)

  real, allocatable :: soiltempk_davg(:)
  real, allocatable :: soiltempk_dmin(:)
  real, allocatable :: soiltempk_dmax(:)

Contains

!===============================================================================

  subroutine alloc_average_vars(mwa, mwsfc, mland)

    use oname_coms, only: nl
    implicit none

    integer, intent(in) :: mwa, mwsfc, mland

  ! Daily values - atmosphere

    if (nl%ioutput_davg == 1) then

       allocate( press_davg(mwa))
       allocate(   vxe_davg(mwa))
       allocate(   vye_davg(mwa))
       allocate(   vze_davg(mwa))
       allocate(rshort_davg(mwa))

       allocate(tempk_davg(mwa))
       allocate(tempk_dmin(mwa))
       allocate(tempk_dmax(mwa))

       allocate(accpmic_dtot(mwa))
       allocate(accpcon_dtot(mwa))

       allocate(press_ul_davg(mwa))
       allocate(  vxe_ul_davg(mwa))
       allocate(  vye_ul_davg(mwa))
       allocate(  vze_ul_davg(mwa))

  ! Daily values - sfc grid

       allocate(airtempk_davg(mwsfc))
       allocate(cantempk_davg(mwsfc))

       allocate(sfluxt_davg(mwsfc))
       allocate(sfluxr_davg(mwsfc))

  ! Daily values - land

       allocate(vegtempk_davg(mland))
       allocate(soiltempk_davg(mland))

    endif

    ! Always allocate the following because they are used in subroutine flux_accum

    allocate(airtempk_dmin(mwsfc))
    allocate(airtempk_dmax(mwsfc))

    allocate(cantempk_dmin(mwsfc))
    allocate(cantempk_dmax(mwsfc))

    allocate(vegtempk_dmin(mland))
    allocate(vegtempk_dmax(mland))

    allocate(soiltempk_dmin(mland))
    allocate(soiltempk_dmax(mland))

  end subroutine alloc_average_vars

!===============================================================================

  subroutine reset_davg_vars()

    use oname_coms, only: nl
    implicit none

    if (nl%ioutput_davg == 1) then

         npoints_davg = 0

        press_davg(:) = 0.0
          vxe_davg(:) = 0.0
          vye_davg(:) = 0.0
          vze_davg(:) = 0.0
       rshort_davg(:) = 0.0

       tempk_davg(:) = 0.0
       tempk_dmin(:) = 10000.0
       tempk_dmax(:) = 0.0

       accpmic_dtot(:) = 0.0
       accpcon_dtot(:) = 0.0

       press_ul_davg(:) = 0.0
         vxe_ul_davg(:) = 0.0
         vye_ul_davg(:) = 0.0
         vze_ul_davg(:) = 0.0

       ! Daily values - sfc grid

       airtempk_davg(:) = 0.0
       cantempk_davg(:) = 0.0

         sfluxt_davg(:) = 0.0
         sfluxr_davg(:) = 0.0

       ! Daily values - land

         vegtempk_davg(:) = 0.0
        soiltempk_davg(:) = 0.0

    endif

    ! Always reset the following because they are used in subroutine flux_accum

     airtempk_dmin(:) = 10000.0
     airtempk_dmax(:) = 0.0

     cantempk_dmin(:) = 10000.0
     cantempk_dmax(:) = 0.0

       vegtempk_dmin(:) = 10000.0
       vegtempk_dmax(:) = 0.0

      soiltempk_dmin(:) = 10000.0
      soiltempk_dmax(:) = 0.0

  end subroutine reset_davg_vars

End Module mem_average_vars
