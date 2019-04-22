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
subroutine inithh()

  ! Horizontally-homogeneous initialization of model fields

  use misc_coms, only: itsflg
  implicit none

  ! Arrange the input sounding

  if (itsflg == 3) then
     call arrsnd_thil()
  elseif (itsflg == 2) then
     call arrsnd_theta()
  else
     call arrsnd_tair()
  endif

  ! Compute the 1-D reference state variables
  call refs1d()

  ! Initialize the 3-D model fields
  call fldshhi()

end subroutine inithh

!===============================================================================

subroutine arrsnd_thil()

use consts_coms, only: pio180, eps_virt, p00k, rocp, cpor, p00i, &
                       t00, eps_vap, rdryog, alvlocp, gocp
use misc_coms,   only: nsndg, ipsflg, irtsflg, iusflg, ps, ts, rts, us, vs, &
                       hs, thds, p_sfc
use therm_lib,   only: eslf

implicit none

integer :: ksndg, iterate
real    :: dir, spd, thavg, esat
real    :: tvirt1, tvirt2, press, theta, tair, qt, qv, qc, qsat
real    :: qvap(nsndg)

real, external :: rts2qt

! Convert winds from spd/direction to U and V components

if (iusflg /= 0) then
   do ksndg = 1, nsndg
      dir = us(ksndg)
      spd = vs(ksndg)
      us(ksndg) = -spd * sin(pio180*dir)
      vs(ksndg) = -spd * cos(pio180*dir)
   enddo
endif

if (ipsflg == 0) then

   ! PS array is pressure in mb with HS(1) = height of lowest sounding level

   do ksndg = 1, nsndg
      ps(ksndg) = ps(ksndg) * 100.
   enddo

   ! Compute temperature and humidity. Need to iterate in case of saturation
   do ksndg = 1, nsndg

      theta = ts(ksndg)
      tair  = theta * (ps(ksndg) * p00i) ** rocp
      qt    = rts2qt(irtsflg, rts(ksndg), tair, ps(ksndg))
      qv    = qt
      qc    = 0.0

      do iterate = 1, 20
         theta = 0.7 * theta &
               + 0.3 * ts(ksndg) * (1. + alvlocp * qc / max(tair,253.))

         tair  = theta * (ps(ksndg) * p00i) ** rocp

         qt    = rts2qt(irtsflg, rts(ksndg), tair, ps(ksndg))
         esat  = eslf(tair - t00)
         qsat  = eps_vap * esat / (ps(ksndg) - esat * (1. - eps_vap))

         qc    = max(0., qt - qsat)
         qv    = qt - qc
      enddo

      rts (ksndg) = qt
      ts  (ksndg) = tair
      thds(ksndg) = theta
      qvap(ksndg) = qv
   enddo

   do ksndg = 2, nsndg
      hs(ksndg) = hs(ksndg-1) - rdryog * 0.5  &
           * ( ts(ksndg)   * (1. + eps_virt * qvap(ksndg  ))   &
             + ts(ksndg-1) * (1. + eps_virt * qvap(ksndg-1)))  &
           * log( ps(ksndg) / ps(ksndg-1) )
   enddo

elseif (ipsflg == 1) then

   ! PS array is height in meters with P_SFC = surface pressure

   hs(1:nsndg) = ps(1:nsndg)
   ps(1)       = p_sfc * 100.

   ! We know pressure at the surface; compute temperature and humidity

   ksndg = 1

   theta = ts(ksndg)
   tair  = theta * (ps(ksndg) * p00i) ** rocp

   qt    = rts2qt(irtsflg, rts(ksndg), tair, ps(ksndg))
   qv    = qt
   qc    = 0.0

   do iterate = 1, 20
      theta = 0.7 * theta &
            + 0.3 * ts(ksndg) * (1. + alvlocp * qc / max(tair,253.))

      tair  = theta * (ps(ksndg) * p00i) ** rocp

      qt    = rts2qt(irtsflg, rts(ksndg), tair, ps(ksndg))
      esat  = eslf(tair - t00)
      qsat  = eps_vap * esat / (ps(ksndg) - esat * (1. - eps_vap))

      qc    = max(0., qt - qsat)
      qv    = qt - qc
   enddo

   rts (ksndg) = qt
   ts  (ksndg) = tair
   thds(ksndg) = theta
   qvap(ksndg) = qv

   ! Above surface, we need to compute temperature, humidity, and pressure.
   ! Need to iterate in case of saturation

   do ksndg = 2, nsndg

      ! virtual temp correction from layer below
      tvirt2 = 1. + eps_virt * qvap(ksndg-1)

      theta = ts(ksndg)
      thavg = 0.5 * (theta + thds(ksndg-1) * tvirt2)
      press = ( ps(ksndg-1) ** rocp &
              + gocp * (hs(ksndg-1) - hs(ksndg)) * p00k / thavg ) ** cpor
      tair  = theta * (press * p00i) ** rocp
      qt    = rts2qt(irtsflg, rts(ksndg), tair, press)
      qv    = qt
      qc    = 0.0

      do iterate = 1, 20

         tvirt1 = 1. + eps_virt * qv
         thavg  = 0.5 * (theta * tvirt1 + thds(ksndg-1) * tvirt2)

         press  = ( ps(ksndg-1) ** rocp &
                  + gocp * (hs(ksndg-1) - hs(ksndg)) * p00k / thavg ) ** cpor

         tair   = theta * (press * p00i)**rocp

         qt     = rts2qt(irtsflg, rts(ksndg), tair, press)
         esat   = eslf(tair - t00)
         qsat   = eps_vap * esat / (press - esat * (1. - eps_vap))

         qc     = max(0., qt - qsat)
         qv     = qt - qc

         theta = 0.7 * theta &
               + 0.3 * ts(ksndg) * (1. + alvlocp * qc / max(tair,253.))

      enddo

      rts (ksndg) = qt
      ts  (ksndg) = tair
      thds(ksndg) = theta
      ps  (ksndg) = press
      qvap(ksndg) = qv
   enddo

endif

end subroutine arrsnd_thil

!===============================================================================

subroutine arrsnd_theta()

use consts_coms, only: pio180, eps_virt, p00k, rocp, cpor, p00i, &
                       t00, eps_vap, rdryog, gocp
use misc_coms,   only: nsndg, ipsflg, irtsflg, iusflg, ps, ts, rts, us, vs, &
                       hs, thds, p_sfc
use therm_lib,   only: eslf

implicit none

integer :: ksndg, iterate
real    :: dir, spd, thavg, esat
real    :: tvirt1, tvirt2, press, theta, tair, qt, qv, qc, qsat
real    :: qvap(nsndg)

real, external :: rts2qt

! Convert winds from spd/direction to U and V components

if (iusflg /= 0) then
   do ksndg = 1, nsndg
      dir = us(ksndg)
      spd = vs(ksndg)
      us(ksndg) = -spd * sin(pio180*dir)
      vs(ksndg) = -spd * cos(pio180*dir)
   enddo
endif

if (ipsflg == 0) then

   ! PS array is pressure in mb with HS(1) = height of lowest sounding level
   ! Compute temperature and humidity

   do ksndg = 1, nsndg
      ps  (ksndg) = ps(ksndg) * 100.
      thds(ksndg) = ts(ksndg)
      ts  (ksndg) = ts(ksndg) * (ps(ksndg) * p00i) ** rocp
      rts (ksndg) = rts2qt(irtsflg, (rts(ksndg)), ts(ksndg), ps(ksndg))

      esat        = eslf(ts(ksndg) - t00)
      qsat        = eps_vap * esat / (ps(ksndg) - esat * (1. - eps_vap))
      qvap(ksndg) = rts(ksndg) - max(0., rts(ksndg) - qsat)
   enddo

   do ksndg = 2, nsndg
      hs(ksndg) = hs(ksndg-1) - rdryog * 0.5  &
           * ( ts(ksndg)   * (1. + eps_virt * qvap(ksndg  ))   &
             + ts(ksndg-1) * (1. + eps_virt * qvap(ksndg-1)))  &
           * log( ps(ksndg) / ps(ksndg-1) )
   enddo


elseif (ipsflg == 1) then

   ! PS array is height in meters with P_SFC = surface pressure

   hs(1:nsndg) = ps(1:nsndg)
   ps(1)       = p_sfc * 100.

   ! We know pressure at the surface; compute temperature and humidity

   ksndg = 1

   thds(ksndg) = ts(ksndg)
   ts  (ksndg) = ts(ksndg) * (ps(ksndg) * p00i) ** rocp
   rts (ksndg) = rts2qt(irtsflg, (rts(ksndg)), ts(ksndg), ps(ksndg))

   esat        = eslf(ts(ksndg) - t00)
   qsat        = eps_vap * esat / (ps(ksndg) - esat * (1. - eps_vap))
   qvap(ksndg) = rts(ksndg) - max(0., rts(ksndg) - qsat)

   ! Above surface, we need to compute temperature, humidity, and pressure.
   ! Need to iterate in case of saturation

   do ksndg = 2, nsndg

      ! virtual temp correction from layer below
      tvirt2 = 1. + eps_virt * qvap(ksndg-1)

      theta = ts(ksndg)
      thavg = 0.5 * (theta + thds(ksndg-1) * tvirt2)
      press = ( ps(ksndg-1) ** rocp &
              + gocp * (hs(ksndg-1) - hs(ksndg)) * p00k / thavg ) ** cpor
      tair  = theta * (press * p00i) ** rocp
      qt    = rts2qt(irtsflg, rts(ksndg), tair, press)
      qv    = qt
      qc    = 0.0

      do iterate = 1, 20

         tvirt1 = 1. + eps_virt * qv
         thavg  = 0.5 * (theta * tvirt1 + thds(ksndg-1) * tvirt2)
         press  = ( ps(ksndg-1) ** rocp &
                + gocp * (hs(ksndg-1) + hs(ksndg)) * p00k / thavg ) ** cpor

         tair   = theta * (press * p00i)**rocp

         qt     = rts2qt(irtsflg, rts(ksndg), tair, press)
         esat   = eslf(tair - t00)
         qsat   = eps_vap * esat / (press - esat * (1. - eps_vap))

         qc     = max(0., qt - qsat)
         qv     = qt - qc

      enddo

      rts (ksndg) = qt
      ts  (ksndg) = tair
      thds(ksndg) = theta
      ps  (ksndg) = press
      qvap(ksndg) = qv
   enddo

endif

end subroutine arrsnd_theta

!===============================================================================

subroutine arrsnd_tair()

use consts_coms, only: pio180, eps_virt, p00k, rocp, p00, t00, &
                       eps_vap, rdryog, gordry
use misc_coms,   only: nsndg, ipsflg, irtsflg, iusflg, ps, ts, rts, &
                       us, vs, hs, thds, p_sfc
use therm_lib,   only: eslf

implicit none

integer :: ksndg, iterate
real    :: dir, spd, tavg, esat
real    :: tvirt1, tvirt2, press, qt, qv, qc, qsat
real    :: qvap(nsndg)

real, external :: rts2qt

! Convert winds from spd/direction to U and V components

if (iusflg /= 0) then
   do ksndg = 1, nsndg
      dir = us(ksndg)
      spd = vs(ksndg)
      us(ksndg) = -spd * sin(pio180*dir)
      vs(ksndg) = -spd * cos(pio180*dir)
   enddo
endif

! Convert temperature to Kelvin

if (irtsflg == 0) then
   do ksndg = 1, nsndg
      ts(ksndg) = ts(ksndg) + t00
   enddo
endif

if (ipsflg == 0) then

   ! PS array is pressure in mb with HS(1) = height of lowest sounding level
   ! Compute potential temperature and humidity

   do ksndg = 1, nsndg
      ps  (ksndg) = ps(ksndg) * 100.
      thds(ksndg) = ts(ksndg) * (p00 / ps(ksndg)) * rocp
      rts (ksndg) = rts2qt(irtsflg, (rts(ksndg)), ts(ksndg), ps(ksndg))

      esat        = eslf(ts(ksndg) - t00)
      qsat        = eps_vap * esat / (ps(ksndg) - esat * (1. - eps_vap))
      qvap(ksndg) = rts(ksndg) - max(0., rts(ksndg) - qsat)
   enddo

   do ksndg = 2, nsndg
      hs(ksndg) = hs(ksndg-1) - rdryog * 0.5  &
           * ( ts(ksndg)   * (1. + eps_virt * qvap(ksndg  ))   &
             + ts(ksndg-1) * (1. + eps_virt * qvap(ksndg-1)))  &
           * log( ps(ksndg) / ps(ksndg-1) )
   enddo

elseif (ipsflg == 1) then

   ! PS array is height in meters with P_SFC = surface pressure

   hs(1:nsndg) = ps(1:nsndg)
   ps(1)       = p_sfc * 100.

   ! We know pressure at the surface; compute potential temperature and humidity

   ksndg = 1

   thds(ksndg) = ts(ksndg) * (p00 / ps(ksndg)) ** rocp
   rts (ksndg) = rts2qt(irtsflg, (rts(ksndg)), ts(ksndg), ps(ksndg))

   esat        = eslf(ts(ksndg) - t00)
   qsat        = eps_vap * esat / (ps(ksndg) - esat * (1. - eps_vap))
   qvap(ksndg) = rts(ksndg) - max(0., rts(ksndg) - qsat)

   ! Above surface, we need to compute potential temperature, humidity,
   ! and pressure. Need to iterate in case of saturation

   do ksndg = 2, nsndg

      ! virtual temp correction from layer below
      tvirt2 = 1. + eps_virt * qvap(ksndg-1)

      tavg  = 0.5 * (ts(ksndg) + ts(ksndg-1) * tvirt2)
      press = ps(ksndg-1) * exp(gordry * (hs(ksndg-1) - hs(ksndg)) / tavg)

      qt    = rts2qt(irtsflg, rts(ksndg), ts(ksndg), press)
      qv    = qt
      qc    = 0.0

      do iterate = 1, 20

         tvirt1 = 1. + eps_virt * qv

         tavg  = 0.5 * (ts(ksndg) * tvirt1 + ts(ksndg-1) * tvirt2)
         press = ps(ksndg-1) * exp(gordry * (hs(ksndg-1) - hs(ksndg)) / tavg)

         qt     = rts2qt(irtsflg, rts(ksndg), ts(ksndg), press)
         esat   = eslf(ts(ksndg) - t00)
         qsat   = eps_vap * esat / (press - esat * (1. - eps_vap))

         qc     = max(0., qt - qsat)
         qv     = qt - qc

      enddo

      rts (ksndg) = qt
      thds(ksndg) = ts(ksndg) * (p00 / press) ** rocp
      ps  (ksndg) = press
      qvap(ksndg) = qv

   enddo

endif

end subroutine arrsnd_tair

!===============================================================================

real function rts2qt(irtsflg, rts, tair, press)

  use therm_lib,   only: eslf
  use consts_coms, only: eps_vap, t00

  implicit none

  integer, intent(in) :: irtsflg
  real,    intent(in) :: rts
  real,    intent(in) :: tair
  real,    intent(in) :: press
  real                :: vapor_press

  if (irtsflg == 2) then

     ! RTS is given as specific humidity in g/kg
     rts2qt = rts * 1.e-3

  else

     if (irtsflg == 0) then

        ! RTS is dew point in degrees C
        vapor_press = eslf(rts)

     elseif (irtsflg == 1) then

        ! RTS is Dew point in Kelvin
        vapor_press = eslf(rts - t00)

     elseif (irtsflg == 3) then

        ! RTS is given as relative humidity in percent
        vapor_press = .01 * rts * eslf(tair - t00)

     elseif (irtsflg == 4) then

        ! RTS is given as dew point depression in Kelvin
        vapor_press = eslf(tair - rts - t00)

     elseif (irtsflg /= 2) then

        write(*,*) irtsflg
        stop 'illegal value of irtsflg'

     endif

     ! Do not allow vapor pressure to exceed ambient pressure

     vapor_press = min(press, vapor_press)

     ! Compute specific humidity from vapor pressure and ambient pressure

     rts2qt = eps_vap * vapor_press  / (press - vapor_press * (1.0 - eps_vap))

  endif

end function rts2qt

!===============================================================================

subroutine refs1d()
! +---------------------------------------------------------------------
! \   This routine computes the reference state sounding on the model
! \     levels from input sounding defined on pressure levels.
! +---------------------------------------------------------------------

use misc_coms,   only: io6, nsndg, hs, thds, us, vs, ts, rts, ps, &
                       pr01d, dn01d, rt01d, th01d, u01d, v01d
use consts_coms, only: cvocp, p00k, rdry, eps_virt, grav, gravo2, p00, &
                       rocp, t00, p00i, eps_vap
use mem_grid,    only: mza, zm, zt, dzt_top, dzt_bot
use micro_coms,  only: miclevel
use therm_lib,   only: eslf

implicit none

integer :: k, iter
real    :: dens1, exner, temp
real    :: esat, qsat, qv, tair

write(io6,*) 'Beginning refs1d '

if (zt(mza) > hs(nsndg)) then
   write(io6,*) ' !!! Input sounding is not high enough !!!'
   write(io6,*) ' !!! Sounding top (m): ',hs(nsndg)
   write(io6,*) ' !!! Model top (m):    ',zt(mza)
   stop 'refs1d'
endif

call htint(nsndg,thds,hs,mza,th01d,zt)
call htint(nsndg,  us,hs,mza, u01d,zt)
call htint(nsndg,  vs,hs,mza, v01d,zt)
if (miclevel >= 1) then
   call htint(nsndg,rts,hs,mza,rt01d,zt)
else
   rt01d(1:mza) = 0.
endif

! New interpolation of pressure prior to iterative hydrostatic integration

call htint(nsndg,ps,hs,mza,pr01d,zt)

esat  = eslf(ts(1) - t00)
qsat  = eps_vap * esat / (ps(1) - esat * (1. - eps_vap))
qv    = rts(1) - max(0., rts(1) - qsat)
dens1 = ps(1) ** cvocp * p00k / (rdry * thds(1) * (1. + eps_virt * qv))

! Use iterative method for hydrostatic integration

do iter = 1,100

   tair  = th01d(1) * (pr01d(1) * p00i) ** rocp
   esat  = eslf(tair - t00)
   qsat  = eps_vap * esat / (pr01d(1) - esat * (1. - eps_vap))
   qv    = rt01d(1) - max(0., rt01d(1) - qsat)

   dn01d(1) = pr01d(1) ** cvocp * p00k  &
            / (rdry * th01d(1) * (1. + eps_virt * qv))
   pr01d(1) = ps(1) - gravo2 * (dens1 + dn01d(1)) * (zt(1) - hs(1))

   do k = 2,mza

      tair  = th01d(k) * (pr01d(k) * p00i) ** rocp
      esat  = eslf(tair - t00)
      qsat  = eps_vap * esat / (pr01d(k) - esat * (1. - eps_vap))
      qv    = rt01d(k) - max(0., rt01d(k) - qsat)

      dn01d(k) = pr01d(k) ** cvocp * p00k  &
               / (rdry * th01d(k) * (1. + eps_virt * qv))

      ! Impose minimum value of 1 Pa to avoid overshoot to negative values
      ! during iteration
      pr01d(k) = max(1.,  &
           pr01d(k-1) - grav * (dn01d(k-1) * dzt_top(k-1) + dn01d(k) * dzt_bot(k)))
   enddo
enddo

rt01d(1) = rt01d(2)
th01d(1) = th01d(2)

! Print out initial state column

write(io6,*) ' '
write(io6,*) '============================================================================'
write(io6,*) '                    OLAM INITIAL STATE COLUMN (hhi)'
write(io6,*) '============================================================================'
write(io6,*) '   zm(m)    k     zt(m)   pr01d(Pa) dn01d(kg/m3) th01d(K)  temp  rt01d(g/kg)'
write(io6,*) '============================================================================'
write(io6,*) ' '

do k = mza,2,-1

   exner = (pr01d(k) / p00) ** rocp  ! exner WITHOUT CP factor
   temp  = th01d(k) * exner

   write(io6, '(f10.2,1x,6(''-----------''))') zm(k)
   write(io6, '(10x,i5,f10.2,f10.2,f10.4,2f10.2,f10.4)')  &
       k,zt(k),pr01d(k),dn01d(k),th01d(k),temp,rt01d(k)*1.e3
enddo

write(io6, '(f10.2,1x,9(''-------''))') zm(1)
write(io6,*) ' '

end subroutine refs1d

!===============================================================================

subroutine fldshhi()

use mem_basic,   only: theta, thil, tair, press, rho, wc, wmc, &
                       vc, vp, vmp, vmc, sh_w, sh_v
use mem_micro,   only: sh_c
use micro_coms,  only: miclevel
use mem_ijtabs,  only: jtab_w, jtab_v, itab_v, jtv_init, jtw_init, jtv_wall
use misc_coms,   only: mdomain, th01d, pr01d, dn01d, rt01d, u01d, v01d, &
                       iparallel
use consts_coms, only: cvocp, p00k, rdry, rvap, p00, p00i, rocp, alvlocp, &
                       grav, erad
use mem_grid,    only: mza, lpv, lpw, vnx, vny, vnz, xev, yev, zev, &
                       dzt_top, dzt_bot
use olam_mpi_atm,only: mpi_send_w, mpi_recv_w, &
                       mpi_send_v, mpi_recv_v
use obnd,        only: lbcopy_v, lbcopy_w
use therm_lib,   only: rhovsl

implicit none

integer :: j,iw,k,ka,iv,iter,iw1,iw2,kbc,mrl

real :: temp, exner
real :: uv01dx, uv01dy, uv01dz, uv01dr, raxis, rhovs
real :: pkhyd

! Choose as an internal pressure boundary condition the pressure level at or
! below (in elevation) the 49900 Pa surface.  Find the k index of this level.

kbc = mza
do while (pr01d(kbc) < 49900.)
   kbc = kbc - 1
enddo

!----------------------------------------------------------------------
do j = 1,jtab_w(jtw_init)%jend(1); iw = jtab_w(jtw_init)%iw(j)
!---------------------------------------------------------------------

   ka = lpw(iw)

! Fill arrays prior to iteration

   wc (1:mza,iw) = 0.
   wmc(1:mza,iw) = 0.

   do k = ka, mza
      theta(k,iw) = th01d(k)
      thil(k,iw)  = theta(k,iw)
      press(k,iw) = pr01d(k)
      rho(k,iw)   = dn01d(k)

      if (miclevel == 0) then
         sh_w(k,iw) = 0.
         sh_v(k,iw) = 0.
      elseif (miclevel == 1) then
         sh_w(k,iw) = rt01d(k)
         sh_v(k,iw) = rt01d(k)
      else
         sh_w(k,iw) = rt01d(k)
         sh_v(k,iw) = rt01d(k)
      endif
   enddo

   do iter = 1,100

!  Compute density for all grid levels

      do k = ka, mza

         if (miclevel == 0) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k / (rdry * theta(k,iw))
         elseif (miclevel == 1) then
            rho(k,iw) = press(k,iw) ** cvocp * p00k &
               / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))
         else
            exner = (press(k,iw) * p00i) ** rocp  ! Defined WITHOUT CP factor
            temp = exner * theta(k,iw)
            rhovs = rhovsl(temp-273.15)

            sh_c(k,iw) = max(0., sh_w(k,iw) - rhovs/real(rho(k,iw)))
            sh_v(k,iw) = sh_w(k,iw) - sh_c(k,iw)

            rho(k,iw) = press(k,iw) ** cvocp * p00k &
               / (theta(k,iw) * (rdry * (1. - sh_w(k,iw)) + rvap * sh_v(k,iw)))

            thil(k,iw) = theta(k,iw) / (1. + alvlocp * sh_c(k,iw) / max(temp,253.))
         endif

      enddo

! Integrate hydrostatic equation upward and downward from kbc level
! Impose minimum value of 0.1 Pa to avoid overshoot to negative values
! during iteration.  Use weighting to damp oscillations

      do k = kbc+1,mza
         pkhyd = press(k-1,iw) &
               - grav * (rho(k-1,iw) * dzt_top(k-1) + rho(k,iw) * dzt_bot(k))
         press(k,iw) = .05 * press(k,iw) + .95 * max(.1,pkhyd)
      enddo

      do k = kbc-1,ka,-1
         pkhyd = press(k+1,iw) &
               + grav * (rho(k+1,iw) * dzt_bot(k+1) + rho(k,iw) * dzt_top(k))
         press(k,iw) = .05 * press(k,iw) + .95 * max(.1,pkhyd)
      enddo

   enddo

   do k = ka,mza
      tair(k,iw) = theta(k,iw) * (press(k,iw) * p00i) ** rocp
   enddo

   do k = 1, ka-1
      thil(k,iw) = thil(ka,iw)
   enddo

enddo

! LBC copy (THETA and TAIR will be copied later with the scalars)

mrl = 1

if (iparallel == 1) then
   call mpi_send_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc,rvara2=wmc,rvara3=thil)
   call mpi_recv_w(mrl, dvara1=press, dvara2=rho, &
                   rvara1=wc,rvara2=wmc,rvara3=thil)
endif

call lbcopy_w(1, a1=wc, a2=wmc, a3=thil, d1=press, d2=rho)

! Initialize VMC, VC

!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_init)%jend(1); iv = jtab_v(jtv_init)%iv(j)
   iw1 = itab_v(iv)%iw(1); iw2 = itab_v(iv)%iw(2)
!----------------------------------------------------------------------

   if (iw1 == 1) iw1 = iw2
   if (iw2 == 1) iw2 = iw1

   ka = lpv(iv)

! If sounding winds are to be interpreted as eastward (U) and
! northward (V) components, rotate winds from geographic to
! polar stereographic orientation

! V point coordinates and normal vector components

   do k = ka,mza

      if (mdomain <= 1) then  ! Model uses "earth" coordinates
         raxis = sqrt(xev(iv) ** 2 + yev(iv) ** 2)  ! dist from earth axis

         if (raxis > 1.e3) then
            uv01dr = -v01d(k) * zev(iv) / erad  ! radially outward from axis

            uv01dx = (-u01d(k) * yev(iv) + uv01dr * xev(iv)) / raxis
            uv01dy = ( u01d(k) * xev(iv) + uv01dr * yev(iv)) / raxis
            uv01dz =   v01d(k) * raxis / erad

            vc(k,iv) = uv01dx * vnx(iv) + uv01dy * vny(iv) + uv01dz * vnz(iv)
         else
            vc(k,iv) = 0.
         endif

      else
         vc(k,iv) = u01d(k) * vnx(iv) + v01d(k) * vny(iv)
      endif

      vmc(k,iv) = vc(k,iv) * .5 * (rho(k,iw1) + rho(k,iw2))

   enddo

! For below-ground points, set VC to 0

   vc (1:ka-1,iv) = 0.0
   vmc(1:ka-1,iv) = 0.0

enddo

! Set VMC, VC = 0 at channel (non-topo) walls

!----------------------------------------------------------------------
do j = 1,jtab_v(jtv_wall)%jend(1); iv = jtab_v(jtv_wall)%iv(j)
!----------------------------------------------------------------------

   vmc(:,iv) = 0.
   vc (:,iv) = 0.

enddo

! MPI parallel send/recv of V group

if (iparallel == 1) then
   call mpi_send_v(mrl, rvara1=vmc, rvara2=vc)
   call mpi_recv_v(mrl, rvara1=vmc, rvara2=vc)
endif

! LBC copy of VMC, VC

call lbcopy_v(1, vmc=vmc, vc=vc)

! Set VMP and VP

if (allocated(vmp)) vmp(:,:) = vmc(:,:)
if (allocated(vp )) vp (:,:) = vc (:,:)

end subroutine fldshhi
