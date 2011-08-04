!===============================================================================
! OLAM version 4.0

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

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================

real function eslf(t)

!     This function calculates the saturation vapor pressure over liquid water
!     as a function of Celsius temperature following Flatau et al. (JAM, 1992).

implicit none
real, intent(in) :: t
real             :: x
real, parameter  :: c0 = .6105851e+03, c1 = .4440316e+02, c2 = .1430341e+01
real, parameter  :: c3 = .2641412e-01, c4 = .2995057e-03, c5 = .2031998e-05
real, parameter  :: c6 = .6936113e-08, c7 = .2564861e-11, c8 =-.3704404e-13

x = max(-80.,t)
eslf = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
end function eslf

!===============================================================================

real function esif(t)

!     This function calculates the saturation vapor pressure over ice
!     as a function of Celsius temperature following Flatau et al. (JAM, 1992).

implicit none
real, intent(in) :: t
real             :: x
real, parameter  :: c0 = .6114327e+03, c1 = .5027041e+02, c2 = .1875982e+01
real, parameter  :: c3 = .4158303e-01, c4 = .5992408e-03, c5 = .5743775e-05
real, parameter  :: c6 = .3566847e-07, c7 = .1306802e-09, c8 = .2152144e-12

x = max(-80.,t)
esif = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

return
end function esif

!===============================================================================

real function eslpf(t)

!     This function calculates the partial derivative of liquid saturation vapor
!     pressure with respect to temperature as a function of Celsius temperature
!     following Flatau et al. (JAM, 1992).

implicit none
real, intent(in) :: t
real             :: x
real, parameter  :: d0 = .4443216e+02, d1 = .2861503e+01, d2 = .7943347e-01
real, parameter  :: d3 = .1209650e-02, d4 = .1036937e-04, d5 = .4058663e-07
real, parameter  :: d6 =-.5805342e-10, d7 =-.1159088e-11, d8 =-.3189651e-14

x = max(-80.,t)
eslpf = d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))

return
end function eslpf

!===============================================================================

real function esipf(t)

!     This function calculates the partial derivative of ice saturation vapor
!     pressure with respect to temperature as a function of Celsius temperature
!     following Flatau et al. (JAM, 1992).

implicit none
real, intent(in) :: t
real             :: x
real, parameter  :: d0 = .5036342e+02, d1 = .3775758e+01, d2 = .1269736e+00
real, parameter  :: d3 = .2503052e-02, d4 = .3163761e-04, d5 = .2623881e-06
real, parameter  :: d6 = .1392546e-08, d7 = .4315126e-11, d8 = .5961476e-14

x = max(-80.,t)
esipf = d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))

return
end function esipf

!===============================================================================

real function rslf(p,t)

!     This function calculates the liquid saturated mixing ratio as a
!     function of pressure (Pa) and Kelvin temperature using the
!     saturated vapor pressure curves of Flatau et al. (JAM, 1992).

use consts_coms, only: t00, eps_vap
implicit none

real, intent(in) :: p, t
real             :: x, eslf
real, parameter  :: c0 = .6105851e+03, c1 = .4440316e+02, c2 = .1430341e+01
real, parameter  :: c3 = .2641412e-01, c4 = .2995057e-03, c5 = .2031998e-05
real, parameter  :: c6 = .6936113e-08, c7 = .2564861e-11, c8 =-.3704404e-13

x = max(-80.,t-t00)
eslf = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
eslf = min(eslf,0.9*p)
rslf = eps_vap*eslf/(p-eslf)

return
end function rslf

!===============================================================================

real function rsif(p,t)

!     This function calculates the ice saturated mixing ratio as a
!     function of pressure (Pa) and Kelvin temperature using the
!     saturated vapor pressure curves of Flatau et al. (JAM, 1992).

use consts_coms, only: t00, eps_vap
implicit none

real, intent(in) :: p, t
real             :: x, esif
real, parameter  :: c0 = .6114327e+03, c1 = .5027041e+02, c2 = .1875982e+01
real, parameter  :: c3 = .4158303e-01, c4 = .5992408e-03, c5 = .5743775e-05
real, parameter  :: c6 = .3566847e-07, c7 = .1306802e-09, c8 = .2152144e-12

x = max(-80.,t-t00)
esif = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
esif = min(esif,0.9*p)
rsif = eps_vap*esif/(p-esif)

return
end function rsif

!===============================================================================

real function rslfp(p,t)

!     This function calculates the partial derivative of liquid saturated
!     mixing ratio with respect to temperature as a function of pressure (Pa)
!     and Kelvin temperature

use consts_coms, only: t00, eps_vap
implicit none

real, intent(in) :: p, t
real             :: x, eslpf
real, parameter  :: d0 = .4443216e+02, d1 = .2861503e+01, d2 = .7943347e-01
real, parameter  :: d3 = .1209650e-02, d4 = .1036937e-04, d5 = .4058663e-07
real, parameter  :: d6 =-.5805342e-10, d7 =-.1159088e-11, d8 =-.3189651e-14

x = max(-80.,t-t00)
eslpf = d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))
eslpf = min(eslpf,0.9*p)
rslfp = eps_vap*eslpf/(p-eslpf)

return
end function rslfp

!===============================================================================

real function rsifp(p,t)

!     This function calculates the partial derivative of ice saturated
!     mixing ratio with respect to temperature as a function of pressure (Pa)
!     and Kelvin temperature

use consts_coms, only: t00, eps_vap
implicit none

real, intent(in) :: p, t
real             :: x, esipf
real, parameter  :: d0 = .5036342e+02, d1 = .3775758e+01, d2 = .1269736e+00
real, parameter  :: d3 = .2503052e-02, d4 = .3163761e-04, d5 = .2623881e-06
real, parameter  :: d6 = .1392546e-08, d7 = .4315126e-11, d8 = .5961476e-14

x = max(-80.,t-t00)
esipf = d0+x*(d1+x*(d2+x*(d3+x*(d4+x*(d5+x*(d6+x*(d7+x*d8)))))))
esipf = min(esipf,0.9*p)
rsifp = eps_vap*esipf/(p-esipf)

return
end function rsifp

!===============================================================================

real function rhovsl(tc)

!     This function calculates the density of water vapor at saturation
!     over liquid as a function of Celsius temperature

use consts_coms, only: t00
implicit none
real, intent(in) :: tc

real, parameter :: c0 = .6105851e+03 ,c1 = .4440316e+02 ,c2 =  .1430341e+01
real, parameter :: c3 = .2641412e-01 ,c4 = .2995057e-03 ,c5 =  .2031998e-05
real, parameter :: c6 = .6936113e-08 ,c7 = .2564861e-11 ,c8 = -.3704404e-13
real, parameter :: rvap = 461.
real :: esl,x

x = max(-80.,tc)
esl = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rhovsl = esl / (rvap * (tc + t00))

return
end function rhovsl

!===============================================================================

real function rhovsi(tc)

!     This function calculates the density of water vapor at saturation
!     over ice as a function of Celsius temperature

use consts_coms, only: t00
implicit none
real, intent(in) :: tc

real, parameter :: c0 = .6114327e+03 ,c1 = .5027041e+02 ,c2 = .1875982e+01
real, parameter :: c3 = .4158303e-01 ,c4 = .5992408e-03 ,c5 = .5743775e-05
real, parameter :: c6 = .3566847e-07 ,c7 = .1306802e-09 ,c8 = .2152144e-12
real, parameter :: rvap = 461.
real :: esi,x

x = max(-80.,tc)
esi = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))
rhovsi = esi / (rvap * (tc + t00))

return
end function rhovsi

!===============================================================================

real function rhovsil(tc)

implicit none

real, intent(in) :: tc
real, external :: rhovsl,rhovsi

!     This function calculates the density of water vapor at saturation,
!     over liquid or ice depending on temperature, as a function of
!     Celsius temperature

if (tc >= 0.) then
   rhovsil = rhovsl(tc)
else
   rhovsil = rhovsi(tc)
endif

return
end function rhovsil

!===============================================================================

real function td(p,rs)

!     This function calculates the dew point temperature (K) using Bolton's
!     formula (MWR, 1980) given the atmospheric pressure p (Pa) and water
!     vapor mixing ratio rs (kg/kg)
!
use consts_coms, only: eps_vap, t00
implicit none

real, intent(in) :: p, rs
real             :: rr, es, esln, desdt, dt
integer          :: iter
real, external   :: eslpf, eslf

rr = max(rs, 1.e-8)
es = p*rr / (eps_vap+rr)
esln = log(es)

td = (29.65*esln-5016.7042) / (esln-24.0852)

! We can refine the estimate for td using the eslpf funcion and
! Newton iteration:
!
!if (td-t00 > -70.) then
!   do iter=1,5
!      desdt = eslpf(td-t00)
!      dt = (es - eslf(td-t00))/desdt
!      td = td + dt
!      if (abs(dt).lt. 0.01) exit
!   enddo
!endif

return
end function td

!===============================================================================

subroutine thetae(p,t,rv,the)

use consts_coms, only: cp, gordry, alvlocp, cpog, gocp, p00, eps_virt, rocp
implicit none

real, intent(in)  :: p, t, rv
real, intent(out) :: the

real    :: pit,tupo,ttd,dz,tupn,tmn
integer :: iter
real, external :: td

pit  = p
tupo = t
ttd  = td(p,rv)
dz   = cpog*(t-ttd)
iter = 1

if (dz > 0.0) then

   do iter=1,50
      tupn = t - gocp*dz
      tmn  = 0.5*(tupn+t) * (1.0+eps_virt*rv)
      pit  = p * exp(-gordry*dz/tmn)

      if (abs(tupn-tupo) < 0.001) exit

      ttd  = td(pit,rv)
      tupo = tupn
      dz   = dz + cpog*(tupn-ttd)
   enddo

endif

the = tupo * (p00/pit)**rocp * exp(alvlocp*rv/tupo)

if (iter == 51) then
   print *, "WARNING:"
   print *, "Theta_e computation failed to converge after 50 iterations"
   print *, "in subroutine thetae."
   print *, "P, T, R, theta_e: ", p, t, rv, the
endif

return
end subroutine thetae

!===============================================================================

real function tw(rvp,thet,p)
use consts_coms, only: p00, p00i, rocp, alvlocp
implicit none

real, intent(in) :: rvp, thet, p
real             :: piter, t, x, a, tq, d
integer          :: id
real, external   :: rslf, td

piter = p

do id=1,10
   t = thet * (piter*p00i)**rocp
   x = 0.02 * (td(rvp,piter) - t)
   if ( abs(x) < 0.01 ) exit
   piter = piter * 2.**x
enddo

t  = thet * (piter*p00i)**rocp
a  = thet * exp(alvlocp*rslf(t,piter)/t)
tq = 253.16
d  = 120.0

do id= 1,12
   d = d/2.
   ! if the temperature difference,x, is small,exit this loop
   x = a*exp(-alvlocp*rslf(tq,p)/tq) - tq*(p00/p)**rocp
   if ( abs(x) < 0.01 ) exit
   tq = tq + sign(d,x)
enddo

tw = tq

return
end function tw

!===============================================================================

subroutine the2t(the,p,th,t,r)

use consts_coms, only: p00i, rocp, cp, alvlocp
implicit none

real, intent(in)  :: the, p
real, intent(out) :: th, t, r

real    :: pi, to
integer :: iter
real, external :: rslf

r  = 0.012  ! first guess
to = 295.0  ! first guess

pi = (p*p00i)**rocp
to = the*pi*exp(-alvlocp*r/to)

do iter=1,50
   r  = rslf(p,to)
   th = the*exp(-alvlocp*r/to)
   t  = th*pi
   if (abs(to-t) < 0.005) return
   to = to + (t-to)*.3
enddo

print *,            'WARNING:'
print *,            ' Iteration did not converge in routine the2t.'
print '(a,5e15.6)', ' the,p,to,th,r: ', the, p, to, th, r

return
end subroutine the2t

!===============================================================================

subroutine qtk(q,tempk,fracliq)
use consts_coms, only: t00, alli, allii, cicei, cliqi
implicit none

real, intent(in)  :: q
real, intent(out) :: tempk, fracliq

!     Input:
!        q        internal energy [J/kg]
!     Outputs:
!       tempk    temperature [K]
!       fracliq  liquid fraction [dimensionless]

if (q <= 0.) then
   fracliq = 0.
   tempk = q * cicei + t00
elseif (q >= alli) then
   fracliq = 1.
   tempk = q * cliqi + 193.36
else
   fracliq = q * allii
   tempk = t00
endif

return
end subroutine qtk

!===============================================================================

subroutine qtc(q,tempc,fracliq)

use consts_coms, only: alli, allii, cicei, cliqi
implicit none

real, intent(in)  :: q
real, intent(out) :: tempc, fracliq

!     Input:
!        q        internal energy [J/kg]
!     Outputs:
!        tempc    temperature [C]
!        fracliq  liquid fraction [dimensionless]

if (q <= 0.) then
   fracliq = 0.
   tempc = q * cicei
elseif (q >= alli) then
   fracliq = 1.
   tempc = q * cliqi - 80.
else
   fracliq = q * allii
   tempc = 0.
endif

return
end subroutine qtc

!===============================================================================

subroutine qwtk(qw,w,dryhcap,tempk,fracliq)

use consts_coms, only: t00, alli, allii, cice, cliq
implicit none

real, intent(in)  :: qw, w, dryhcap
real, intent(out) :: tempk, fracliq
real              :: qwliq0

!     Inputs:
!        qw       internal energy [J/m^2] or [J/m^3]
!        w        mass [kg/m^2] or [kg/m^3]
!        dryhcap  heat capacity of nonwater part [J/(m^2 K)] or [J/(m^3 K)]
!     Outputs:
!        tempk    temperature [K]
!        fracliq  liquid fraction [dimensionless]

qwliq0 = w * alli
if (qw <= 0.) then
   fracliq = 0.
   tempk = qw / (cice * w + dryhcap) + t00
elseif (qw >= qwliq0) then
   fracliq = 1.
   tempk = (qw - qwliq0) / (cliq * w + dryhcap) + t00
else
   fracliq = qw / qwliq0
   tempk = t00
endif

return
end subroutine qwtk

