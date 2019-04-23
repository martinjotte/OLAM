module buoyancy

contains

! This routine determines the coefficients A and B that are used to compute
! the buoyancy flux or the buoyancy gradient as a function of the ice-liquid
! potential temperature and total water fluxes:
!
! < w'theta_v' > = A * < w'theta_il' > + B * < w'q_t' >
!
! dB / dz = A d(theta_il) / dz + B d(Q_t) / dz
!
! The formulation of A and B is taken from Cuijpers and Duynkerke (1993, JAS)
! and is based on earlier work by Deardorff (1976, 1980). For unsaturated
! air the coefficients are given by ADRY and BDRY and is based on the standard
! definition of virtual potential temperature; for saturated air the
! coefficients are given as AWET and BWET and take into account latent heating.
!
! The coefficients A and B are computed at flux levels by weighting the wet and
! dry coefficients using the cloud fraction. When the cloud fraction between two
! layers is different, the mixture of air between the two adjacent layers is
! examined to determine if the mixture is saturated. If the mixture is saturated
! then the maximum of the two adjacent cloud fractions is used, else if dry
! the minimum of the two adjacent cloud fractions is used in the weighting.
!
! Subgrid saturation (currently only by the convective schemes) is included
! in the cloud fraction and the calculation of the wet coefficients.
!
! The scheme returns the wet, dry, and/or total A and B coeficients and
! optionally the bouyancy gradient.

subroutine comp_buoy(iw, buoy, a, b, adry, bdry, awet, bwet)

  use mem_grid,    only: mza, dzim, lpw
  use mem_basic,   only: thil, theta, tair, sh_v, sh_w, rho
  use consts_coms, only: t00, eps_vapi, eps_virt, alvl, rdry, rvap, alvlocp
  use therm_lib,   only: rhovsl
  use mem_radiate, only: cloud_frac
  use mem_cuparm,  only: iactcu, qwcon

  implicit none

  integer,        intent( in) :: iw
  real, optional, intent(out) :: buoy(mza) ! vertical gradient of buoyancy
                                           ! variable (dB/dz)
  real, optional, intent(out) :: a   (mza), b   (mza)
  real, optional, intent(out) :: awet(mza), bwet(mza)
  real, optional, intent(out) :: adry(mza), bdry(mza)

  real, parameter      :: alvlordry = alvl / rdry
  real, parameter      :: alvlorvap = alvl / rvap
  real, parameter      :: c_aw      = alvlocp * alvlorvap

  integer              :: k
  real                 :: qsat(mza), qliq(mza), frac(mza)
  real                 :: qt, qs, qc, aa, bb, ad, bd, aw, bw
  real                 :: th, ta, dthl, dqw, thl, xsat, fracm

  ! loop over T layers
  do k = lpw(iw), mza
     qsat(k) = rhovsl(tair(k,iw)-t00) / real(rho(k,iw))

     qliq(k) = sh_w(k,iw) - sh_v(k,iw)
     if (iactcu(iw) > 0) qliq(k) = qliq(k) + qwcon(k,iw)

     if (qliq(k) < 1.e-8) then
        frac(k) = 0.0
     else
        frac(k) = max(0.1, cloud_frac(k,iw))
     endif
  enddo

  ! loop over W levels
  do k = lpw(iw), mza-1

     ! Values at flux levels

     qt  = 0.5 * (sh_w (k,iw) + sh_w (k+1,iw))
     th  = 0.5 * (theta(k,iw) + theta(k+1,iw))
     thl = 0.5 * (thil (k,iw) + thil (k+1,iw))
     qc  = 0.5 * (qliq (k)    + qliq (k+1)   )
     qs  = 0.5 * (qsat (k)    + qsat (k+1)   )

     dthl = thil(k+1,iw) - thil(k,iw)
     dqw  = sh_w(k+1,iw) - sh_w(k,iw)

     ta = 0.5 * ( tair(k  ,iw) / theta(k  ,iw) &
                + tair(k+1,iw) / theta(k+1,iw) ) * th

     ! Dry coefficients

     ad = eps_virt * qt + 1.0
     bd = eps_virt * th

     ! Wet coefficients

     aw = (1.0 - qt + qs * ( eps_vapi + alvlordry / ta )) &
        / (1.0 + c_aw  * qs / ta**2)
     bw = th * (aw * alvlocp / ta - 1.0)

     ! Is the mixture of air between the adjacent layers saturated?

     xsat = qc * (ad * alvlocp - eps_vapi * thl) &
          / max(abs( (ad - aw) * dthl + (bd - bw) * dqw ), 1.e-7)

     if (xsat > 0.5) then
        fracm = max( frac(k), frac(k+1) )
     else
        fracm = min( frac(k), frac(k+1) )
     endif

     ! Weight the wet and dry coefficients to get the total buoyancy

     aa = fracm * aw + (1. - fracm) * ad
     bb = fracm * bw + (1. - fracm) * bd

     if (present(a)) a(k) = aa
     if (present(b)) b(k) = bb

     if (present(awet)) awet(k) = aw
     if (present(bwet)) bwet(k) = bw

     if (present(adry)) adry(k) = ad
     if (present(bdry)) bdry(k) = bd

     if (present(buoy)) then
        buoy(k) = ( aa * dthl + bb * dqw ) * dzim(k)
     endif

  enddo

  if (present(a)) a(mza) = a(mza-1)
  if (present(b)) b(mza) = b(mza-1)

  if (present(adry)) adry(mza) = adry(mza-1)
  if (present(bdry)) bdry(mza) = bdry(mza-1)

  if (present(awet)) awet(mza) = awet(mza-1)
  if (present(bwet)) bwet(mza) = bwet(mza-1)

  if (present(buoy)) buoy(mza) = buoy(mza-1)

end subroutine comp_buoy

end module buoyancy
