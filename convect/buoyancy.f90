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
  use mem_basic,   only: theta, tair, rr_v, rr_w, rho
  use consts_coms, only: t00, eps_vapi, eps_virt, alvl, rdry, rvap, &
                         alvlocp, alviocp
  use therm_lib,   only: rhovsl
  use mem_radiate, only: cloud_frac
  use mem_cuparm,  only: iactcu, qwcon, kcutop, kcubot
  use mem_micro,   only: rr_c, rr_p
  use micro_coms,  only: miclevel
  use oname_coms,  only: nl

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

  real                 :: qsat(mza), qliq(mza), qice(mza), frac(mza)
  real                 :: thlo(mza), qto(mza)
  real                 :: qt, qs, qc, aa, bb, ad, bd, aw, bw
  real                 :: th, ta, dthl, dqw, thl, xsat, fracm, fac
  integer              :: k

!  moist_buoy = 0  buoyancy does not take into account condensate
!               1  resolved cloud droplets only
!               2  resolved could droplets and pristine ice
!               3  resolved and convective could droplets and pristine ice

  qliq(:) = 0.
  qice(:) = 0.

  if (nl%moist_buoy >= 1) then
     if (miclevel < 3) then
        do k = lpw(iw), mza
           qliq(k) = rr_w(k,iw) - rr_v(k,iw)
        enddo
     else
        if (allocated(rr_c)) then
           do k = lpw(iw), mza
              qliq(k) = rr_c(k,iw)
           enddo
        endif
        if (nl%moist_buoy >= 2 .and. allocated(rr_p)) then
           do k = lpw(iw), mza
              qice(k) = rr_p(k,iw)
           enddo
        endif
     endif

     if (nl%moist_buoy >= 3 .and. iactcu(iw) > 0) then
        do k = kcubot(iw), kcutop(iw)
           qliq(k) = qliq(k) + qwcon(k,iw)
        enddo
     endif
  endif

  ! loop over T levels
  do k = lpw(iw), mza

     fac = 1.0 / (1.0 + rr_v(k,iw))

     qliq(k) = qliq(k) * fac
     qice(k) = qice(k) * fac
     qsat(k) = rhovsl(tair(k,iw)-t00) / real(rho(k,iw)) * fac

     if (qliq(k) + qice(k) < 1.e-8) then
        frac(k) = 0.0
     else
        frac(k) = max(0.1, cloud_frac(k,iw))
     endif

     thlo(k) = theta(k,iw) / ( 1.0 + (alvlocp * qliq(k) + alviocp * qice(k)) / max(tair(k,iw),253.0))
     qto (k) = rr_v(k,iw) * fac + qliq(k) + qice(k)
  enddo

  ! loop over W levels
  do k = lpw(iw), mza-1

     thl = 0.5 * (thlo(k) + thlo(k+1))
     qs  = 0.5 * (qsat(k) + qsat(k+1))
     qt  = 0.5 * (qto (k) + qto (k+1))
     th  = 0.5 * (theta(k,iw) + theta(k+1,iw))
     qc  = 0.5 * (qliq(k) + qice(k) + qliq(k+1) + qice(k+1))

     dthl = thlo(k+1) - thlo(k)
     dqw  = qto (k+1) - qto (k)

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

     aa = ad
     bb = bd

     if (qc > 1.e-7) then
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
     endif

     ! Output variables

     if (present(a)) a(k) = aa
     if (present(b)) b(k) = bb

     if (present(awet)) awet(k) = aw
     if (present(bwet)) bwet(k) = bw

     if (present(adry)) adry(k) = ad
     if (present(bdry)) bdry(k) = bd

     if (present(buoy)) buoy(k) = (aa * dthl + bb * dqw) * dzim(k)

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
