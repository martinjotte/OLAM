subroutine cuparm_emanuel(iw, dtlong)

!iu1, iu2, iu3, thtend, qtend, pcprate, &
!                          u1tend, u2tend, u3tend)

  use mem_grid,    only: mwa, mza, lpu, lpw, volt, mua, zm, zt
  use mem_tend,    only: thilt, sh_wt
  use mem_basic,   only: theta, press, rho, sh_v, vxe, vye, vze
  use consts_coms, only: p00i, t00, rocp, grav
  use misc_coms,   only: confrq
  use mem_cuparm , only: thsrc, rtsrc, conprr, cbmf

!!, mza_ut, &
!!                         mza_st, cbmf, qcsubg, nqparm, confrq, condelay, &
!!                         eman_mixu, eman_mixscal
  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: dtlong
  
!!
!!  real,    intent(out) :: thtend(mza), qtend(mza), pcprate
!!  real,    intent(out) :: u1tend(mza_ut), u2tend(mza_ut), u3tend(mza_ut)

  real, dimension(mza) :: tc, qc, qsc, vx, vy, vz, pc, pfc, qcldc, exner
  real, dimension(mza) :: tt, qt, vxt, vyt, vzt, qcldw

  integer :: k, n, ka, k2, kc, ka1, ka2, ka3, koff, na, nd, nl
  real, external :: rhovsl, rhovsil

  ! temporary for scalars
  integer, parameter :: ntra = 0
  real :: scalc(mza,ntra), scalt(mza,ntra)

  real    :: pcprate, wprime, tprime, qprime, pb, pt
  integer :: iflag, kcbase, kctop

  scalt = 0.
  ka = lpw(iw)

  nd = mza - ka
  na = nd + 1
  nl = nd - 1

  do k=ka,mza-1

     ! Height on convective grid (relative to surface)
     kc = k - ka + 1
     
     ! zc(kc)  = zt(k) - zm(ka-1)
     ! ztc(kc) = zt(k) - zt(ka)
     
     exner(k) = (p00i * press(k,iw)) ** rocp

     tc(kc)  = theta(k,iw) * exner(k)
     qc(kc)  = sh_v(k,iw)
!    qsc(kc) = rhovsl(tc(kc)-t00) / rho(k,iw)
     qsc(kc) = rhovsil(tc(kc)-t00) / rho(k,iw)
     pc(kc)  = press(k,iw)*0.01

     vx(kc) = vxe(k,iw)
     vy(kc) = vye(k,iw)
     vz(kc) = vze(k,iw)
  enddo

  ! pressure at full (W) layers - compute hydrostatically

  pfc(1) =  0.01 * (press(ka,iw) + (zt(ka)-zm(ka-1))*rho(ka,iw)*grav)
  do k = ka, mza-2
     kc = k - ka + 1
     pb = press(k,iw)   + (zt(k)  -zm(k))*rho(k,iw)  *grav
     pt = press(k+1,iw) + (zt(k+1)-zm(k))*rho(k+1,iw)*grav
     pfc(kc+1) = 0.005 * (pb + pt)
  enddo
  pfc(na) = 0.01 * (press(mza-1,iw) + (zt(mza-1)-zm(mza-1))*rho(mza-1,iw)*grav)

  !
  ! Future - add conserved tracer mixing
  ! 
!!  ntra = naddsc
!!  do n = 1, ntra
!!     do kc = ka, mza-1
!!        scalc(kc,n) = 
!!     enddo
!!  enddo
  
  call convect(                                                  &
       tc,  qc,   qsc,    vx,  vy, vz,       scalc,  pc,     pfc, &
       nd,  na,   nl,     ntra, dtlong,   iflag,  tt, qt, vxt, &
       vyt, vzt, scalt,   pcprate, wprime, tprime, qprime, cbmf(iw), qcldc, &
       kcbase, kctop)

!!  if (iflag == 4) then
!!     write(*,'(A)')    " CFL Error in the Emanuel cumulus scheme"
!!     write(*,'(A,I0)') " at iw = ", iw
!!  endif

  if (iflag == 1 .or. iflag == 4) then

     conprr(iw) = pcprate

     do k=ka,mza-1
        kc = k - ka + 1
        thsrc(k,iw) = tt(kc) / exner(k)
        rtsrc(k,iw) = qt(kc)
     enddo

  else

     thsrc(:,iw) = 0.0
     rtsrc(:,iw) = 0.0
     conprr (iw) = 0.0

  endif


     ! convective mixing???

!!     if (eman_mixu == 1) then
!!        do k=ka,mza-1
!!           kc = k - ka + 1
!!           u1tend(k) = vxt(kc)
!!           u2tend(k) = vyt(kc)
!!           u3tend(k) = vzt(kc)
!!        enddo
!!     endif

!!     if (naddsc > 0 .and. eman_mixscal == 1) then
!!        do n=1,naddsc
!!           do k=ka,mza-1
!!           enddo
!!        enddo
!!     endif

!!     if (frac_cld == 1) then
!!        ! qcldw(k)  = qcldc(kc)
!!     enddo


  return
end subroutine cuparm_emanuel
