subroutine cuparm_emanuel(iw, dtlong)

  use mem_grid,    only: mwa, mza, lpu, lpw, volt, mua, zm, zt
  use mem_tend,    only: thilt, sh_wt
  use mem_basic,   only: theta, tair, press, rho, sh_v, vxe, vye, vze
  use consts_coms, only: t00, grav
  use misc_coms,   only: confrq
  use mem_cuparm , only: thsrc, rtsrc, conprr, cbmf
  implicit none

  integer, intent(in)  :: iw
  real,    intent(in)  :: dtlong

  real, dimension(mza) :: tc, qc, qsc, vx, vy, vz, pc, pfc, qcldc
  real, dimension(mza) :: tt, qt, vxt, vyt, vzt, gz

  integer :: k, ka, kc, nd, na, nl
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

  do k = ka, mza-1
     kc = k - ka + 1 ! relative to surface

     gz(kc) = grav * (zt(k) - zt(ka)) ! geopotential height relative to 1st level
     
     tc(kc)  = tair(k,iw)
     qc(kc)  = sh_v(k,iw)
     pc(kc)  = press(k,iw)*0.01

     vx(kc) = vxe(k,iw)
     vy(kc) = vye(k,iw)
     vz(kc) = vze(k,iw)
  enddo

  do k = ka, mza-1
     kc = k - ka + 1
     qsc(kc) = rhovsil(tc(kc)-t00) / rho(k,iw)
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
       tc,  qc,   qsc,    vx,  vy, vz,       scalc,  pc,  pfc, gz, &
       nd,  nl,     ntra, dtlong,   iflag,  tt, qt, vxt, &
       vyt, vzt, scalt,   pcprate, wprime, tprime, qprime, cbmf(iw), qcldc, &
       kcbase, kctop)

!!  if (iflag == 4) then
!!     write(*,'(A)')    " CFL Error in the Emanuel cumulus scheme"
!!     write(*,'(A,I0)') " at iw = ", iw
!!  endif

  if (iflag == 1 .or. iflag == 4) then

     conprr(iw) = pcprate

     do k = ka, mza-1
        kc = k - ka + 1
        thsrc(k,iw) = tt(kc) * theta(k,iw) / tair(k,iw)
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
