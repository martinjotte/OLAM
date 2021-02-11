MODULE module_cu_g3_tracer

  private
  public :: grell_mix_driver

CONTAINS

  subroutine grell_mix_driver(iw, dtl)

    use mem_cuparm,  only: kcutop, kcubot, kudbot, cbmf, iactcu, &
                           cddf, kddtop, kddmax, kstabi
    use var_tables,  only: scalar_tab, num_cumix, cumix_map
    use mem_basic,   only: rho
    use mem_grid,    only: lpw, arw, arw0, arw0i, volt

    implicit none

    integer, intent(in) :: iw
    real,    intent(in) :: dtl

    real, allocatable :: dens(:), area(:), voa(:)
    real, allocatable :: tr(:,:), dtr(:,:)

    real              :: a0, xmb, edt
    integer           :: ktf, ktop, kbcon, k22, kstab, jmin, kdet
    integer           :: ka, k, kc, n, ns, nt, itype
    logical           :: idd

    ! exit if no convection
    if (iactcu(iw) < 1 .or. cbmf(iw) < 1.e-7) return

    ! exit if no scalars to mix
    if (num_cumix < 1) return

    ka = lpw(iw)

    ktop  = kcutop(iw) - ka + 1
    kbcon = kcubot(iw) - ka + 1
    k22   = kudbot(iw) - ka + 1
    kstab = kstabi(iw) - ka + 1
    xmb   = cbmf(iw)

    if (iactcu(iw) == 1 .and. cddf(iw) >= 1.e-7) then
       idd  = .true.
       edt  = cddf(iw) / xmb
       jmin = kddtop(iw) - ka + 1
       kdet = kddmax(iw) - ka + 1
    else
       idd = .false.
       edt   = 0.0
       jmin  = 1
       kdet  = 1
    endif

    ktf = ktop + 2

    allocate(dens(ktf))
    allocate(area(ktf))
    allocate(voa (ktf))
    allocate( tr (ktf,num_cumix))
    allocate(dtr (ktf,num_cumix))

    a0 = arw0(iw)
    nt = num_cumix

    do kc = 1, ktf
       k = kc + ka - 1
       dens(kc) = real(rho(k,iw))
       area(kc) = min(arw(k-1,iw), a0)
       voa (kc) = volt(k,iw) * arw0i(iw)
    enddo
    area(1) = area(2)

    if     (iactcu(iw) == 1 .or. iactcu(iw) == 2) then
       itype = 1
    elseif (iactcu(iw) == 3) then
       itype = 2
    else
       return
    endif

    do ns = 1, num_cumix
       n = cumix_map(ns)
       do kc = 1, ktf
          k = kc + ka - 1
          tr(kc,ns) = scalar_tab(n)%var_p(k,iw) &
                    + scalar_tab(n)%var_t(k,iw) * dtl / dens(kc)
       enddo
    enddo

    call cup_tracer_mix(tr,dtr,idd,dens,voa,area,a0,xmb,k22,kbcon,ktop, &
                        kstab,edt,jmin,kdet,ktf,nt,itype)

    ! Compute tendencies due to convective fluxes

    do ns = 1, num_cumix
       n = cumix_map(ns)
       do kc = 1, ktop
          k = kc + ka - 1
          scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) &
                                    + dtr(k,ns) * dens(kc)
       enddo
    enddo

  end subroutine grell_mix_driver


  subroutine cup_tracer_mix(tr,dtr,idd,dens,voa,arw,arw0,xmb,k22,kbcon, &
             ktop,kstabi,edt,jmin,kdet,ktf,nt,itype)

    use module_cu_g3

    implicit none

    logical, intent(in)  :: idd
    real,    intent(in)  :: xmb, edt
    integer, intent(in)  :: k22, kbcon,ktop,kstabi,jmin,kdet,ktf,nt,itype
    real,    intent(in)  :: voa(ktf), arw(ktf), dens(ktf), arw0
    real,    intent(in)  ::  tr(ktf,nt)
    real,    intent(out) :: dtr(ktf,nt)

    integer              :: ierr, k, kfull
    real                 :: entr_rate, detr_rate, arw0m, rad
    real, dimension(ktf) :: cdd, cd, etd, en, zu, zd
    real, dimension(ktf) :: facu1, facu2, facd1, facd2
    real                 :: tr_cup(ktf,nt), tc(ktf,nt), tcd(ktf,nt)

    ierr = 0
!
!--- gross entrainment rate (these may be changed later on in the
!--- program, depending what your detrainment is!!)
!
    if (itype == 1) then

       rad = sqrt(arw0)
       if (rad < radius_deep) then
          entr_rate = 0.24 / rad
          detr_rate = detr_deep + (entr_rate - entr_deep)
       else
          entr_rate  = entr_deep
          detr_rate  = detr_deep
       endif

    elseif (itype == 2) then

       entr_rate = entr_shal
       detr_rate = detr_shal

    else

       return

    endif
!
!--- Special for shaved cells
!
    arw0m = nearest(arw0,-1.)
    do k = 1, ktf-1
       if (arw(k) >= arw0m) exit
    enddo
    kfull = k
!
!--- Compute entrainment and detrainment
!
    call cup_up_en_cd(en,cd,facu1,facu2,voa,arw,entr_rate,detr_rate,kfull,ierr,ktf)
!
!--- Increase detrainment in stable layers for deep/mid convection
!
    if (itype == 1 .and. kstabi <= ktop) then
       do k = kstabi, ktop
          cd   (k) = min(cd(k-1)/voa(k-2)+.15*entr_rate,1.5*entr_rate)*voa(k-1)
          facu2(k) = en(k)  / (1.0 + en(k) - 0.5*cd(k))
          facu1(k) = 1.0 - facu2(k)
       enddo
    endif
!
!--- Normalized updraft mass flux profile
!
    call cup_up_nms(zu,cd,en,ktop,ierr,k22,facu1,facu2,arw,arw0,ktf)
!
!--- Scale updrafts
!
    do k = k22, ktop
       zu(k) = zu(k) * xmb
    enddo
!
!--- Tracer concentration at cloud (mid-point) levels
!
    call cup_env_clev_tr(tr,tr_cup,ktf,nt)
!
!--- Tracer concentration in updraft
!
    call cup_up_tracer(k22,facu1,facu2,tr_cup,tr,tc,ktop,ierr,ktf,nt)

    if (idd) then
!
!--- Downdraft mass flux and entrainment/detrainment rates
!
       call cup_dd_nms(zd,voa,cdd,etd,facd1,facd2,entr_rate, &
                       jmin,ierr,arw,arw0,kdet,kfull,ktf)
!
!--- Scale downdrafts
!
       do k = 2, jmin
          zd(k) = zd(k) * xmb * edt
       enddo
!
!--- Tracer concentration in downdraft
!
       call cup_dd_tracer(tcd,facd1,facd2,jmin,ierr,entr_rate,tr,tr_cup,nt,ktf)
!
!--- Compute tendencies with downdrafts
!
      call cup_dellas_tr(ierr,voa,dtr,zu,tr,tr_cup,tc,dens, &
           ktop,k22,cd,en,tcd,zd,cdd,etd,jmin,nt,ktf)

   else
!
!--- Compute tendencies without downdrafts
!
      call cup_dellas_tr_nodd(ierr,voa,dtr,zu,tr,tr_cup,tc,dens, &
           ktop,k22,cd,en,nt,ktf)

    endif

  end subroutine cup_tracer_mix


  SUBROUTINE cup_up_tracer(k22,fact1,fact2,tr_cup,tr,tc,ktop,ierr,ktf,nt)

    implicit none
!
!  on input
!
    integer                                                           &
       ,intent (in   )                   ::                           &
       ktf,nt,k22,ktop,ierr
    real,    dimension (ktf,nt)                                       &
       ,intent (in   )                   ::                           &
      tr,tr_cup
    real,    dimension (ktf)                                          &
      ,intent (in   )                    ::                           &
      fact1,fact2
!
! input and output
!
    real,    dimension (ktf,nt)                                       &
       ,intent (out  )                   ::                           &
       tc
!
!  local variables in this routine
!
    integer                              ::                           &
       k,n

    tc = 0.
    if (ierr.ne.0) return

    do n = 1, nt
       tc(k22,n) = tr_cup(k22,n)
    enddo

    do n = 1, nt
       do k = k22+1, ktop+1
          tc(k,n) = tc(k-1,n) * fact1(k) + tr(k-1,n) * fact2(k)
       enddo
    enddo

  END SUBROUTINE cup_up_tracer


  SUBROUTINE cup_env_clev_tr(tr,tr_cup,ktf,nt)
    implicit none

    integer, intent(in)  :: ktf, nt
    real,    intent(in)  :: tr    (ktf,nt)
    real,    intent(out) :: tr_cup(ktf,nt)
    integer              :: k, n

    do n = 1, nt
       tr_cup(1,n) = tr(1,n)
       do k = 2, ktf
          tr_cup(k,n) = 0.5 * (tr(k-1,n) + tr(k,n))
       enddo
    enddo

  END SUBROUTINE cup_env_clev_tr



  SUBROUTINE cup_dellas_tr(ierr,dz_cup,della,zu,tr,tr_cup,tc,dens,  &
             ktop,k22,cd,en,tcd,zd,cdd,etd,jmin,ntr,ktf)

    implicit none

    integer                                                           &
       ,intent (in   )                   ::                           &
       ktf,ntr,ierr,ktop,k22,jmin
    real,    dimension (ktf)                                          &
       ,intent (in  )                    ::                           &
       dz_cup,cd,en,zu,dens,cdd,etd,zd
    real,    dimension (ktf,ntr)                                      &
       ,intent (in   )                   ::                           &
       tr,tr_cup,tc,tcd
    real,    dimension (ktf,ntr)                                      &
       ,intent (out  )                   ::                           &
       della
!
!  local variables in this routine
!
    integer :: k, n
    real, dimension(ktf) :: subin,subdn,entra,detup,detdo,rhodzm

    della = 0.

    IF (ierr.ne.0) return

    do k = 1, ktop
       rhodzm(k) = 1.0 / (dens(k) * dz_cup(k))
    enddo

    do n = 1, ntr
       della(1,n) = zd(2) * (.5*(tcd(1,n)+tcd(2,n)) - tr_cup(2,n)) * rhodzm(1)
    enddo
!
!--- SPECIFY DETRAINMENT OF DOWNDRAFT, HAS TO BE CONSISTENT
!--- WITH ZD CALCULATIONS IN SOUNDD.
!
    do k = 2, ktop
       subin(k) = zu(k+1) - zd(k+1)

       subdn(k) = zu(k)   - zd(k)

       entra(k) = en(k+1) * zu(k) + etd(k) * zd(k+1)

       detup(k) = cd(k+1) * zu(k)

       detdo(k) = cdd(k)  * zd(k+1)
    enddo

    do n = 1, ntr

       do k = 2, ktop
          della(k,n) = ( detup(k)*.5*(tc(k+1,n)+tc(k,n)) &
                       + detdo(k)*.5*(tcd(k+1,n)+tcd(k,n)) &
                       - entra(k)*tr(k,n) &
                       + subin(k)*tr_cup(k+1,n) &
                       - subdn(k)*tr_cup(k,n) &
                       ) * rhodzm(k)
       enddo

       ! Special at jmin
       della(jmin,n) = della(jmin,n) - zd(jmin) * tcd(jmin+1,n) * rhodzm(jmin)

       ! Special at k22
       della(k22-1,n) = della(k22-1,n) - zu(k22) * tc(k22,n) * rhodzm(k22-1)

    enddo

  END SUBROUTINE cup_dellas_tr



  SUBROUTINE cup_dellas_tr_nodd(ierr,dz_cup,della,zu,tr,tr_cup,tc,dens,  &
                                ktop,k22,cd,en,ntr,ktf)

    implicit none

    integer                                                           &
       ,intent (in   )                   ::                           &
       ktf,ntr,ierr,ktop,k22
    real,    dimension (ktf)                                          &
       ,intent (in  )                    ::                           &
       dz_cup,cd,en,zu,dens
    real,    dimension (ktf,ntr)                                      &
       ,intent (in   )                   ::                           &
       tr,tr_cup,tc
    real,    dimension (ktf,ntr)                                      &
       ,intent (out  )                   ::                           &
       della
!
!  local variables in this routine
!
    integer :: k, n
    real, dimension(ktf) :: entra, detup, rhodzm

    della = 0.

    IF (ierr.ne.0) return
!
!--- SPECIFY DETRAINMENT OF DOWNDRAFT, HAS TO BE CONSISTENT
!--- WITH ZD CALCULATIONS IN SOUNDD.
!
    do k = k22, ktop
       entra (k) = en(k+1) * zu(k)
       detup (k) = cd(k+1) * zu(k)
       rhodzm(k) = 1.0 / (dens(k) * dz_cup(k))
    enddo

    do n = 1, ntr
       do k = k22, ktop
          della(k,n) = ( detup(k)*.5*(tc(k+1,n)+tc(k,n)) &
                       - entra(k)*tr(k,n) &
                       + zu(k+1)*tr_cup(k+1,n) &
                       - zu(k)*tr_cup(k,n) &
                       ) * rhodzm(k)
       enddo
    enddo

  END SUBROUTINE cup_dellas_tr_nodd


  SUBROUTINE cup_dd_tracer(tcd,fact1,fact2,jmin,ierr,entr,tr,tr_cup,ntr,ktf)

    implicit none
!
!  on input
!
    integer                                                           &
       ,intent (in   )                   ::                           &
       ktf,ntr,jmin,ierr
    real,    dimension (ktf,ntr)                                      &
       ,intent (in   )                   ::                           &
       tr,tr_cup
    real,    dimension (ktf)                                          &
       ,intent (in   )                   ::                           &
       fact1,fact2
    real                                                              &
       ,intent (in   )                   ::                           &
       entr
!
! input and output
!
    real,    dimension (ktf,ntr)                                      &
       ,intent (out  )                   ::                           &
       tcd
!
!  local variables in this routine
!
    integer                              ::                           &
       k,n

    tcd = 0.0
    IF (ierr.ne.0) return

    if (entr <= 1.e-20) then

       do n = 1, ntr
          do k = 1, jmin+1
             tcd(k,n) = tr_cup(jmin,n)
          enddo
       enddo

    else

       do n = 1, ntr
          tcd(jmin+1,n) = tr_cup(jmin,n)
          tcd(jmin  ,n) = tr_cup(jmin,n)

          do k = jmin-1, 2, -1
             tcd(k,n) = tcd(k+1,n) * fact1(k) + tr(k,n) * fact2(k)
          enddo

          tcd(1,n) = tcd(2,n)
       enddo

    endif

  END SUBROUTINE cup_dd_tracer


end MODULE module_cu_g3_tracer
