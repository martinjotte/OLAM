MODULE module_cu_g3_tracer

  private
  public :: grell_mix_driver

CONTAINS

  subroutine grell_mix_driver(iw, dtl, ktf)

    use mem_cuparm,  only: kcutop, kcubot, kudbot, cbmf, iactcu, &
                           cddf, kddtop, kddmax, kstabi
    use var_tables,  only: scalar_tab, num_cumix, cumix_map
    use mem_basic,   only: rho
    use mem_grid,    only: lpw, arw, arw0, volti
    use mem_ijtabs,  only: itab_w

    implicit none

    integer, intent(in) :: iw, ktf
    real,    intent(in) :: dtl

    real                :: dens(ktf), area(ktf), voa(ktf), aov(ktf)
    real                ::  tr(ktf,num_cumix)
    real                :: dtr(ktf,num_cumix)

    real                :: a0, xmb, edt
    integer             :: ktop, kbcon, k22, kstab, jmin, kdet
    integer             :: ka, k, kc, n, ns, nt, itype
    logical             :: idd

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

    a0 = arw0(iw)
    nt = num_cumix

    do kc = 1, ktf
       k = kc + ka - 1
       dens(kc) = real(rho(k,iw))
       area(kc) = min(arw(k-1,iw), a0)
       aov (kc) = volti(k,iw) * arw0(iw) ! effective 1/dz
       voa (kc) = 1.0 / aov(kc)          ! effective dz
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

    call cup_tracer_mix(tr,dtr,idd,voa,aov,area,a0,xmb,k22,kbcon, &
                        ktop,kstab,edt,jmin,kdet,ktf,nt,itype,itab_w(iw)%iwglobe)

    ! Compute tendencies due to convective fluxes

    do ns = 1, num_cumix
       n = cumix_map(ns)
       do kc = 1, ktop
          k = kc + ka - 1
          scalar_tab(n)%var_t(k,iw) = scalar_tab(n)%var_t(k,iw) + dtr(kc,ns)
       enddo
    enddo

  end subroutine grell_mix_driver


  subroutine cup_tracer_mix(tr,dtr,idd,voa,aov,arw,arw0,xmb,k22,kbcon, &
             ktop,kstabi,edt,jmin,kdet,ktf,nt,itype,iw)

    use module_cu_g3, only: cup_up_en_cd, cup_up_nms, cup_dd_nms, cup_up_tracer, &
                            cup_dd_tracer, cup_dellas, cup_dellas_nodd, &
                            radius_deep, entr_deep, detr_deep, entr_shal, detr_shal
    implicit none

    logical, intent(in)  :: idd
    real,    intent(in)  :: xmb, edt
    integer, intent(in)  :: k22,kbcon,ktop,kstabi,jmin,kdet,ktf,nt,itype,iw
    real,    intent(in)  :: voa(ktf), aov(ktf), arw(ktf), arw0
    real,    intent(in)  ::  tr(ktf,nt)
    real,    intent(out) :: dtr(ktf,nt)

    integer              :: ierr, k, n, kfull
    real                 :: entr_rate, detr_rate, arw0m, rad
    real, dimension(ktf) :: cdd, cd, etd, en, zu, zd
    real, dimension(ktf) :: facu1, facu2, facd1, facd2
    real                 :: tc(ktf), tcd(ktf)

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

    endif

    !--- loop over tracers

    do n = 1, nt

       !--- Tracer concentration in updraft

       call cup_up_tracer(k22,tc,tr(:,n),facu1,facu2,ierr,ktop,ktf)

       if (idd) then

          !--- Tracer concentration in downdraft

          call cup_dd_tracer(tcd,tr(:,n),facd1,facd2,ierr,jmin,ktf)

          !--- Compute tendencies with downdrafts

          call cup_dellas(ierr,aov,dtr(:,n),tr(:,n),tc,zu,tcd,zd,ktop,ktf)

       else

          !--- Compute tendencies without downdrafts

          call cup_dellas_nodd(ierr,aov,dtr(:,n),tr(:,n),tc,zu,k22,ktop,ktf)

       endif

    enddo

  end subroutine cup_tracer_mix


end MODULE module_cu_g3_tracer
