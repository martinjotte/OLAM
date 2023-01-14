subroutine cuparm_driver()

  use mem_grid,         only: vxn_ew, vxn_ns, vxn_ew, vyn_ns, vyn_ew, vzn_ns, &
                              mza, lpw
  use module_cu_g3,     only: grell_driver
  use module_cu_kfeta,  only: cuparm_kfeta, kf_lutab
  use module_cu_tiedtke,only: cuparm_tiedtke
  use module_cu_emanuel,only: cuparm_emanuel
  use misc_coms,        only: io6, time_istp8, time_istp8p, nqparm, confrq, &
                              dtlm, mstp, runtype, iparallel
  use mem_ijtabs,       only: itab_w, jtab_w, mrl_begl, istp, mrls, jtw_prog, jtw_wadj
  use mem_cuparm,       only: thsrc, rtsrc, aconpr, conprr, kstabi, &
                              kcutop, kcubot, qwcon, iactcu, cbmf, cddf, kddtop, &
                              kudbot, kddmax, cu_pwa, cu_pev, umsrc, vmsrc
  use mem_tend,         only: thilt, rr_wt, vmxet, vmyet, vmzet
  use mem_basic,        only: rho, rr_w, theta, tair, thil
  use olam_mpi_atm,     only: mpi_send_w, mpi_recv_w
  use oname_coms,       only: nl
  use var_tables,       only: num_cumix
  use consts_coms,      only: alvlocp, r8
  use olam_mpi_atm,     only: mpi_recv_w, mpi_send_w

  use module_cu_g3_tracer, only: grell_mix_driver

  implicit none

  integer, save :: init_kf = 0
  integer       :: j, iw, k, mrl, mrlw
  real          :: dthmax, dtlong4, dtli, confrq4
  integer       :: iwqmax, kqmax, km, km1
  real          :: rt, ravail

! For KF_eta parameterization, initialize scheme if needed

  if ( any(nqparm(1:mrls) == 3) ) then
     if (init_kf == 0) then
        init_kf = 1
        call kf_lutab()
     endif
  endif

  dtlong4 = real(dtlm)
  dtli    = 1.0 / dtlong4

! Check whether it is time to update cumulus parameterization tendencies

  if ((istp == 1 .and. mod(time_istp8p, confrq) < dtlm) .or. &
      (istp == 1 .and. mstp == 0 .and. runtype == 'HISTADDGRID')) then

! Print message that cumulus parameterization is being computed

     write(io6, '(A,f0.2,A)') &
           ' Convective tendencies updated at ', time_istp8/3600., &
           ' hrs into simulation'

     confrq4  = real(confrq)

! Loop over all IW grid cells where cumulus parameterization may be done

!----------------------------------------------------------------------
     !$omp parallel do private(iw,mrlw,km,km1,k) &
     !$omp    schedule( guided )
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

! MRL for current IW column

        mrlw = itab_w(iw)%mrlw

! Zero out previous convective terms

        thsrc(:,iw) = 0.0
        rtsrc(:,iw) = 0.0
        qwcon(:,iw) = 0.0

        if (allocated(umsrc)) umsrc(:,iw) = 0.0
        if (allocated(vmsrc)) vmsrc(:,iw) = 0.0

        if (allocated(cu_pwa)) cu_pwa(:,iw) = 0.0
        if (allocated(cu_pev)) cu_pev(:,iw) = 0.0

        conprr(iw) =  0.0
        kcutop(iw) = -1
        kddtop(iw) = -1
        kcubot(iw) = -1
        kudbot(iw) = -1
        kddmax(iw) = -1
        kstabi(iw) = -1
        iactcu(iw) =  0
        cddf  (iw) =  0.0

! Select cumulus convection scheme based on MRL of current IW column

        if (nqparm(mrlw) == 1) then

! Tiedtke convective parameterization

           km = mza + 1 - lpw(iw) ! km  = # of T levels in cuparm_tiedtke
           km1 = km + 1           ! km1 = # of W levels in cuparm_tiedtke

           call cuparm_tiedtke(iw,km,km1,confrq4)

        elseif (nqparm(mrlw) == 2 .or. nqparm(mrlw) == 5) then

! Grell deep and shallow convection

           call grell_driver( iw, dtlong4 )

        elseif (nqparm(mrlw) == 3) then

! Kain-Fritsch deep convection

           call cuparm_kfeta(iw, dtlong4)

        elseif (nqparm(mrlw) == 4) then

! Emanuel convective parameterization

           call cuparm_emanuel(iw, dtlong4)

        endif

! Specify a minimum value of convective cloud water

        if (iactcu(iw) > 0) then
           do k = kcubot(iw), kcutop(iw)
              qwcon(k,iw) = max(qwcon(k,iw),1.e-5)
           enddo
        else
           cbmf(iw) = 0.0
           cddf(iw) = 0.0
        endif

     enddo
     !$omp end parallel do

! Print maximum heating rate in current parallel sub-domain

     dthmax = 0.
     iwqmax = 0
     kqmax  = 0

     ! this will need to be changed to work with OpenMP
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
        if (iactcu(iw) > 0) then
           do k = lpw(iw), mza
              if (thsrc(k,iw) > dthmax) then
                 dthmax = thsrc(k,iw)
                 iwqmax = iw
                 kqmax  = k
              endif
           enddo
        endif
     enddo

     write(io6, '(A,I0,A,I0,A,F0.3,A,F0.4)') " MAX CONVECTIVE HEATING RATE AT IW=",  &
           iwqmax, " K=", kqmax, " IS ", dthmax*86400., " K/DAY"

     if (iparallel == 1) then
        call mpi_send_w(i1dvara1=kcutop, i1dvara2=kcubot, i1dvara3=iactcu, &
                        r1dvara1=cbmf, r1dvara2=conprr)
     endif

     ! Update cloud fraction if convection was updated
     call calc_3d_cloud_fraction()

  endif

! Add current value of convective tendencies to thilt and rr_wt arrays
! every long timestep (whether cumulus parameterization is updated
! this timestep or not)

!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0) then

     !$omp parallel do private(iw,k,ravail,rt) schedule(guided)
     do j = 1,jtab_w(jtw_prog)%jend; iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

        ! Skip if no convection in this column

        if (iactcu(iw) < 1) cycle

        ! Skip if not doing convection on this mrl

        if (nqparm(itab_w(iw)%mrlw) == 0) then
           iactcu(iw) = 0
           conprr(iw) = 0.0
           cycle
        endif

        aconpr(iw) = aconpr(iw) + conprr(iw) * dtlong4

        do k = lpw(iw), min(kcutop(iw)+1,mza)
           if (tair(k,iw) > 253.) then
              thilt(k,iw) = thilt(k,iw) + thsrc(k,iw) * theta(k,iw) / tair(k,iw)
           else
              thilt(k,iw) = thilt(k,iw) + thsrc(k,iw) * thil(k,iw) / tair(k,iw)
           endif

           ravail = 0.95 * max(rr_w(k,iw), 0.0) * real(rho(k,iw))
           rt     = max( rtsrc(k,iw), -ravail * dtli )

           rr_wt(k,iw) = rr_wt(k,iw) + rt
        enddo

        if ( any( nqparm(itab_w(iw)%mrlw) == (/1,2,4,5/) ) ) then
           do k = lpw(iw), min(kcutop(iw)+1,mza)

              vmxet(k,iw) = vmxet(k,iw) + umsrc(k,iw) * vxn_ew(iw) &
                                        + vmsrc(k,iw) * vxn_ns(iw)

              vmyet(k,iw) = vmyet(k,iw) + umsrc(k,iw) * vyn_ew(iw) &
                                        + vmsrc(k,iw) * vyn_ns(iw)

              vmzet(k,iw) = vmzet(k,iw) + vmsrc(k,iw) * vzn_ns(iw)

           enddo
        endif

        ! Compute tracer mixing due to convection if using Grell scheme

        if ( (nl%conv_tracer_mix > 0 .and. num_cumix > 0) ) then
           if (nqparm(itab_w(iw)%mrlw) == 2) then
              call grell_mix_driver( iw, dtlong4 )
           endif
        endif

     enddo
     !$omp end parallel do

  endif

  if (iparallel == 1) then
     if ((istp == 1 .and. mod(time_istp8p, confrq) < dtlm) .or. &
         (istp == 1 .and. mstp == 0 .and. runtype == 'HISTADDGRID')) then

        call mpi_recv_w(i1dvara1=kcutop, i1dvara2=kcubot, i1dvara3=iactcu, &
                        r1dvara1=cbmf, r1dvara2=conprr)
     endif
  endif

end subroutine cuparm_driver



subroutine reset_cuparm()

  use mem_ijtabs, only: jtab_w, jtw_wadj, itab_w
  use mem_cuparm, only: conprr, kcutop, kcubot, iactcu, cbmf, cddf, kudbot, &
                        kddtop, kddmax
  use misc_coms,  only: nqparm

  implicit none

  integer :: j, iw, mrlw

  !$omp parallel do private(iw,mrlw)
  do j = 1,jtab_w(jtw_wadj)%jend; iw = jtab_w(jtw_wadj)%iw(j)

     ! MRL for current IW column
     mrlw = itab_w(iw)%mrlw

     if (nqparm(mrlw) == 0) then

        iactcu(iw) =  0

        if (allocated( conprr )) conprr(iw) =  0.0
        if (allocated( kcutop )) kcutop(iw) = -1
        if (allocated( kcubot )) kcubot(iw) = -1
        if (allocated( kudbot )) kudbot(iw) = -1
        if (allocated( kddtop )) kddtop(iw) = -1
        if (allocated( kddmax )) kddmax(iw) = -1
        if (allocated( cbmf   ))   cbmf(iw) =  0.0
        if (allocated( cddf   ))   cddf(iw) =  0.0

     endif

  enddo
  !$omp end parallel do

end subroutine reset_cuparm
