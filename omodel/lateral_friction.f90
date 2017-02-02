subroutine lateral_friction(mrl)

  use mem_ijtabs, only: jtab_w, jtw_prog, itab_w, itabg_w
  use mem_grid,   only: mwa, mva, mza, nsw_max, dnu, dnv, dzt, zfact, &
                        arv, arw0, lpw, lsw, arw, volti, lpv
  use mem_basic,  only: vxe, vye, vze, rho
  use mem_tend,   only: vmxet, vmyet, vmzet
  use consts_coms,only: vonk

  implicit none

  integer, intent(in)        :: mrl
  integer                    :: j, iw, iv, k, ks, n, ksmax
  real                       :: awind, dist
  real,    allocatable, save :: area_log(:,:)
  real,    allocatable       :: area_sum(:)
  integer, allocatable       :: lsv(:)
  integer, allocatable, save :: lsv_max(:)
  logical,              save :: firstime = .true.
  real,    parameter         :: z0 = 8.5

  if (firstime) then
     firstime = .false.

     allocate(lsv    (mva))
     allocate(lsv_max(mwa))

     ! Find level of highest partially blocked lateral V face

     do iv = 1, mva
        do k = lpv(iv), mza-1
           if ( arv(k,iv) > 0.999 * dnu(iv) * dzt(k) * zfact(k) ) exit
        enddo
        lsv(iv) = k - 1
     enddo
     
     ! Find level of higest partially blocked face in each W column

     lsv_max = 0
     ksmax   = 1
     do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
        do n  = 1, itab_w(iw)%npoly
           iv = itab_w(iw)%iv(n)
           lsv_max(iw) = max(lsv_max(iw),lsv(iv))
        enddo
        ksmax = max(ksmax, lsv_max(iw)-lpw(iw)+1)
     enddo

     allocate(area_log(ksmax,mwa))
     allocate(area_sum(ksmax))

     area_log(:,:) = 0.0

     do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        area_sum(:)   = 0.0

        ! Sum the area of lateral faces that are closed

        do n  = 1, itab_w(iw)%npoly
           iv = itab_w(iw)%iv(n)

           do k  = lpw(iw), lsv(iv)
              ks = k - lpw(iw) + 1
              area_sum(ks) = area_sum(ks) + (dnu(iv) * dzt(k) * zfact(k) - arv(k,iv))
           enddo
        enddo

        ! Compute AREA / VOL and multiply by neutral log drag term
        ! Estimate distance to lateral obstacle as sqrt(ARW) / 2

        do k  = lpw(iw), lsv_max(iw)
           ks = k - lpw(iw) + 1

           if (k > lpw(iw)) then
              dist = 0.5 * sqrt(0.5*(arw(k,iw)+arw(k-1,iw)))
           else
              dist = 0.5 * sqrt(arw(k,iw))
           endif

           dist = max(dist, 2.0 * z0)
           
           area_log(ks,iw) = area_sum(ks) * volti(k,iw) * vonk * vonk / (log(dist/z0))**2
        enddo

     enddo

     deallocate(lsv)
     deallocate(area_sum)

  endif

  ! lateral stress flux divergence = AREA / VOL * rho * vx * |V| * (vonk / log(d/z0)) ^ 2

  if (mrl > 0) then     
                                                                                 
     !$omp parallel do private(iw,k,ks,awind)
     do j = 1, jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

        do k  = lpw(iw), lsv_max(iw)
           ks = k - lpw(iw) + 1
           
           awind = sqrt( vxe(k,iw)**2 + vye(k,iw)**2 + vze(k,iw)**2 ) * area_log(ks,iw) * rho(k,iw)

           vmxet(k,iw) = vmxet(k,iw) - awind * vxe(k,iw)
           vmyet(k,iw) = vmyet(k,iw) - awind * vye(k,iw)
           vmzet(k,iw) = vmzet(k,iw) - awind * vze(k,iw)
        enddo

     enddo

  endif

end subroutine lateral_friction
