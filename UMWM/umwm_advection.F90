module umwm_advection

  implicit none

  ! Subroutines for calculation of wave energy advection in geographical (x, y)
  ! and directional (theta) space. 1st order upstream finite differences.

  real, allocatable :: cgc_dnu(:,:,:)
  real, allocatable :: rot    (:,:,:)
  real, allocatable :: dtr_max(:,:)

contains

!===============================================================================

subroutine init_propagation

  use umwm_module, only: om, pm, cg0, cthp_dirv
  use mem_sfcg,    only: mvsfc, sfcg, itab_wsfc, itab_vsfc
  use mem_para,    only: myrank
  use mem_sea,     only: omsea

  implicit none

  integer :: ivsfc, iw1, iw2, isea1, isea2, p, o
  real    :: cgv

  allocate( cgc_dnu(pm,om-2,mvsfc) )

  !$omp parallel do private(iw1,iw2,isea1,isea2,p,o,cgv)
  do ivsfc = 2, mvsfc

     cgc_dnu(:,:,ivsfc) = 0.

     iw1 = itab_vsfc(ivsfc)%iwn(1)
     if (iw1 < 2) cycle

     iw2 = itab_vsfc(ivsfc)%iwn(2)
     if (iw2 < 2) cycle

     ! restrict to V points adjacent to a primary W cell
     if (itab_wsfc(iw1)%irank /= myrank .and. itab_wsfc(iw2)%irank /= myrank) cycle

     ! restrict to V points adjacent to a sea cell
     if (sfcg%leaf_class(iw1) /= 0 .and. sfcg%leaf_class(iw2) /= 0) cycle

     isea1 = iw1 - omsea
     if (sfcg%leaf_class(iw1) /= 0) isea1 = 1

     isea2 = iw2 - omsea
     if (sfcg%leaf_class(iw2) /= 0) isea2 = 1

     do o = 1, om-2
        cgv = 0.5 * (cg0(o,isea1) + cg0(o,isea2)) * sfcg%dnu(ivsfc)
        do p = 1, pm
           cgc_dnu(p,o,ivsfc) = cgv * cthp_dirv(p,ivsfc)
        enddo
     enddo

  enddo
  !$omp end parallel do

  call check_cfl()

end subroutine init_propagation

!===============================================================================

subroutine check_cfl()

  use mem_sea,     only: msea, omsea
  use mem_sfcg,    only: sfcg, itab_wsfc
  use umwm_module, only: om, pm, dt_olam
  use mem_para,    only: myrank
  use misc_coms,   only: io6, iparallel

#ifdef OLAM_MPI
  use mpi
#endif

  implicit none

  ! 1st order upstream finite difference advection in geographical space.

  integer :: o, p, isea, iwsfc, jv, ivsfc, ierr
  real    :: dirv, cflmax
  real    :: cfl(pm,om-2)

  cflmax = 0.0

  !$omp parallel private(cfl)
  !$omp do private(isea,iwsfc,jv,ivsfc,dirv,p,o) reduction(max:cflmax)
  do isea = 2, msea
     iwsfc = isea + omsea

     ! skip cell that are not primary on this rank
     if (itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! skip cells with sufficient seaice
     if (sfcg%glatw(iwsfc) > 83.) cycle

     ! Compute cfl based on the outgoing flux from this cell
     cfl(1:pm,1:om-2) = 0.

     do jv = 1, itab_wsfc(iwsfc)%npoly
        ivsfc  = itab_wsfc(iwsfc)%ivn (jv)
        dirv   = itab_wsfc(iwsfc)%dirv(jv)

        do o = 1, om-2
           do p = 1, pm

              ! cgc: if < 0, velocity is outward from isea cell
              cfl(p,o) = cfl(p,o) - min( cgc_dnu(p,o,ivsfc) * dirv, 0.0 )

           enddo
        enddo

     enddo

     cflmax = max(cflmax, maxval(cfl) / sfcg%area(iwsfc))

  enddo
  !$omp end do nowait
  !$omp end parallel

#ifdef OLAM_MPI
  if (iparallel == 1) then
     call MPI_Allreduce(MPI_IN_PLACE, cflmax, 1, MPI_REAL, MPI_MAX, MPI_COMM_WORLD, ierr)
  endif
#endif

  write(io6,*) "UM Wave Model max propagation CFL, DT: ", cflmax * dt_olam, dt_olam

end subroutine check_cfl

!===============================================================================

subroutine propagation(isea, iwsfc)

  use mem_sea,     only: omsea
  use mem_sfcg,    only: sfcg, itab_wsfc!, nswmzons
  use umwm_module, only: evsf, evs, om, pm, oc, fice_uth, dta

  implicit none

  ! 1st order upstream finite difference advection in geographical space.

  integer, intent(in) :: isea, iwsfc
  integer             :: o, p, jv, ivsfc, iwnsfc, isean
  real                :: cdv, cvc, dirv, dtoa
  real                :: flux(pm,om-2)

  integer,  parameter :: nswmzons = 0

  ! Compute sum of fluxes through sides

  flux(:,1:oc(isea)) = 0.0

  do jv = 1, itab_wsfc(iwsfc)%npoly
     ivsfc  = itab_wsfc(iwsfc)%ivn (jv)
     iwnsfc = itab_wsfc(iwsfc)%iwn (jv)
     dirv   = itab_wsfc(iwsfc)%dirv(jv)

     isean = iwnsfc - omsea
     if (isean < 2) isean = 1  ! evs(1) = 0, radiative LBC

     cvc = 0.
     if (nswmzons > 0) cvc = sfcg%vc(ivsfc) * sfcg%dnu(ivsfc) * dirv

     ! Group velocity and current velocity normal to cell edge
     ! (if > 0, velocity is inward toward isea cell)
     do o = 1, oc(isea)
        do p = 1, pm

           cdv = cgc_dnu(p,o,ivsfc) * dirv + cvc

           flux(p,o) = flux(p,o) + max(cdv,0.) * evs(p,o,isean) &
                                 + min(cdv,0.) * evs(p,o,isea )
        enddo
     enddo
  enddo

  ! Integrate in time

  dtoa = dta / sfcg%area(iwsfc)

  do o = 1, oc(isea)
     do p = 1, pm
        evsf(p,o,isea) = evsf(p,o,isea) + dtoa * flux(p,o)
     enddo
  enddo

end subroutine propagation

!===============================================================================

subroutine init_refraction()

  use mem_sea,     only: msea, omsea
  use mem_sfcg,    only: sfcg, itab_wsfc
  use umwm_module, only: om, pm, cp0, sthp_dirv, dth
  use mem_para,    only: myrank

  implicit none

  ! 1st order upstream finite difference advection in
  ! directional space -- bottom- and current-induced refraction.

  integer :: isea, iwsfc, ivsfc, iwnsfc, iwnsea
  integer :: jv, o, p, pl
  real    :: dirv, fact, cpv
  real    :: vort(pm,om-2)

  allocate( rot (pm,om-2,msea) )
  allocate( dtr_max(om-2,msea) )

  !$omp parallel private(vort)
  !$omp do private(iwsfc,fact,jv,ivsfc,iwnsfc,dirv,iwnsea,cpv,o,p,pl)
  do isea = 2, msea
     iwsfc = isea + omsea

     ! skip border cells
     if (itab_wsfc(iwsfc)%irank /= myrank) cycle

     ! skip points close to N. Pole
     if (sfcg%glatw(iwsfc) > 83.) cycle

     vort(:,:) = 0.

     fact = 0.5 / sfcg%area(iwsfc)

     do jv = 1, itab_wsfc(iwsfc)%npoly

        ivsfc  = itab_wsfc(iwsfc)%ivn(jv)
        iwnsfc = itab_wsfc(iwsfc)%iwn(jv)
        dirv   = itab_wsfc(iwsfc)%dirv(jv) * sfcg%dnu(ivsfc)

        if (sfcg%leaf_class(iwnsfc) == 0) then
           iwnsea = iwnsfc - omsea
        else
           iwnsea = 1
        endif

        ! Vorticity from phase velocity
        do o = 1, om-2
           cpv = fact * (cp0(o,isea) + cp0(o,iwnsea)) * dirv
           do p = 1, pm
              vort(p,o) = vort(p,o) + cpv * sthp_dirv(p,ivsfc)
           enddo
        enddo

     enddo  ! jv

     ! evaluate rotation at cell edges:
     do o = 1, om-2
        do p = 1, pm
           pl = p + 1
           if (p == pm) pl = 1
           rot(p,o,isea) = 0.5 * (vort(p,o) + vort(pl,o))
        enddo
     enddo

     do o = 1, om-2
        dtr_max(o,isea) = dth / maxval( abs(rot(:,o,isea)) )
     enddo

  enddo
  !$omp end do nowait
  !$omp end parallel

end subroutine init_refraction

!===============================================================================

subroutine refraction(isea, iwsfc)

  use mem_sfcg,    only: itab_wsfc, sfcg !, nswmzons
  use umwm_module, only: om, pm, oc, evs, evsf, dth, oneovdth, dta

  implicit none

  ! 1st order upstream finite difference advection in
  ! directional space -- bottom- and current-induced refraction.

  integer, intent(in) :: isea, iwsfc
  integer             :: jv, o, p, imsfc, pl, pr
  real                :: vortvc, r
  real                :: flux(pm,om-2), dtrovdth(om-2)

  integer,  parameter :: nswmzons = 0

  ! Vorticity from swm ocean current

  if (nswmzons > 0) then

     vortvc = 0.
     do jv = 1, itab_wsfc(iwsfc)%npoly
        imsfc  = itab_wsfc(iwsfc)%imn(jv)
        vortvc = vortvc + sfcg%vort(imsfc)
     enddo
     vortvc = vortvc / real(itab_wsfc(iwsfc)%npoly)

  endif

  ! Flux between p and pl

  do o = 1, oc(isea)
     do p = 1, pm

        pl = p + 1
        if (p==pm) pl = 1

        r = rot(p,o,isea)
        if (nswmzons > 0) r = rot(p,o,isea) + vortvc

        flux(p,o) = max(r,0.) * evs(p ,o,isea) &
                  + min(r,0.) * evs(pl,o,isea)
     enddo
  enddo

  ! Make sure timestep is stable

  if (nswmzons > 0) then

     do o = 1, oc(isea)
        dtrovdth(o) = min(dta, dtr_max(o,isea) / &
                               (1. + abs(vortvc) * dtr_max(o,isea) * oneovdth)) * oneovdth
     enddo

  else

     do o = 1, oc(isea)
        dtrovdth(o) = min(dta, dtr_max(o,isea)) * oneovdth
     enddo

  endif

  ! Integrate forward in time

  do o = 1, oc(isea)
     do p = 1, pm

        pr = p - 1
        if (p==1) pr = pm

        evsf(p,o,isea) = evsf(p,o,isea) - dtrovdth(o) * (flux(p,o) - flux(pr,o))

     enddo
  enddo

end subroutine refraction

!===============================================================================

end module umwm_advection
