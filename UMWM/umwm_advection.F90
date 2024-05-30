module umwm_advection

  ! Subroutines for calculation of wave energy advection in geographical (x, y)
  ! and directional (theta) space. 1st order upstream finite differences.

  use umwm_module

  implicit none

contains

!===============================================================================

subroutine propagation()

  use mem_sea,  only: sea, msea, omsea
  use mem_sfcg, only: sfcg, itab_wsfc

  ! 1st order upstream finite difference advection in geographical space.

  integer :: o, p, isea, iwsfc, npoly, jv, ivsfc, iwnsfc, iwnsea_cg, iwnsea_e
  real :: cg, cc, cgc, dirv, areai

  real :: flux(om,pm,msea)

  !$omp parallel do private(iwsfc, npoly, jv, ivsfc, iwnsfc, dirv, &
  !$omp                     iwnsea_cg, iwnsea_e, o, p, cg, cc, cgc, areai)
  do isea = 2,msea
     iwsfc = isea + omsea
     npoly = itab_wsfc(iwsfc)%npoly
     flux(:,:,isea) = 0.

     do jv = 1,npoly
        ivsfc  = itab_wsfc(iwsfc)%ivn(jv)
        iwnsfc = itab_wsfc(iwsfc)%iwn(jv)
        dirv   = itab_wsfc(iwsfc)%dirv(jv)
        if (sfcg%leaf_class(iwnsfc) == 0) then
           iwnsea_cg = iwnsfc - omsea
           iwnsea_e  = iwnsfc - omsea
        else
           iwnsea_cg = isea
           iwnsea_e = 1 ! universal "boundary" cell, where evs = 0.
        endif

        do p = 1,pm
           do o = 1,oc(isea)

              !  group velocity and current velocity normal to cell edge
              !  (if > 0, velocity is inward toward isea cell)
              cg = 0.5 * (cg0(o,isea) + cg0(o,iwnsea_cg)) * cthp_dirv(p,ivsfc) * dirv
              cc = sfcg%vc(ivsfc)                                              * dirv
              cgc = cg + cc

              ! advective energy flux through cell edges
              flux(o,p,isea) = flux(o,p,isea) &
                             + (cgc + abs(cgc)) * sfcg%dnu(ivsfc) * evs(o,p,iwnsea_e) &
                             + (cgc - abs(cgc)) * sfcg%dnu(ivsfc) * evs(o,p,isea)

           enddo
        enddo
     enddo

     ! integrate in time
     areai = 1.0 / sfcg%area(iwsfc)
     do p = 1,pm
        do o = 1,oc(isea)
           evsf(o,p,isea) = evsf(o,p,isea) + 0.5 * dta * flux(o,p,isea) * areai

           if (sea%seaicec(isea) > fice_uth) then
              evsf(o,p,isea) = 0.0
           endif
        enddo
     enddo
  enddo
  !$omp end parallel do

end subroutine propagation

!===============================================================================

subroutine refraction()

  use mem_sea, only: sea, msea, omsea
  use mem_sfcg, only: sfcg, itab_wsfc

  ! 1st order upstream finite difference advection in
  ! directional space -- bottom- and current-induced refraction.

  integer :: isea, iwsfc, imsfc, ivsfc, iwnsfc, iwnsea
  integer :: npoly, jv, o, p
  real :: dtr_cfl, dirv, circcp, vortvc

  real :: flux(om,pm,msea)

  flux(:,:,1) = 0.

  !$omp parallel do private(iwsfc, npoly, vortvc, jv, imsfc, o, p, &
  !$omp                     circcp, ivsfc, iwnsfc, dirv, iwnsea)
  do isea = 2,msea
     iwsfc = isea + omsea
     npoly = itab_wsfc(iwsfc)%npoly
     vortvc = 0.

     ! Vorticity from swm ocean current
     do jv = 1,npoly
        imsfc  = itab_wsfc(iwsfc)%imn(jv)
        vortvc = vortvc + sfcg%vort(imsfc)
     enddo
     vortvc = vortvc / real(npoly)

     do p = 1,pm
        do o = 1,oc(isea)

           ! Vorticity from phase velocity
           circcp = 0
           do jv = 1,npoly
              ivsfc  = itab_wsfc(iwsfc)%ivn(jv)
              iwnsfc = itab_wsfc(iwsfc)%iwn(jv)
              dirv   = itab_wsfc(iwsfc)%dirv(jv)
              if (sfcg%leaf_class(iwnsfc) == 0) then
                 iwnsea = iwnsfc - omsea
              else
                 iwnsea = isea
              endif
              circcp = circcp + 0.5 * (cp0(o,isea) + cp0(o,iwnsea)) * sthp_dirv(p,ivsfc) &
                                    * sfcg%dnu(jv) * itab_wsfc(iwsfc)%dirv(jv)
           enddo
           flux(o,p,isea) = vortvc + circcp / sfcg%area(isea) ! Add both vorticity contributions

        enddo
        flux(oc(isea)+1:om,p,isea) = 0.
     enddo

     ! evaluate rotation at cell edges:
     do p = 1,pm
        do o = 1,oc(isea)
           rotl(o,p,isea) = 0.5 * (flux(o,p,isea) + flux(o,pl(p),isea))
           rotr(o,p,isea) = 0.5 * (flux(o,p,isea) + flux(o,pr(p),isea))
        enddo
     enddo
  enddo
  !$omp end parallel do

  dtr_cfl = dth / maxval(flux)

  if (dta > dtr_cfl) then
     write(6,'(a,3f10.1,e12.3)') 'dta exceeds dtr_cfl ', dta, dtr_cfl, dth, maxval(flux)
     write(6,'(a)')              'Reducing dta to dtr_cfl in subroutine refraction'
  endif

  dtr = min(dtr_cfl, dta)

  !$omp parallel do private(o, p)
  do isea = 2,msea
     do p = 1,pm
        do o = 1,oc(isea)

           ! compute tendencies
           flux(o,p,isea) = 0.5 * oneovdth                                &
              * ((rotl(o,p,isea) + abs(rotl(o,p,isea))) * evs(o,p,    isea) &
              + ( rotl(o,p,isea) - abs(rotl(o,p,isea))) * evs(o,pl(p),isea) &
              - ( rotr(o,p,isea) + abs(rotr(o,p,isea))) * evs(o,pr(p),isea) &
              - ( rotr(o,p,isea) - abs(rotr(o,p,isea))) * evs(o,p,    isea))

           ! integrate
           evsf(o,p,isea) = evsf(o,p,isea) - dtr * flux(o,p,isea)

           if (sea%seaicec(isea) > fice_uth) then
              evsf(o,p,isea) = 0.0
           endif
        enddo
     enddo
  enddo
  !$omp end parallel do

end subroutine refraction

end module umwm_advection
