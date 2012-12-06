subroutine ed_rad_wrapper(pass_flag)

  use grid_coms, only: ngrids
  use ed_state_vars, only: edgrid_g, polygontype
  use mem_leaf, only: land
  use canopy_radiation_coms, only: fvis_beam_def, fvis_diff_def, fnir_beam_def, fnir_diff_def

  implicit none

  integer :: ifm,ipy,isi,iwl
  integer, intent(in) :: pass_flag
  type(polygontype), pointer :: cpoly

  if(pass_flag == 1)then
     do ifm = 1, ngrids
        call flag_stable_cohorts(edgrid_g(ifm))
        call ed_rad_normalization(edgrid_g(ifm),1)
        call radiate_driver(edgrid_g(ifm))
        call ed2land_radiation(edgrid_g(ifm))
     enddo
  elseif(pass_flag == 2)then
     do ifm = 1, ngrids
        do ipy = 1, edgrid_g(ifm)%npolygons
           cpoly => edgrid_g(ifm)%polygon(ipy)
           iwl = edgrid_g(ifm)%iwl(ipy)
           do isi = 1, cpoly%nsites
              cpoly%met(isi)%par_beam = (land%rshort(iwl) - land%rshort_diffuse(iwl)) *   &
                   fvis_beam_def
              cpoly%met(isi)%par_diffuse = land%rshort_diffuse(iwl) * fvis_diff_def
              cpoly%met(isi)%nir_beam = (land%rshort(iwl) - land%rshort_diffuse(iwl)) *   &
                   fnir_beam_def
              cpoly%met(isi)%nir_diffuse = land%rshort_diffuse(iwl) * fnir_diff_def
              cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse +   &
                   cpoly%met(isi)%nir_diffuse
              cpoly%met(isi)%rshort = cpoly%met(isi)%rshort_diffuse +   &
                   cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam
              cpoly%met(isi)%rlong = land%rlong(iwl)

              call olam_normed_radiation(cpoly)

           enddo
           edgrid_g(ifm)%met(ipy)%par_beam = cpoly%met(1)%par_beam
           edgrid_g(ifm)%met(ipy)%par_diffuse = cpoly%met(1)%par_diffuse
           edgrid_g(ifm)%met(ipy)%nir_beam = cpoly%met(1)%nir_beam
           edgrid_g(ifm)%met(ipy)%nir_diffuse = cpoly%met(1)%nir_diffuse
        enddo
        call update_rad_avg(edgrid_g(ifm))
        call update_olam_land_rad(edgrid_g(ifm))
     enddo
  elseif(pass_flag == 3)then
     do ifm = 1, ngrids
        call int_met_avg(edgrid_g(ifm))
     enddo
  else
     print*,'Invalid pass_flag in ed_rad_wrapper',pass_flag
     stop
  endif

  return
end subroutine ed_rad_wrapper

!========================================================================

subroutine ed2land_radiation(cgrid)
  use ed_state_vars, only: edtype,polygontype,sitetype,patchtype
  use mem_leaf, only: land
  implicit none

  type(edtype), target :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype), pointer :: csite
  type(patchtype), pointer :: cpatch
  integer :: iwl,isi,ipa,ipy,ico,k
  real :: rshort_beam
  
  do ipy = 1, cgrid%npolygons
     iwl = cgrid%iwl(ipy)
     cpoly => cgrid%polygon(ipy)
     land%albedo_beam(iwl) = 0.
     land%albedo_diffuse(iwl) = 0.
     land%rlongup(iwl) = 0.
     land%rlong_albedo(iwl) = 0.
     do isi = 1, cpoly%nsites
        rshort_beam = cpoly%met(isi)%rshort - cpoly%met(isi)%rshort_diffuse
        csite => cpoly%site(isi)
        do ipa = 1, csite%npatches
           land%albedo_beam(iwl) = land%albedo_beam(iwl) + csite%albedo_beam(ipa) *   &
                csite%area(ipa) * cpoly%area(isi)
           land%albedo_diffuse(iwl) = land%albedo_diffuse(iwl) + csite%albedo_diffuse(ipa) *   &
                csite%area(ipa) * cpoly%area(isi)
           land%rlongup(iwl) = land%rlongup(iwl) + csite%rlongup(ipa) *   &
                csite%area(ipa) * cpoly%area(isi)
           land%rlong_albedo(iwl) = land%rlong_albedo(iwl) + csite%rlong_albedo(ipa) *   &
                csite%area(ipa) * cpoly%area(isi)
  
           cpatch => csite%patch(ipa)
           do ico = 1, cpatch%ncohorts
              if(cpatch%leaf_resolvable(ico) .or. cpatch%wood_resolvable(ico))then
                 cpatch%par_l_beam(ico) = cpatch%par_l_beam(ico) / cpoly%met(isi)%par_beam
                 cpatch%par_l_diffuse(ico) = cpatch%par_l_diffuse(ico) / cpoly%met(isi)%par_diffuse
                 cpatch%par_l(ico) = cpatch%par_l_beam(ico) + cpatch%par_l_diffuse(ico)

                 cpatch%rshort_l_beam(ico) = cpatch%rshort_l_beam(ico) / rshort_beam
                 cpatch%rshort_l_diffuse(ico) = cpatch%rshort_l_diffuse(ico) /   &
                      cpoly%met(isi)%rshort_diffuse
                 cpatch%rshort_l(ico) = cpatch%rshort_l_beam(ico) + cpatch%rshort_l_diffuse(ico)

                 cpatch%rshort_w_beam(ico) = cpatch%rshort_w_beam(ico) / rshort_beam
                 cpatch%rshort_w_diffuse(ico) = cpatch%rshort_w_diffuse(ico) /   &
                      cpoly%met(isi)%rshort_diffuse
                 cpatch%rshort_w(ico) = cpatch%rshort_w_beam(ico) + cpatch%rshort_w_diffuse(ico)
              endif
           enddo

           csite%par_l_beam_max(ipa) = csite%par_l_beam_max(ipa) / cpoly%met(isi)%par_beam
           csite%par_l_diffuse_max(ipa) = csite%par_l_diffuse_max(ipa) / cpoly%met(isi)%par_diffuse
           csite%par_l_max(ipa) = csite%par_l_beam_max(ipa) + csite%par_l_diffuse_max(ipa)

           csite%rshort_g_beam(ipa) = csite%rshort_g_beam(ipa) / rshort_beam
           csite%rshort_g_diffuse(ipa) = csite%rshort_g_diffuse(ipa) / cpoly%met(isi)%rshort_diffuse
           csite%rshort_g(ipa) = csite%rshort_g_beam(ipa) + csite%rshort_g_diffuse(ipa)

           do k = 1, csite%nlev_sfcwater(ipa)
              csite%rshort_s_beam(k,ipa) = csite%rshort_s_beam(k,ipa) / rshort_beam
              csite%rshort_s_diffuse(k,ipa) = csite%rshort_s_diffuse(k,ipa) /   &
                   cpoly%met(isi)%rshort_diffuse
              csite%rshort_s(k,ipa) = csite%rshort_s_beam(k,ipa) + csite%rshort_s_diffuse(k,ipa)
           enddo

      enddo
     enddo
  enddo

  

  return
end subroutine ed2land_radiation

!============================================================================
subroutine ed_rad_normalization(cgrid, pass_flag)

  use ed_state_vars, only: edtype, polygontype
  implicit none

  integer, intent(in) :: pass_flag
  type(edtype), target :: cgrid
  type(polygontype), pointer :: cpoly
  integer :: ipy, isi

  do ipy = 1,cgrid%npolygons
     cpoly => cgrid%polygon(ipy)
     do isi = 1,cpoly%nsites
        cpoly%met(isi)%par_beam = 1.0
        cpoly%met(isi)%par_diffuse = 2.0
        cpoly%met(isi)%nir_beam = 3.0
        cpoly%met(isi)%nir_diffuse = 4.0
        cpoly%met(isi)%rshort_diffuse = cpoly%met(isi)%par_diffuse +   &
             cpoly%met(isi)%nir_diffuse
        cpoly%met(isi)%rshort = cpoly%met(isi)%rshort_diffuse +   &
             cpoly%met(isi)%par_beam + cpoly%met(isi)%nir_beam
        cpoly%met(isi)%rlong = 1.0
     enddo
  enddo


  return
end subroutine ed_rad_normalization

!========================================================================

subroutine olam_normed_radiation(cpoly)
  use ed_state_vars, only: polygontype,sitetype,patchtype
  use mem_leaf, only: land
  implicit none

  type(polygontype), target :: cpoly
  type(sitetype), pointer :: csite
  type(patchtype), pointer :: cpatch
  integer :: isi,ipa,ico,k
  real :: rshort_beam
  
  do isi = 1, cpoly%nsites
     rshort_beam = cpoly%met(isi)%rshort - cpoly%met(isi)%rshort_diffuse
     csite => cpoly%site(isi)
     do ipa = 1, csite%npatches
        cpatch => csite%patch(ipa)
        do ico = 1, cpatch%ncohorts
           if(cpatch%leaf_resolvable(ico) .or. cpatch%wood_resolvable(ico))then
              cpatch%par_l_beam(ico) = cpatch%par_l_beam(ico) * cpoly%met(isi)%par_beam
              cpatch%par_l_diffuse(ico) = cpatch%par_l_diffuse(ico) * cpoly%met(isi)%par_diffuse
              cpatch%par_l(ico) = cpatch%par_l_beam(ico) + cpatch%par_l_diffuse(ico)
              
              cpatch%rshort_l_beam(ico) = cpatch%rshort_l_beam(ico) * rshort_beam
              cpatch%rshort_l_diffuse(ico) = cpatch%rshort_l_diffuse(ico) *   &
                   cpoly%met(isi)%rshort_diffuse
              cpatch%rshort_l(ico) = cpatch%rshort_l_beam(ico) + cpatch%rshort_l_diffuse(ico)

              cpatch%rlong_l_incid(ico) = cpatch%rlong_l_incid(ico) * cpoly%met(isi)%rlong
              cpatch%rlong_l(ico) = cpatch%rlong_l_incid(ico) + cpatch%rlong_l_surf(ico)
              
              cpatch%rshort_w_beam(ico) = cpatch%rshort_w_beam(ico) * rshort_beam
              cpatch%rshort_w_diffuse(ico) = cpatch%rshort_w_diffuse(ico) *   &
                   cpoly%met(isi)%rshort_diffuse
              cpatch%rshort_w(ico) = cpatch%rshort_w_beam(ico) + cpatch%rshort_w_diffuse(ico)

              cpatch%rlong_w_incid(ico) = cpatch%rlong_w_incid(ico) * cpoly%met(isi)%rlong
              cpatch%rlong_w(ico) = cpatch%rlong_w_incid(ico) + cpatch%rlong_w_surf(ico)
              
           endif
        enddo
!print*,sum(cpatch%rshort_l(1:cpatch%ncohorts)),sum(cpatch%rlong_l(1:cpatch%ncohorts)),sum(cpatch%par_l(1:cpatch%ncohorts))
        csite%par_l_beam_max(ipa) = csite%par_l_beam_max(ipa) * cpoly%met(isi)%par_beam
        csite%par_l_diffuse_max(ipa) = csite%par_l_diffuse_max(ipa) * cpoly%met(isi)%par_diffuse
        csite%par_l_max(ipa) = csite%par_l_beam_max(ipa) + csite%par_l_diffuse_max(ipa)
        
        csite%rshort_g_beam(ipa) = csite%rshort_g_beam(ipa) * rshort_beam
        csite%rshort_g_diffuse(ipa) = csite%rshort_g_diffuse(ipa) * cpoly%met(isi)%rshort_diffuse
        csite%rshort_g(ipa) = csite%rshort_g_beam(ipa) + csite%rshort_g_diffuse(ipa)

        do k = 1, csite%nlev_sfcwater(ipa)
           csite%rshort_s_beam(k,ipa) = csite%rshort_s_beam(k,ipa) * rshort_beam
           csite%rshort_s_diffuse(k,ipa) = csite%rshort_s_diffuse(k,ipa) *   &
                cpoly%met(isi)%rshort_diffuse
           csite%rshort_s(k,ipa) = csite%rshort_s_beam(k,ipa) + csite%rshort_s_diffuse(k,ipa)
        enddo

        csite%rlong_s_incid(ipa) = csite%rlong_s_incid(ipa) * cpoly%met(isi)%rlong
        csite%rlong_g_incid(ipa) = csite%rlong_g_incid(ipa) * cpoly%met(isi)%rlong

        csite%rlong_s(ipa) = csite%rlong_s_surf(ipa) + csite%rlong_s_incid(ipa)
        csite%rlong_g(ipa) = csite%rlong_g_surf(ipa) + csite%rlong_g_incid(ipa)

!print*,csite%rshort_g(ipa),csite%rlong_g(ipa)
!print*,sum(csite%rshort_s(1:csite%nlev_sfcwater(ipa),ipa)),csite%rlong_s(ipa)

!print*,cpoly%met(isi)%rshort_diffuse*csite%albedo_diffuse(ipa)+rshort_beam*csite%albedo_beam(ipa),csite%rlongup(ipa)+csite%rlong_albedo*cpoly%met(isi)%rlong

!print*,cpoly%met(isi)%rshort,cpoly%met(isi)%rlong
        
     enddo
  enddo

  return
end subroutine olam_normed_radiation

!=========================================================================

subroutine update_olam_land_rad(cgrid)
  use ed_state_vars, only: edtype,polygontype,sitetype,patchtype
  use mem_leaf, only: land
  use leaf_coms, only: nzs
  implicit none

  type(edtype), target :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype), pointer :: csite
  type(patchtype), pointer :: cpatch
  integer :: ipy, iwl, isi, ipa, k, ico
 
  do ipy = 1, cgrid%npolygons
     iwl = cgrid%iwl(ipy)
     cpoly => cgrid%polygon(ipy)

     land%rshort_s(1:nzs,iwl) = 0.
     land%rshort_g(iwl) = 0.
     land%rshort_v(iwl) = 0.
     land%rlong_g(iwl) = 0.
     land%rlong_s(iwl) = 0.
     land%rlong_v(iwl) = 0.
     land%snowfac(iwl) = 0.
     land%vf(iwl) = 0.
     land%cosz(iwl) = cgrid%cosz(ipy)

     do isi = 1, cpoly%nsites
        csite => cpoly%site(isi)
        do ipa = 1, csite%npatches
           do k = 1, csite%nlev_sfcwater(ipa)
              land%rshort_s(k,iwl) = land%rshort_s(k,iwl) + csite%rshort_s(k,ipa) *   &
                   csite%area(ipa) * cpoly%area(isi)
           enddo
           land%rshort_g(iwl) = land%rshort_g(iwl) + csite%rshort_g(ipa) *   &
                csite%area(ipa) * cpoly%area(isi)
           land%rlong_g(iwl) = land%rlong_g(iwl) + csite%rlong_g(ipa) *   &
                csite%area(ipa) * cpoly%area(isi)
           land%rlong_s(iwl) = land%rlong_s(iwl) + csite%rlong_s(ipa) *   &
                csite%area(ipa) * cpoly%area(isi)
           cpatch => csite%patch(ipa)
           do ico = 1, cpatch%ncohorts
              land%rshort_v(iwl) = land%rshort_v(iwl) + (cpatch%rshort_l(ico) +  &
                   cpatch%rshort_w(ico)) * csite%area(ipa) * cpoly%area(isi)
              land%rlong_v(iwl) = land%rlong_v(iwl) + (cpatch%rlong_l(ico) +  &
                   cpatch%rlong_w(ico)) * csite%area(ipa) * cpoly%area(isi)
           enddo
        enddo
     enddo
  enddo

  return
end subroutine update_olam_land_rad
  
!======================================================================

subroutine ed_stars_wrapper(iwl, zts, vels, rhos, airtemp, sh_vs, atmco2, &
     vkmsfc, sfluxt, sfluxr, sfluxc, ustar0, zeta, rib, ggbare)

  use grid_coms, only: ngrids
  use ed_state_vars, only: edgrid_g, polygontype, edtype, sitetype
  use mem_leaf, only: land

  implicit none

  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype), pointer :: csite
  integer, intent(in) :: iwl
  real, intent(in) :: zts, vels, rhos, airtemp, sh_vs, atmco2
  real, intent(out) :: vkmsfc, sfluxt, sfluxr, ustar0, sfluxc, &
       zeta, rib, ggbare

  integer :: ifm, ipy, isi, ipa
  real :: patch_vkmsfc, patch_sfluxt, patch_sfluxr, patch_ustar0, &
       patch_sfluxc, patch_zeta, patch_rib, patch_ggbare

  ifm = land%ed_ifm(iwl)
  ipy = land%ed_ipy(iwl)

  cgrid => edgrid_g(ifm)
  cpoly => cgrid%polygon(ipy)

  vkmsfc = 0.
  sfluxt = 0.
  sfluxr = 0.
  sfluxc = 0.
  ustar0 = 0.
  zeta = 0.
  rib = 0.
  ggbare = 0.

  do isi = 1, cpoly%nsites
     csite => cpoly%site(isi)
     do ipa = 1, csite%npatches
        call stars_co2(zts, csite%veg_rough(ipa), vels, rhos, airtemp, &
             sh_vs, atmco2, csite%can_temp(ipa), csite%can_shv(ipa), &
             csite%can_co2(ipa), patch_vkmsfc, patch_sfluxt, patch_sfluxr,   &
             patch_sfluxc, patch_ustar0, patch_zeta, patch_rib, patch_ggbare)
        vkmsfc = vkmsfc + patch_vkmsfc * csite%area(ipa) * cpoly%area(isi)
        sfluxt = sfluxt + patch_sfluxt * csite%area(ipa) * cpoly%area(isi)
        sfluxr = sfluxr + patch_sfluxr * csite%area(ipa) * cpoly%area(isi)
        sfluxc = sfluxc + patch_sfluxc * csite%area(ipa) * cpoly%area(isi)
        ustar0 = ustar0 + patch_ustar0 * csite%area(ipa) * cpoly%area(isi)
        zeta = zeta + patch_zeta * csite%area(ipa) * cpoly%area(isi)
        rib = rib + patch_rib * csite%area(ipa) * cpoly%area(isi)
        ggbare = ggbare + patch_ggbare * csite%area(ipa) * cpoly%area(isi)
     enddo
!     cpoly%met(isi)%vels = vels
     cpoly%met(isi)%geoht = zts
!     cpoly%met(isi)%vels_unstab = vels
!     cpoly%met(isi)%vels_stab = vels
  enddo

  cgrid%met(ipy)%geoht = zts

  return
end subroutine ed_stars_wrapper

!===================================================================================

subroutine copy_cuparm_to_ed()

  use ed_state_vars, only: edgrid_g, edtype, polygontype
  use leaf_coms, only: mwl
  use mem_leaf, only: land

  implicit none

  integer :: iwl, ifm, ipy, isi
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly

  do iwl = 2, mwl
     if(land%ed_flag(iwl) == 1)then
        ifm = land%ed_ifm(iwl)
        cgrid => edgrid_g(ifm)
        ipy = land%ed_ipy(iwl)
        cpoly => cgrid%polygon(ipy)
        do isi = 1, cpoly%nsites
           cpoly%met(isi)%pcpg = land%pcpg(cgrid%iwl(ipy))
           cpoly%met(isi)%qpcpg = land%qpcpg(cgrid%iwl(ipy))
           cpoly%met(isi)%dpcpg = land%dpcpg(cgrid%iwl(ipy))
        enddo
        cgrid%met(ipy)%pcpg = land%pcpg(cgrid%iwl(ipy))
        cgrid%met(ipy)%qpcpg = land%qpcpg(cgrid%iwl(ipy))
        cgrid%met(ipy)%dpcpg = land%dpcpg(cgrid%iwl(ipy))
     endif
  enddo

  return
end subroutine copy_cuparm_to_ed

!===============================================================

subroutine copy_micro_to_ed()

  use ed_state_vars, only: edgrid_g, edtype, polygontype
  use leaf_coms, only: mwl
  use mem_leaf, only: land

  implicit none

  integer :: iwl, ifm, ipy, isi
  type(edtype), pointer :: cgrid
  type(polygontype), pointer :: cpoly

  do iwl = 2, mwl
     if(land%ed_flag(iwl) == 1)then
        ifm = land%ed_ifm(iwl)
        cgrid => edgrid_g(ifm)
        ipy = land%ed_ipy(iwl)
        cpoly => cgrid%polygon(ipy)
        do isi = 1, cpoly%nsites
           cpoly%met(isi)%pcpg = land%pcpg(cgrid%iwl(ipy))
           cpoly%met(isi)%qpcpg = land%qpcpg(cgrid%iwl(ipy))
           cpoly%met(isi)%dpcpg = land%dpcpg(cgrid%iwl(ipy))
        enddo
        cgrid%met(ipy)%pcpg = land%pcpg(cgrid%iwl(ipy))
        cgrid%met(ipy)%qpcpg = land%qpcpg(cgrid%iwl(ipy))
        cgrid%met(ipy)%dpcpg = land%dpcpg(cgrid%iwl(ipy))
     endif
  enddo

  return
end subroutine copy_micro_to_ed

!===============================================================

subroutine ed_biophys_wrapper()
  use ed_misc_coms, only: integration_scheme, idoutput, imoutput, iqoutput, &
       current_time, dtlsm, frqsum, outputMonth, iyoutput, unitfast, &
       frqfast, ifoutput, itoutput, imontha, iyeara, unitstate, frqstate, &
       isoutput, outfast, outstate, nrec_fast, nrec_state
  use grid_coms, only: ngrids, time, timmax
  use ed_state_vars, only: edgrid_g, filltables, edtype, polygontype, sitetype
  use rk4_driver, only: rk4_timestep
  use ed_consts_coms, only: day_sec
  use rk4_coms, only: record_err, errmax_fout, integ_err, nerr, sanity_fout, &
       reset_integ_err
  use ed_node_coms, only: nnodetot, mynum
  use leaf_coms, only: nzg
  implicit none

  integer :: ipy, isi,ipa
  type(sitetype), pointer :: csite
  type(polygontype), pointer :: cpoly
  type(edtype), pointer :: cgrid

  character(len=32) :: fmtcntr
  integer :: ifm
  logical :: writing_dail, writing_mont, writing_dcyc, new_day
  logical, save :: past_one_day=.false.
  logical, save :: past_one_month=.false.
  logical :: new_month, new_year
  logical :: mont_analy_time, dail_analy_time, dcyc_analy_time, reset_time, &
       annual_time, writing_year, the_end, analysis_time, dcycle_time, &
       history_time
  integer :: ndays, nn
  integer, external :: num_days

  if(record_err)then
     write(fmtcntr,fmt='(a,i3.3,a)')  '(i4.4,1x,2(i3.2,1x),',nerr,'(i13,1x))'
  endif

  !----- Solve the photosynthesis and biophysics. -------------------------------------!                  
  select case (integration_scheme)
  case (0)
     do ifm=1,ngrids
        call euler_timestep(edgrid_g(ifm))
     end do
  case (1)
     do ifm=1,ngrids
        call rk4_timestep(edgrid_g(ifm),ifm)
     end do
  case (2)
     do ifm=1,ngrids
        call heun_timestep(edgrid_g(ifm))
     end do
  end select

  !------------------------------------------------------------------------------------!                  

  writing_dail      = idoutput > 0
  writing_mont      = imoutput > 0
  writing_dcyc      = iqoutput > 0
  writing_year      = iyoutput > 0
!  out_time_fast     = current_time
!  out_time_fast%month = -1
  
  !---------------------------------------------------------------------------------------!
  !     Checking if the user has indicated a need for any of the fast flux diagnostic     !
  ! variables, these are used in conditions of ifoutput,idoutput and imoutput conditions. !
  ! If they are not >0, then set the logical, fast_diagnostics to false.                  !
  !---------------------------------------------------------------------------------------!
!  fast_diagnostics = checkbudget   .or. ifoutput /= 0 .or. idoutput /= 0 .or.             &
!       imoutput /= 0 .or. iqoutput /= 0 .or. itoutput /= 0
  


  !------------------------------------------------------------------------------------!
  !     Update the daily averages if daily or monthly analysis are needed.             !
  !------------------------------------------------------------------------------------!
  if (writing_dail .or. writing_mont .or. writing_dcyc) then
     do ifm=1,ngrids
        call integrate_ed_daily_output_state(edgrid_g(ifm))
     end do
  end if
  !------------------------------------------------------------------------------------!
  
  
  
  !------------------------------------------------------------------------------------!
  !     Update the model time.                                                         !
  !------------------------------------------------------------------------------------!
  time=time+dble(dtlsm)
  call update_model_time_dm(current_time, dtlsm)
  !------------------------------------------------------------------------------------!

  !----- Check whether it is some special time... -------------------------------------!
  new_day         = current_time%time < dtlsm
  if (.not. past_one_day .and. new_day) past_one_day=.true.
  
  new_month       = current_time%date == 1  .and. new_day
  if (.not. past_one_month .and. new_month) past_one_month=.true.
  new_year        = current_time%month == 1 .and. new_month
  mont_analy_time = new_month .and. writing_mont
  dail_analy_time = new_day   .and. writing_dail
  dcyc_analy_time = new_month .and. writing_dcyc
  reset_time      = mod(time,dble(frqsum)) < dble(dtlsm)
  the_end         = mod(time,timmax) < dble(dtlsm)
  annual_time     = new_month .and. writing_year .and.                                 &
       current_time%month == outputMonth

  !----- Check whether this is time to write fast analysis output or not. -------------!
  select case (unitfast)
  case (0,1) !----- Now both are in seconds -------------------------------------------!
     analysis_time   = mod(current_time%time, frqfast) < dtlsm .and.                   &
          (ifoutput /= 0 .or. itoutput /=0)
     dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
  case (2)   !----- Months, analysis time is at the new month -------------------------!
     analysis_time   = new_month .and. (ifoutput /= 0 .or. itoutput /=0) .and.         &
          mod(real(12+current_time%month-imontha),frqfast) == 0.
     dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
  case (3) !----- Year, analysis time is at the same month as initial time ------------!
     analysis_time   = new_month .and. (ifoutput /= 0 .or. itoutput /= 0) .and.        &
          current_time%month == imontha .and.                             &
          mod(real(current_time%year-iyeara),frqfast) == 0.
     dcycle_time     = mod(current_time%time, frqfast) < dtlsm .and. iqoutput /= 0
  end select

  !----- Check whether this is time to write restart output or not. -------------------!
  select case(unitstate)
  case (0,1) !----- Now both are in seconds -------------------------------------------!
     history_time   = mod(current_time%time, frqstate) < dtlsm .and. isoutput /= 0
  case (2)   !----- Months, history time is at the new month --------------------------!
     history_time   = new_month .and. isoutput /= 0 .and.                              &
          mod(real(12+current_time%month-imontha),frqstate) == 0.
  case (3) !----- Year, history time is at the same month as initial time -------------!
     history_time   = new_month .and. isoutput /= 0 .and.                              &
          current_time%month == imontha .and.                              &
          mod(real(current_time%year-iyeara),frqstate) == 0.
  end select

  !------------------------------------------------------------------------------------!
  !    Update nrec_fast and nrec_state if it is a new month and outfast/outstate are   !
  ! monthly and frqfast/frqstate are daily or by seconds.                              !
  !------------------------------------------------------------------------------------!
  if (new_month) then
     ndays=num_days(current_time%month,current_time%year)
     if (outfast  == -2.) nrec_fast  = ndays*ceiling(day_sec/frqfast)
     if (outstate == -2.) nrec_state = ndays*ceiling(day_sec/frqstate)
  end if


  !----- Check if this is the beginning of a new simulated day. -----------------------!
  if (new_day) then
     
     if (record_err) then
        
        open (unit=77,file=trim(errmax_fout),form='formatted',access='append'          &
             ,status='old')
        write (unit=77,fmt=fmtcntr) current_time%year,current_time%month               &
             ,current_time%date,(integ_err(nn,1),nn=1,nerr)
        close(unit=77,status='keep')
        
        open (unit=78,file=trim(sanity_fout),form='formatted',access='append'          &
             ,status='old')
        write (unit=78,fmt=fmtcntr) current_time%year,current_time%month               &
             ,current_time%date,(integ_err(nn,2),nn=1,nerr)
        close(unit=78,status='keep')
        
        call reset_integ_err()
     end if
     
     !----- Do phenology, growth, mortality, recruitment, disturbance. ----------------!
     call vegetation_dynamics(new_month,new_year)
     !----- First day of a month. -----------------------------------------------------!
     if (new_month) then

            !------------------------------------------------------------------------------!
            !      On the monthly timestep we have performed various fusion/fission calls. !
            ! Therefore the var-table's pointer vectors must be updated, and the global    !
            ! definitions of the total numbers must be exported to all nodes.              !
            !      Also, if we do not need to fill the tables until we do I/O, so instead  !
            ! of running this routine every time the demographics change, we set this flag !
            ! and run the routine when the next IO occurs.                                 !
            !------------------------------------------------------------------------------!
!            if (nnodetot > 1) then
!               if (mynum == 1) write(unit=*,fmt='(a)')                                     &
!                                               '-- Monthly node synchronization - waiting'
!               if (whos_slow ) then
!                  write(unit=*,fmt='(a,1x,i5,1x,a,1x,f7.1)') 'Node',mynum                  &
!                                                            ,'time', walltime(wtime_start)
!               end if
!               call MPI_Barrier(MPI_COMM_WORLD,ierr)
!               if (mynum == 1) write(unit=*,fmt='(a)') '-- Synchronized.'
!            end if

        filltables=.true.   ! call filltab_alltypes

        !----- Re-allocate integration buffer. ----------------------------------------!
        call initialize_rk4patches(.false.)
     end if
  end if
  !------------------------------------------------------------------------------------!


  !------------------------------------------------------------------------------------!
  !     Call the model output driver.                                                  !
  !------------------------------------------------------------------------------------!
  call ed_output(analysis_time,new_day,dail_analy_time,mont_analy_time,dcyc_analy_time &
       ,annual_time,writing_dail,writing_mont,writing_dcyc,history_time       &
       ,dcycle_time,the_end)
  !------------------------------------------------------------------------------------!



  !------------------------------------------------------------------------------------!
  !      Reset time happens every frqsum.  This is to avoid variables to build up when !
  ! history and analysis are off.  This should be done outside ed_output so I have a   !
  ! chance to copy some of these to BRAMS structures.                                  !
  !------------------------------------------------------------------------------------!
  if (reset_time) then
     do ifm=1,ngrids
        call reset_averaged_vars(edgrid_g(ifm))
     end do
  end if
  !------------------------------------------------------------------------------------!

  !      Update the meteorological driver, and the hydrology parameters.               !
  !------------------------------------------------------------------------------------!
  do ifm=1,ngrids
!     call send_meteorol_to_ed(edgrid_g(ifm))
  end do
  if (new_day .and. new_month) then
     do ifm = 1,ngrids
        call updateHydroParms(edgrid_g(ifm))
     end do
  end if
  !------------------------------------------------------------------------------------!



  !------------------------------------------------------------------------------------!
  !      Update the yearly variables.                                                  !
  !------------------------------------------------------------------------------------!
  if (analysis_time .and. new_month .and. new_day .and. current_time%month == 6) then
     do ifm = 1,ngrids
        call update_ed_yearly_vars(edgrid_g(ifm))
     end do
  end if
  !------------------------------------------------------------------------------------!

  !------------------------------------------------------------------------------------!
  !      Update lateral hydrology.                                                     !
  !------------------------------------------------------------------------------------!
  call calcHydroSubsurface()
  call calcHydroSurface()
  call writeHydro()
  !------------------------------------------------------------------------------------!

  do ifm = 1,ngrids
     cgrid => edgrid_g(ifm)
     do ipy = 1, cgrid%npolygons
        cpoly => cgrid%polygon(ipy)
        do isi = 1, cpoly%nsites
           csite => cpoly%site(isi)
           do ipa = 1, csite%npatches
              call format_soil_energy_for_olam(csite,ipa,cpoly%ntext_soil(1:nzg,isi))
           enddo
        enddo
     enddo
  enddo

  return
end subroutine ed_biophys_wrapper

!===========================================================================
subroutine format_soil_energy_for_ed(csite,ipa,ntext_soil)
  use ed_state_vars, only: sitetype
  use leaf_coms, only: nzg
  use consts_coms, only: alli
  use soil_coms, only: soil
  use ed_consts_coms, only: tsupercool, allii, t3ple, cliq, cicei, cliqi, qicet3, cice
  implicit none

  type(sitetype), target :: csite
  integer, intent(in) :: ipa
  integer :: k
  real :: qwliq0
  integer, dimension(nzg), intent(in) :: ntext_soil
  real :: dryhcap, w

  do k = 1, nzg
     w = 1000.0 * csite%soil_water(k,ipa)
     dryhcap = soil(ntext_soil(k))%slcpd
     qwliq0 = w * alli
     if(csite%soil_energy(k,ipa) <= 0.)then
        csite%soil_energy(k,ipa) = csite%soil_tempk(k,ipa) *   &
             (cice * w + dryhcap)
     elseif(csite%soil_energy(k,ipa) >= qwliq0)then
        csite%soil_energy(k,ipa) = csite%soil_tempk(k,ipa) *   &
             (dryhcap + w * cliq) -  &
             w * cliq * tsupercool
     else
        csite%soil_energy(k,ipa) = csite%soil_fracliq(k,ipa) *   &
             w / allii +   &
             (dryhcap+w*cice)*t3ple
     endif
  enddo
  
  do k = 1, csite%nlev_sfcwater(ipa)
     if(csite%sfcwater_energy(k,ipa) <= 0.)then
        csite%sfcwater_energy(k,ipa) = csite%sfcwater_tempk(k,ipa) / cicei
     elseif(csite%sfcwater_energy(k,ipa) >= alli)then
        csite%sfcwater_energy(k,ipa) = (csite%sfcwater_tempk(k,ipa) - tsupercool)/cliqi
     else
        csite%sfcwater_energy(k,ipa) = (csite%sfcwater_fracliq(k,ipa) / allii) + qicet3
     endif
  enddo

  return
end subroutine format_soil_energy_for_ed

!===========================================================================

subroutine format_soil_energy_for_olam(csite,ipa,ntext_soil)
  use ed_state_vars, only: sitetype
  use leaf_coms, only: nzg
  use ed_consts_coms, only: cice, t3ple, alli, qicet3, qliqt3
  use soil_coms, only: soil
  use consts_coms, only: t00, cliq, allii

  implicit none

  type(sitetype), target :: csite
  integer, intent(in) :: ipa
  integer, dimension(nzg), intent(in) :: ntext_soil
  real :: dryhcap, w
  integer :: k
  real :: qwfroz, qwmelt

  do k = 1, nzg
     w = 1000.0 * csite%soil_water(k,ipa)
     dryhcap = soil(ntext_soil(k))%slcpd

     qwfroz = (dryhcap + w * cice) * t3ple
     qwmelt = qwfroz + w * alli

     if(csite%soil_energy(k,ipa) < qwfroz)then
        csite%soil_energy(k,ipa) = (csite%soil_tempk(k,ipa) - t00) *  &
             (2106. * w + dryhcap)
     elseif(csite%soil_energy(k,ipa) > qwmelt)then
        csite%soil_energy(k,ipa) = (csite%soil_tempk(k,ipa) - t00) *   &
             (cliq * w + dryhcap) + w * alli
     else
        csite%soil_energy(k,ipa) = csite%soil_fracliq(k,ipa) * w * alli
     endif
  enddo
  
  do k = 1, csite%nlev_sfcwater(ipa)
     if(csite%sfcwater_energy(k,ipa) <= qicet3)then
        csite%sfcwater_energy(k,ipa) = (csite%sfcwater_tempk(k,ipa) - t00) * 2106.
     elseif(csite%sfcwater_energy(k,ipa) >= qliqt3)then
        csite%sfcwater_energy(k,ipa) = (csite%sfcwater_tempk(k,ipa) - 193.36) * cliq
     else
        csite%sfcwater_energy(k,ipa) = csite%sfcwater_fracliq(k,ipa) / allii
     endif
  enddo

  return
end subroutine format_soil_energy_for_olam



!=============================================================================== 

subroutine stars_co2(zts, rough, vels, rhos, airtemp, sh_vs, atmco2,       &
                 cantemp, can_shv, can_co2, vkmsfc, sfluxt, sfluxr, sfluxc, &
                 ustar0, ed_zeta, ed_rib, ed_ggbare)

! Subroutine stars computes surface heat and vapor fluxes and momentum drag
! coefficient from Louis (1981) equations

use consts_coms, only: vonk, grav
use misc_coms,   only: io6

implicit none

! Input variables

real, intent(in) :: zts       ! height above surface of {vels, ths, sh_vs} [m]
real, intent(in) :: rough     ! surface roughness height [m]
real, intent(in) :: vels      ! atmos near-surface wind speed [m/s]
real, intent(in) :: airtemp   ! atmos near-surface temp [K]
real, intent(in) :: sh_vs     ! atmos near-surface vapor spec hum [kg_vap/m^3]
real, intent(in) :: cantemp   ! canopy air temp [K]
real, intent(in) :: can_shv   ! canopy air vapor spec hum [kg_vap/m^3]

real, intent(in) :: atmco2    ! atmos near-surface CO2 [ppm]
real, intent(in) :: can_co2   ! canopy co2 [ppm]

real(kind=8), intent(in) :: rhos  ! atmos near-surface density [kg/m^3]

! Output variables

real, intent(out) :: vkmsfc   ! surface drag coefficient for this flux cell
real, intent(out) :: sfluxt   ! surface sensible heat flux for this flux cell
real, intent(out) :: sfluxr   ! surface vapor flux for this flux cell
real, intent(out) :: sfluxc   ! surface CO2 flux for this flux cell
real, intent(out) :: ustar0   ! surface friction velocity for this flux cell

real, intent(out) :: ed_zeta
real, intent(out) :: ed_rib
real, intent(out) :: ed_ggbare

! Local parameters

real, parameter :: b = 5.
real, parameter :: csm = 7.5
real, parameter :: csh = 5.
real, parameter :: d = 5.
real, parameter :: ustmin = .1  ! lower bound on ustar (friction velocity)
real, parameter :: ubmin = .25  ! lower bound on wind speed (should use 1.0 for
                                !   convec case and 0.1 for stable case)
! Local variables

real :: vels0  ! wind speed with minimum imposed [m/s]
real :: a2     ! drag coefficient in neutral conditions, here same for h/m
real :: c1
real :: c2
real :: c3
real :: cm
real :: ch
real :: fm
real :: fh
real :: ri     ! bulk richardson numer, eq. 3.45 in Garratt
real :: tstar  !
real :: rstar  !
real :: cstar  !
real :: vtscr  ! ustar0 times density             

! Routine to compute Louis (1981) surface layer parameterization.

vels0 = max(vels,ubmin)

a2 = (vonk / log(zts / rough)) ** 2
c1 = a2 * vels0
ri = grav * zts * (airtemp - cantemp)  &
   / (.5 * (airtemp + cantemp) * vels0 * vels0)

ed_rib = ri

if (airtemp - cantemp > 0.) then   ! STABLE CASE

   fm = 1. / (1. + (2. * b * ri / sqrt(1. + d * ri)))
   fh = 1. / (1. + (3. * b * ri * sqrt(1. + d * ri)))

else                            ! UNSTABLE CASE

   c2 = b * a2 * sqrt(zts / rough * (abs(ri)))
   cm = csm * c2
   ch = csh * c2
   fm = (1. - 2. * b * ri / (1. + 2. * cm))
   fh = (1. - 3. * b * ri / (1. + 3. * ch))

endif

ustar0 = max(ustmin,sqrt(c1 * vels0 * fm))
c3 = c1 * fh / ustar0
tstar = c3 * (airtemp - cantemp)
rstar = c3 * (sh_vs - can_shv)
cstar = c3 * (atmco2 - can_co2)

ed_zeta = grav * vonk * c3 * (airtemp-cantemp)  &
     /(airtemp * ustar0 * ustar0)
ed_ggbare = c3 * ustar0

vtscr = ustar0 * rhos

vkmsfc = vtscr * ustar0 * zts / vels0
sfluxt = - vtscr * tstar
sfluxr = - vtscr * rstar
sfluxc = - vtscr * cstar

return
end subroutine stars_co2

!===============================================================================

subroutine get_ed2_atm_co2(iwl,co2)
  
  use mem_leaf,      only: land
  use ed_state_vars, only: edgrid_g
  
  integer, intent(in)  :: iwl
  real,    intent(out) :: co2

  integer :: my_ifm, my_ipy

  my_ifm = land%ed_ifm(iwl)
  my_ipy = land%ed_ipy(iwl)
  co2 = edgrid_g(my_ifm)%met(my_ipy)%atm_co2

return
end subroutine get_ed2_atm_co2

  



  
