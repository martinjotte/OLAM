!===============================================================================
! OLAM version 4.0

! Portions of this software are copied or derived from the RAMS software
! package.  The following copyright notice pertains to RAMS and its derivatives,
! including OLAM:  

   !----------------------------------------------------------------------------
   ! Copyright (C) 1991-2006  ; All Rights Reserved ; Colorado State University; 
   ! Colorado State University Research Foundation ; ATMET, LLC 

   ! This software is free software; you can redistribute it and/or modify it 
   ! under the terms of the GNU General Public License as published by the Free
   ! Software Foundation; either version 2 of the License, or (at your option)
   ! any later version. 

   ! This software is distributed in the hope that it will be useful, but
   ! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   ! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
   ! for more details.
 
   ! You should have received a copy of the GNU General Public License along
   ! with this program; if not, write to the Free Software Foundation, Inc.,
   ! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA 
   ! (http://www.gnu.org/licenses/gpl.html) 
   !----------------------------------------------------------------------------

! OLAM was developed at Duke University and the University of Miami, Florida. 
! For additional information, including published references, please contact
! the software authors, Robert L. Walko (rwalko@rsmas.miami.edu)
! or Roni Avissar (ravissar@rsmas.miami.edu).
!===============================================================================
subroutine sfcrad_ed(cosz, cp, ncohort)

  use ed_structure_defs
  use canopy_radiation_coms, only: lai_min, Watts2Ein
  use leaf_coms, only: nzg, nzs, slmsts, emisg
  use consts_coms, only: stefan

  implicit none

  integer, intent(in) :: ncohort
  type(patch) :: cp
  integer :: cohort_count
  type(cohort), pointer :: cc
  real :: fcpct
  real :: alg
  real :: als
  real :: rad
  integer :: k
  real :: fractrans
  real, dimension(nzs) :: fracabs
  real :: absg
  real :: algs
  real :: cosz
  real, dimension(ncohort) :: par_v_beam_array
  real, dimension(ncohort) :: rshort_v_beam_array
  real, dimension(ncohort) :: par_v_diffuse_array
  real, dimension(ncohort) :: rshort_v_diffuse_array
  real(kind=8), dimension(ncohort) :: veg_temp_array
  real(kind=8), dimension(ncohort) ::  lai_array
  integer, dimension(ncohort) :: pft_array
  real :: downward_par_below_beam
  real :: upward_par_above_beam
  real :: downward_nir_below_beam
  real :: upward_nir_above_beam
  integer :: il
  real :: downward_par_below_diffuse
  real :: upward_par_above_diffuse
  real :: downward_nir_below_diffuse
  real :: upward_nir_above_diffuse 
  real :: T_surface
  real :: emissivity
  real, dimension(ncohort) :: lw_v_surf_array
  real, dimension(ncohort) :: lw_v_incid_array
  real :: downward_lw_below_surf
  real :: downward_lw_below_incid
  real :: upward_lw_below_surf
  real :: upward_lw_below_incid
  real :: upward_lw_above_surf
  real :: upward_lw_above_incid
  real :: total_rshort_v
  real :: total_rlong_v
  real :: downward_rshort_below_beam
  real :: downward_rshort_below_diffuse
  real :: surface_absorbed_longwave_surf
  real :: surface_absorbed_longwave_incid


  ! Calculate the total snow depth
  if(cp%nlev_sfcwater == 0)then
     cp%total_snow_depth = 0.0
  else
     cp%total_snow_depth = sum(cp%sfcwater_depth(1:cp%nlev_sfcwater))
  endif

  ! cohort_count is the number of cohorts with leaves that are not 
  ! covered by snow.
  cohort_count = 0

  ! recalc the maximum photosynthetic rates next time around.
  cp%old_stoma_data_max(1:n_pft)%recalc = 1

  ! loop over cohorts
  cc => cp%shortest
  do while(associated(cc))

     ! initialize values
     cc%par_v_beam = 0.0
     cc%rshort_v_beam = 0.0
     cc%par_v_diffuse = 0.0
     cc%rshort_v_diffuse = 0.0
     cc%rlong_v = 0.0
     cc%old_stoma_data%recalc = 1

     ! transfer information from linked lists to arrays
     if(cc%lai > lai_min .and. cc%hite > cp%total_snow_depth)then
        cohort_count = cohort_count + 1
        pft_array(cohort_count) = cc%pft
        lai_array(cohort_count) = dble(cc%lai)
        veg_temp_array(cohort_count) = dble(cc%veg_temp)
        rshort_v_beam_array(cohort_count) = 0.0
        par_v_beam_array(cohort_count) = 0.0
        rshort_v_diffuse_array(cohort_count) = 0.0
        par_v_diffuse_array(cohort_count) = 0.0
        ! long wave doesn't need to be initialized.
     endif
     cc => cc%taller
  enddo
  cp%rshort_s_diffuse = 0.0
  cp%rshort_s_beam = 0.0

  fcpct = cp%soil_water(nzg) / slmsts(cp%ntext_soil(nzg)) ! soil water fraction

  if (fcpct > .5) then
     alg = .14                ! ground albedo
  else
     alg = .31 - .34 * fcpct  ! ground albedo
  endif

  rad = 1.0
  algs = 1.0
  if(cp%nlev_sfcwater == 0)then
     emissivity = emisg(cp%ntext_soil(nzg))
     T_surface = cp%soil_tempk(nzg)
  else

     ! Sfcwater albedo ALS ranges from wet-soil value .14 for all-liquid
     ! to .5 for all-ice
     als = 0.5 - 0.36 * cp%sfcwater_fracliq(cp%nlev_sfcwater)

     rad = 1.0 - als    ! fraction shortwave absorbed into sfcwater + soil

     do k = cp%nlev_sfcwater,1,-1

        ! fractrans is fraction of shortwave entering each sfcwater layer that
        ! gets transmitted through that layer
        fractrans = exp(-20.0 * cp%sfcwater_depth(k))

        ! fracabs(k) is fraction of total incident shortwave (at top of 
        ! top sfcwater layer) that is absorbed in each sfcwater layer
        fracabs(k) = rad * (1.0 - fractrans)

        ! rad is fraction of total incident shortwave (at top of top sfcwater 
        ! layer) that remains at bottom of current sfcwater layer
        rad = rad * fractrans
        
        ! algs will ultimately be the albedo of the soil+sfcwater.  So 
        ! subtract out whatever is getting absorbed by sfcwater.
        algs = algs - fracabs(k)
     enddo
     ! long wave parameter if sfcwater exists
     emissivity = 1.0
     T_surface = cp%sfcwater_tempk(cp%nlev_sfcwater)
  endif

  cp%snowfac = min(.99, cp%total_snow_depth / max(.001,cp%veg_height))

  ! This is the fraction of below-canopy radiation that is absorbed by 
  ! the ground
  absg = (1.0 - alg) * rad

  ! Subtract off ground absorption to obtain the soil+sfcwater albedo.
  algs = algs - absg

  ! Call the radiation parameterizations if there is vegetation
  if(cohort_count > 0)then
     
     ! Long wave first.
     call lw_twostream(cohort_count, emissivity, T_surface,  &
          pft_array(1:cohort_count), lai_array(1:cohort_count),   &
          veg_temp_array(1:cohort_count), lw_v_surf_array(1:cohort_count),  &
          lw_v_incid_array(1:cohort_count), downward_lw_below_surf,  &
          downward_lw_below_incid, upward_lw_below_surf,   &
          upward_lw_below_incid, upward_lw_above_surf, upward_lw_above_incid)
     
     ! Upwelling long wave radiation at the top of the canopy
     cp%rlongup = upward_lw_above_surf
     cp%rlong_albedo = upward_lw_above_incid
     
     ! long wave absorbed by either soil or sfcwater
     surface_absorbed_longwave_surf = downward_lw_below_surf -   &
          upward_lw_below_surf
     surface_absorbed_longwave_incid = downward_lw_below_incid -   &
          upward_lw_below_incid

     ! Compute short wave if it is daytime.
     if(cosz > 0.03)then
        
        ! call the two-stream approximation
        call sw_twostream_clump(algs,  &
             cosz,  &
             cp%lai,  &
             cohort_count,   &
             pft_array(1:cohort_count),  &
             lai_array(1:cohort_count),   &
             par_v_beam_array(1:cohort_count),  &
             par_v_diffuse_array(1:cohort_count),  &
             rshort_v_beam_array(1:cohort_count), &
             rshort_v_diffuse_array(1:cohort_count), &
             downward_par_below_beam,  &
             downward_par_below_diffuse,  &
             upward_par_above_beam,  &
             upward_par_above_diffuse,  &
             downward_nir_below_beam,  &
             downward_nir_below_diffuse,  &
             upward_nir_above_beam,   &
             upward_nir_above_diffuse)        
        
        ! below-canopy downwelling radiation
        downward_rshort_below_beam = downward_par_below_beam +   &
             downward_nir_below_beam
        downward_rshort_below_diffuse = downward_par_below_diffuse +   &
             downward_nir_below_diffuse

        ! soil+sfcwater+veg albedo (different for diffuse and beam radiation)
        cp%albedo_beam = upward_par_above_beam + upward_nir_above_beam
        cp%albedo_diffuse = upward_par_above_diffuse + upward_nir_above_diffuse
        
     else

        ! code expects values for these, even if it is not day.
        downward_rshort_below_beam = 1.0
        downward_rshort_below_diffuse = 1.0
        cp%albedo_beam = algs
        cp%albedo_diffuse = algs

     endif
     
     ! Absorption rates of PAR, rshort, and rlong of the vegetation
     cc => cp%shortest
     il = 0
     do while(associated(cc))
        if(cc%lai > lai_min .and. cc%hite > cp%total_snow_depth)then
           il = il + 1
           cc%rshort_v_beam = rshort_v_beam_array(il)
           cc%rshort_v_diffuse = rshort_v_diffuse_array(il)
           cc%par_v_beam = par_v_beam_array(il) * Watts2Ein
           cc%par_v_diffuse = par_v_diffuse_array(il) * Watts2Ein
           cc%rlong_v_surf = lw_v_surf_array(il)
           cc%rlong_v_incid = lw_v_incid_array(il)
        endif
        cc => cc%taller
     enddo

  else

     ! This is the case where there is no vegetation
     downward_rshort_below_beam = 1.0
     downward_rshort_below_diffuse = 1.0
     surface_absorbed_longwave_surf = - emissivity * stefan * T_surface**4
     surface_absorbed_longwave_incid = emissivity
     cp%albedo_beam = algs
     cp%albedo_diffuse = algs
     cp%rlongup = - surface_absorbed_longwave_surf
     cp%rlong_albedo = 1.0 - surface_absorbed_longwave_incid

  endif

  ! Absorption rate of short wave by the soil
  cp%rshort_g_beam = downward_rshort_below_beam * absg
  cp%rshort_g_diffuse = downward_rshort_below_diffuse * absg

  ! Absorption rate of short wave by the surface water
  do k=1,cp%nlev_sfcwater
     cp%rshort_s_beam(k) = downward_rshort_below_beam * fracabs(k)
     cp%rshort_s_diffuse(k) = downward_rshort_below_diffuse * fracabs(k)
  enddo

  ! Long wave absorption rate at the surface
  if(cp%nlev_sfcwater == 0)then
     cp%rlong_s_surf = 0.0
     cp%rlong_s_incid = 0.0
     cp%rlong_g_surf = surface_absorbed_longwave_surf
     cp%rlong_g_incid = surface_absorbed_longwave_incid
  else
     cp%rlong_s_surf = surface_absorbed_longwave_surf
     cp%rlong_s_incid = surface_absorbed_longwave_incid
     cp%rlong_g_surf = 0.0
     cp%rlong_g_incid = 0.0
  endif

  return
end subroutine sfcrad_ed

!================================================================

subroutine sw_twostream_clump(salb,   &
     scosz,  &
     sLAIm,   &
     ncoh,  &
     pft,  &
     LAI_in,  &
     PAR_beam_flip,  &
     PAR_diffuse_flip,  &
     SW_abs_beam_flip,   &
     SW_abs_diffuse_flip,   &
     DW_vislo_beam,  &
     DW_vislo_diffuse,  &
     UW_vishi_beam,  &
     UW_vishi_diffuse,  &
     DW_nirlo_beam,  &
     DW_nirlo_diffuse,  &
     UW_nirhi_beam,   &
     UW_nirhi_diffuse)

  use pft_coms, only: clumping_factor, n_pft, phenology
  use canopy_radiation_coms, only: diffuse_backscatter_nir,   &
       diffuse_backscatter_vis, &
       leaf_scatter_nir, leaf_scatter_vis,   &
       visible_fraction_dir, visible_fraction_dif

  implicit none

  integer :: il,ncoh,nfcoh,ipft,ncoh2,iband,i,j
  integer, dimension(2*ncoh) :: indx
  real :: salb,scosz,sLAIm,srshort,srshortd,UW_nirhi_beam,UW_nirhi_diffuse
  real :: DW_vislo_diffuse,UW_vishi_beam,DW_nirlo_beam,DW_nirlo_diffuse
  real :: DW_vislo_beam, UW_vishi_diffuse
  integer, dimension(ncoh) :: pft
  real, dimension(ncoh+1) :: SW_abs_beam_flip,PAR_beam_flip
  real, dimension(ncoh+1) :: SW_abs_diffuse_flip,PAR_diffuse_flip
  real(kind=8) :: alb,cosz,LAIm,rshort,rshortd,lambda,lambda_tot,LAI_reduction
  real(kind=8) :: beam_backscatter,eta,zeta,raddiff,diffuse_band
  real(kind=8), dimension(n_pft) :: leaf_scatter
  real(kind=8) :: exk,exki,zetai
  real(kind=8) :: d,rhoo,sigma,source_bot,source_top
  real(kind=8), dimension(n_pft) :: diffuse_backscatter
  real(kind=8) :: cumulative_lai
  
  real(kind=8), dimension(ncoh) :: LAI_in
  real(kind=8), dimension(ncoh) :: expkl_top,expkl_bot,expamk_top,expamk_bot
  real(kind=8), dimension(ncoh) :: expapk_top,expapk_bot,A_top,A_bot,B_top
  real(kind=8), dimension(ncoh) :: B_bot,C_top,C_bot,F_top,F_bot,G_top,G_bot
  real(kind=8), dimension(ncoh) :: H_top,H_bot,beam_bot
  real(kind=8), dimension(ncoh+1) :: upward_vis_beam,upward_vis_diffuse
  real(kind=8), dimension(ncoh+1) :: upward_nir_beam, upward_nir_diffuse
  real(kind=8), dimension(ncoh+1) :: downward_nir_beam, downward_nir_diffuse
  real(kind=8), dimension(ncoh+1) :: downward_vis_beam,downward_vis_diffuse
  real(kind=8), dimension(2*ncoh) :: mastervec_beam,masveccp_beam
  real(kind=8), dimension(2*ncoh) :: mastervec_diffuse,masveccp_diffuse
  real(kind=8), dimension(2*ncoh,2) :: matal
  real(kind=8), dimension(2*ncoh,5) :: mastermat
  real(kind=8), dimension(2*ncoh,2*ncoh) :: masmatcp  
  integer :: ind
  
  ! Convert input variable to double precision
  alb = dble(salb)
  cosz = dble(scosz)
  LAIm = dble(sLAIm)

  ! Calculate factors common for NIR, PAR
  ncoh2 = 2*ncoh
  lambda = 0.5/cosz
  lambda_tot = 0.0
  nfcoh = 0
  do il=1,ncoh
     nfcoh = nfcoh + 1
     ipft = pft(il)
     lambda_tot = lambda_tot + clumping_factor(ipft)
  enddo
  
  lambda_tot = lambda_tot * lambda / float(nfcoh)
  LAI_reduction = dble(min(1.0,sLAIm))
  beam_backscatter = (0.5 + cosz) * (1.0 - cosz*log(1.0+1.0/cosz))
 
  ! Loop over bands
  do iband = 1,2
     if(iband.eq.1)then
        !  First, visible (or PAR)
        do ipft = 1,n_pft
           leaf_scatter(ipft) = dble(leaf_scatter_vis(ipft))
           diffuse_backscatter(ipft) = dble(diffuse_backscatter_vis(ipft))
        enddo
!        raddiff = (rshort -rshortd) * visible_fraction_dir
!        diffuse_band = rshortd * visible_fraction_dif
     elseif(iband.eq.2)then
        !  Then, near infrared (or NIR)
        do ipft = 1,n_pft
           leaf_scatter(ipft) = dble(leaf_scatter_nir)
           diffuse_backscatter(ipft) = dble(diffuse_backscatter_nir)
        enddo
!        raddiff = (rshort -rshortd) * (1.0-visible_fraction_dir)  
!        diffuse_band = rshortd * (1.0-visible_fraction_dif) 
     endif
     
     ! Calculate more factors for this band
     
     ! Calculate the forcings
     beam_bot(ncoh) = exp(-lambda*clumping_factor(pft(ncoh))*LAI_in(ncoh))
     do il=ncoh-1,1,-1
        beam_bot(il) = beam_bot(il+1)  &
             *exp(-lambda*clumping_factor(pft(il))*LAI_in(il))
     enddo
     
     do il=1,ncoh
        ipft = pft(il)
        eta = (1.0 - (1.0-diffuse_backscatter(ipft)) * leaf_scatter(ipft))  &
             * clumping_factor(ipft)
        zeta = leaf_scatter(ipft) * diffuse_backscatter(ipft) *   &
             clumping_factor(ipft)
        exk = sqrt(eta**2 - zeta**2)
        exki = 1.0/exk
        zetai = 1.0/zeta
        ! sources
        source_bot = clumping_factor(ipft)*lambda*leaf_scatter(ipft)  &
             * beam_backscatter * beam_bot(il)

        source_top = source_bot   &
             * exp(lambda*clumping_factor(ipft)*LAI_in(il))
        ! forcing coefficients
        rhoo = - (zeta + eta + clumping_factor(ipft)*lambda)   &
             * clumping_factor(ipft) * lambda    &
             * leaf_scatter(ipft) * beam_backscatter   &
             * beam_bot(il)
        sigma = clumping_factor(ipft) * lambda
        ! calculate exponentials only once
        expkl_bot(il)=1.0
        expkl_top(il)=exp(exk*LAI_in(il))
        expamk_bot(il) = 1.0
        expamk_top(il) = exp((sigma-exk)*LAI_in(il))
        expapk_bot(il) = 1.0
        expapk_top(il) = exp((sigma+exk)*LAI_in(il))
        A_bot(il) = -source_bot*zetai
        A_top(il) = -source_top*zetai &
             +0.5*zetai*(eta*exki-1.0)*expkl_top(il)*rhoo/(sigma-exk) &
             *(expamk_top(il)-expamk_bot(il))  &
             -0.5*zetai*(eta*exki+1.0)/expkl_top(il)*rhoo/(sigma+exk) &
             *(expapk_top(il)-expapk_bot(il))
        B_bot(il) = 0.5*zetai*(eta*exki-1.0)
        B_top(il) = 0.5*zetai*(eta*exki-1.0)*expkl_top(il)
        C_bot(il) = -0.5*zetai*(eta*exki+1.0)
        C_top(il) = -0.5*zetai*(eta*exki+1.0)/expkl_top(il)
        F_bot(il) = 0.0
        F_top(il) = 0.5*exki*expkl_top(il)*rhoo/(sigma-exk)  &
             *(expamk_top(il)-expamk_bot(il))  &
             -0.5*exki/expkl_top(il)*rhoo/(sigma+exk)  &
             *(expapk_top(il)-expapk_bot(il))
        G_bot(il) = 0.5*exki
        G_top(il) = 0.5*exki*expkl_top(il)
        H_bot(il) = -0.5*exki
        H_top(il) = -0.5*exki/expkl_top(il)
     enddo

     
     ! Organize the matrix coefficients
     
     do j=1,ncoh2
        do i=1,ncoh2
           masmatcp(i,j)=0.0
        enddo
        mastervec_beam(j)=0.0
        masveccp_beam(j)=0.0
        mastervec_diffuse(j)=0.0
        masveccp_diffuse(j)=0.0
     enddo
     
     masmatcp(1,1)=G_top(ncoh)
     masmatcp(1,2)=H_top(ncoh)
     mastervec_beam(1)=-F_top(ncoh)
     mastervec_diffuse(1)=1.0
     masveccp_beam(1)=mastervec_beam(1)
     masveccp_diffuse(1)=mastervec_diffuse(1)

     do i=2,ncoh2-2,2
        masmatcp(i,i-1)=G_bot(nint((ncoh2-i+2)*0.5))
        masmatcp(i,i)=H_bot(nint((ncoh2-i+2)*0.5))
        masmatcp(i,i+1)=-G_top(nint((ncoh2-i)*0.5))
        masmatcp(i,i+2)=-H_top(nint((ncoh2-i)*0.5))
        mastervec_beam(i)=-F_bot(nint((ncoh2-i+2)*0.5))+F_top(nint((ncoh2-i)*0.5))
        mastervec_diffuse(i)=0.0
        masveccp_beam(i)=mastervec_beam(i)
        masveccp_diffuse(i)=mastervec_diffuse(i)
     enddo
     do i=3,ncoh2-1,2
        masmatcp(i,i-2)=B_bot(nint((ncoh2-i+3)*0.5))
        masmatcp(i,i-1)=C_bot(nint((ncoh2-i+3)*0.5))
        masmatcp(i,i)=-B_top(nint((ncoh2-i+1)*0.5))
        masmatcp(i,i+1)=-C_top(nint((ncoh2-i+1)*0.5))
        mastervec_beam(i)=-A_bot(nint((ncoh2-i+3)*0.5))+A_top(nint((ncoh2-i+1)*0.5))
        masveccp_beam(i)=mastervec_beam(i)
        mastervec_diffuse(i)=0.0
        masveccp_diffuse(i)=mastervec_diffuse(i)
     enddo
     masmatcp(ncoh2,ncoh2-1)=B_bot(1)-alb*G_bot(1)
     masmatcp(ncoh2,ncoh2)=C_bot(1)-alb*H_bot(1)
     mastervec_beam(ncoh2)= -A_bot(1)+alb*beam_bot(1)
     masveccp_beam(ncoh2)= mastervec_beam(ncoh2)
     mastervec_diffuse(ncoh2)= 0.0
     masveccp_diffuse(ncoh2)= mastervec_diffuse(ncoh2)
     
     ! Prep for inversion
     
     mastermat(1,1)=0.
     mastermat(1,2)=0.
     mastermat(1,3)=masmatcp(1,1)
     mastermat(1,4)=masmatcp(1,2)
     mastermat(1,5)=0.
     do i=2,ncoh2-2,2
        mastermat(i,1)=0.
        mastermat(i,2)=masmatcp(i,i-1)
        mastermat(i,3)=masmatcp(i,i)
        mastermat(i,4)=masmatcp(i,i+1)
        mastermat(i,5)=masmatcp(i,i+2)
     enddo
     do i=3,ncoh2-1,2
        mastermat(i,1)=masmatcp(i,i-2)
        mastermat(i,2)=masmatcp(i,i-1)
        mastermat(i,3)=masmatcp(i,i)
        mastermat(i,4)=masmatcp(i,i+1)
        mastermat(i,5)=0.
     enddo
     mastermat(ncoh2,1)=0.
     mastermat(ncoh2,2)=masmatcp(ncoh2,ncoh2-1)
     mastermat(ncoh2,3)=masmatcp(ncoh2,ncoh2)
     mastermat(ncoh2,4)=0.
     mastermat(ncoh2,5)=0.
     
     ! Invert the matrix
     call bandec(mastermat,ncoh2,2,2,matal,indx,d)

     ! Backsubstitute for beam and diffuse
     call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_beam)
     call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_diffuse)
     
     ! Improve the solution
     call mprove(masmatcp,mastermat,matal,ncoh2,5,2,indx,masveccp_beam,mastervec_beam)
     call mprove(masmatcp,mastermat,matal,ncoh2,5,2,indx,masveccp_diffuse,mastervec_diffuse)
     

     if(iband.eq.1)then
        do i=3,ncoh2-1,2
           ind = nint((ncoh2-i+1)*0.5)+1
           upward_vis_beam(ind) = A_bot(ind) + B_bot(ind) *  &
                mastervec_beam(i-2) + C_bot(ind) * mastervec_beam(i-1)
           upward_vis_diffuse(ind) = B_bot(ind) * mastervec_diffuse(i-2)   &
                + C_bot(ind) * mastervec_diffuse(i-1)
        enddo
        do i=2,ncoh2-2,2
           ind = nint((ncoh2-i)*0.5)+1
           downward_vis_beam(ind) = beam_bot(ind) + F_bot(ind)  &
                + H_bot(ind) * mastervec_beam(i)  &
                + G_bot(ind) * mastervec_beam(i-1)
           downward_vis_diffuse(ind) = H_bot(ind) * mastervec_diffuse(i)  &
                + G_bot(ind) * mastervec_diffuse(i-1)
        enddo
        upward_vis_beam(ncoh+1) = B_top(ncoh) * mastervec_beam(1)  &
             + C_top(ncoh) * mastervec_beam(2) + A_top(ncoh)
        upward_vis_diffuse(ncoh+1) = B_top(ncoh) * mastervec_diffuse(1)  &
             + C_top(ncoh) * mastervec_diffuse(2)
        downward_vis_beam(ncoh+1) = 1.0
        downward_vis_diffuse(ncoh+1) = 1.0
        downward_vis_beam(1) = G_bot(1) * mastervec_beam(ncoh2-1) +   &
             H_bot(1) * mastervec_beam(ncoh2) + F_bot(1) + beam_bot(1)
        downward_vis_diffuse(1) = G_bot(1) * mastervec_diffuse(ncoh2-1) +  &
             H_bot(1) * mastervec_diffuse(ncoh2)
        upward_vis_beam(1) = alb * downward_vis_beam(1)
        upward_vis_diffuse(1) = alb * downward_vis_diffuse(1)
     elseif(iband.eq.2)then
        do i=3,ncoh2-1,2
           ind = nint((ncoh2-i+1)*0.5)+1
           upward_nir_beam(ind) = A_bot(ind) + B_bot(ind) *   &
                mastervec_beam(i-2) + C_bot(ind) * mastervec_beam(i-1)
           upward_nir_diffuse(ind) = B_bot(ind) * mastervec_diffuse(i-2) +  &
                C_bot(ind) * mastervec_diffuse(i-1)
        enddo
        do i=2,ncoh2-2,2
           ind = nint((ncoh2-i)*0.5)+1
           downward_nir_beam(ind) = beam_bot(ind) + F_bot(ind) +   &
                H_bot(ind) * mastervec_beam(i) + G_bot(ind) *   &
                mastervec_beam(i-1)
           downward_nir_diffuse(ind) = H_bot(ind) * mastervec_diffuse(i)  +  &
                G_bot(ind) * mastervec_diffuse(i-1)
        enddo
        upward_nir_beam(ncoh+1) = B_top(ncoh) * mastervec_beam(1) +  &
             C_top(ncoh) * mastervec_beam(2) + A_top(ncoh)
        upward_nir_diffuse(ncoh+1) = B_top(ncoh) * mastervec_diffuse(1) +  &
             C_top(ncoh) * mastervec_diffuse(2)
        downward_nir_beam(ncoh+1) = 1.0
        downward_nir_diffuse(ncoh+1) = 1.0
        downward_nir_beam(1) = G_bot(1) * mastervec_beam(ncoh2-1) +  &
             H_bot(1) * mastervec_beam(ncoh2) + F_bot(1) + beam_bot(1)
        downward_nir_diffuse(1) = G_bot(1) * mastervec_diffuse(ncoh2-1) +  &
             H_bot(1) * mastervec_diffuse(ncoh2)
        upward_nir_beam(1) = alb*downward_nir_beam(1)
        upward_nir_diffuse(1) = alb*downward_nir_diffuse(1)
     endif
  enddo
  
  do il=1,ncoh
     PAR_beam_flip(il) = visible_fraction_dir * sngl(    &
          downward_vis_beam(il+1) - downward_vis_beam(il)  +   &
          upward_vis_beam(il) - upward_vis_beam(il+1)  )
     PAR_diffuse_flip(il) = visible_fraction_dif * sngl(    &
          downward_vis_diffuse(il+1) - downward_vis_diffuse(il)  +   &
          upward_vis_diffuse(il) - upward_vis_diffuse(il+1)  )
     SW_abs_beam_flip(il) = PAR_beam_flip(il) +   &
          (1.0 - visible_fraction_dir) * sngl(  &
          downward_nir_beam(il+1)-downward_nir_beam(il)  &
          +upward_nir_beam(il)-upward_nir_beam(il+1))
     SW_abs_diffuse_flip(il) = PAR_diffuse_flip(il) +   &
          (1.0 - visible_fraction_dif) * sngl(  &
          downward_nir_diffuse(il+1)-downward_nir_diffuse(il)  &
          +upward_nir_diffuse(il)-upward_nir_diffuse(il+1))
     if(PAR_beam_flip(il).lt.0.0)PAR_beam_flip(il)=0.0
     if(PAR_diffuse_flip(il).lt.0.0)PAR_diffuse_flip(il)=0.0
     if(SW_abs_beam_flip(il).lt.0.0)SW_abs_beam_flip(il)=0.0
     if(SW_abs_diffuse_flip(il).lt.0.0)SW_abs_diffuse_flip(il)=0.0
  enddo
  
  DW_vislo_beam = max(0.0,sngl(downward_vis_beam(1))) * visible_fraction_dir
  DW_vislo_diffuse = max(0.0,sngl(downward_vis_diffuse(1))) * visible_fraction_dif
  UW_vishi_beam = max(0.0,sngl(upward_vis_beam(ncoh+1))) * visible_fraction_dir
  UW_vishi_diffuse = max(0.0,sngl(upward_vis_diffuse(ncoh+1))) *visible_fraction_dif
  DW_nirlo_beam = max(0.0,sngl(downward_nir_beam(1))) * (1.0-visible_fraction_dir)
  DW_nirlo_diffuse = max(0.0,sngl(downward_nir_diffuse(1))) * (1.0-visible_fraction_dif)
  UW_nirhi_beam = max(0.0,sngl(upward_nir_beam(ncoh+1))) * (1.0-visible_fraction_dir)
  UW_nirhi_diffuse = max(0.0,sngl(upward_nir_diffuse(ncoh+1))) * (1.0-visible_fraction_dif)
  
  return
end subroutine sw_twostream_clump

!==============================================================
subroutine bandec(a,n,m1,m2,al,indx,d)
implicit none

real(kind=8), parameter :: tiny=1.0e-20
integer :: n,m1,m2
real(kind=8), dimension(n,m1+m2+1) :: a
real(kind=8), dimension(n,m1) :: al
integer, dimension(n) :: indx
real(kind=8) :: d,tvar
integer :: i,j,k,l,mm
real(kind=8) :: dum

!------------------

mm=m1+m2+1
l=m1
do i=1,m1
   do j=m1+2-i,mm
      a(i,j-l)=a(i,j)
   enddo
   l=l-1
   do j=mm-l,mm
      a(i,j)=0.
   enddo
enddo
d=1.0
l=m1
do k=1,n
   dum=a(k,1)
   i=k
   if(l.lt.n)l=l+1
   do j=k+1,l
      if(abs(a(j,1)).gt.abs(dum))then
         dum=a(j,1)
         i=j
      endif
   enddo
   indx(k)=i
   if(dum.eq.0.0)a(k,1)=tiny
   if(i.ne.k)then
      d=-d
      do j=1,mm
         tvar=a(k,j)
         a(k,j)=a(i,j)
         a(i,j)=tvar
      enddo
   endif
   do i=k+1,l
      dum=a(i,1)/a(k,1)
      al(k,i-k)=dum
      do j=2,mm
         a(i,j-1)=a(i,j)-dum*a(k,j)
      enddo
      a(i,mm)=0.0
   enddo
enddo
return
end subroutine bandec


!==========================================================================
subroutine banbks(a,n,m1,m2,al,indx,b)
implicit none

integer :: n,m1,m2,i,k,l,mm
real(kind=8), dimension(n,m1+m2+1) :: a
real(kind=8), dimension(n,m1) :: al
integer, dimension(n) :: indx
real(kind=8), dimension(n) :: b
real(kind=8) :: dum,tvar
!-------------------------

mm=m1+m2+1
l=m1
do k=1,n
   i=indx(k)
   if(i.ne.k)then
      tvar=b(k)
      b(k)=b(i)
      b(i)=tvar
   endif
   if(l.lt.n)l=l+1
   do i=k+1,l
      b(i)=b(i)-al(k,i-k)*b(k)
   enddo
enddo
l=1




do i=n,1,-1
   dum=b(i)
   do k=2,l
      dum=dum-a(i,k)*b(k+i-1)
   enddo
   b(i)=dum/a(i,1)
   if(l.lt.mm)l=l+1
enddo



return
end subroutine banbks

!=======================================================================
subroutine mprove(a,alud,matal,n,np,npp,indx,b,x)
implicit none

integer, parameter :: nmax=100
integer :: n,np,i,j,npp
real(kind=8), dimension(n,n) :: a
real(kind=8), dimension(n,np) :: alud
real(kind=8), dimension(n,npp) :: matal
integer, dimension(n) :: indx
real(kind=8), dimension(n) :: b,x
real(kind=8), dimension(n) :: r
real(kind=8) :: sdp

!-------------------------------------
do i=1,n
   sdp=-b(i)
   do j=1,n
      sdp=sdp+a(i,j)*x(j)
   enddo
   r(i)=sdp
enddo

call banbks(alud,n,npp,npp,matal,indx,r)
do i=1,n
   x(i)=x(i)-r(i)
enddo
return
end subroutine mprove

!=====================================================================
subroutine lw_twostream(ncoh, semgs, sT_grnd, pft, LAI, T_veg,   &
     lw_v_surf, lw_v_incid, downward_lw_below_surf, downward_lw_below_incid,  &
     upward_lw_below_surf, upward_lw_below_incid, upward_lw_above_surf,  &
     upward_lw_above_incid)

  use canopy_radiation_coms, only: emis_v, mubar
  use consts_coms, only: stefan

  implicit none
  
  integer :: ncoh
  real :: semgs
  real :: st_grnd
  integer, dimension(ncoh) :: pft
  real(kind=8), dimension(ncoh) :: LAI
  real(kind=8), dimension(ncoh) :: T_veg
  real, dimension(ncoh) :: lw_v_surf
  real, dimension(ncoh) :: lw_v_incid
  real :: downward_lw_below_surf
  real :: upward_lw_above_surf
  real :: upward_lw_below_surf
  real :: downward_lw_below_incid
  real :: upward_lw_above_incid
  real ::  upward_lw_below_incid
  real(kind=8) :: emgs
  real(kind=8) :: T_grnd
  integer :: ncoh2
  integer :: il
  real(kind=8) :: zeta
  real(kind=8) :: eta
  real(kind=8) :: exk
  real(kind=8) :: zetai
  real(kind=8) :: exki
  real(kind=8), dimension(ncoh) :: source
  real(kind=8), dimension(ncoh) :: forcing
  real(kind=8), dimension(ncoh+1) :: explai
  real(kind=8), dimension(ncoh+1) :: exmlai
  real(kind=8), dimension(ncoh) :: A_dw
  real(kind=8), dimension(ncoh) :: B_dw
  real(kind=8), dimension(ncoh) :: C_dw
  real(kind=8), dimension(ncoh) :: D_dw
  real(kind=8), dimension(ncoh) :: A_uw
  real(kind=8), dimension(ncoh) :: B_uw
  real(kind=8), dimension(ncoh) :: C_uw
  real(kind=8), dimension(ncoh) :: D_uw
  real(kind=8), dimension(ncoh) :: E_uw
  real(kind=8), dimension(ncoh) :: F_uw
  integer :: i
  integer :: j
  real(kind=8), dimension(2*ncoh,2*ncoh) :: masmatcp  
  real(kind=8), dimension(2*ncoh) :: mastervec_surf
  real(kind=8), dimension(2*ncoh) :: mastervec_incid
  real(kind=8), dimension(ncoh+1) :: DW_incid
  real(kind=8), dimension(ncoh+1) :: DW_surf
  real(kind=8), dimension(ncoh+1) :: UW_incid
  real(kind=8), dimension(ncoh+1) :: UW_surf
  real(kind=8), dimension(2*ncoh,5) :: mastermat
  real(kind=8), dimension(2*ncoh,2) :: matal
  integer :: ind
  integer, dimension(2*ncoh) :: indx
  real(kind=8) :: d

  emgs = dble(semgs)
  t_grnd = dble(st_grnd)

  ncoh2 = 2*ncoh

  do il=1,ncoh
     zeta=2.0*(1.0-emis_v(pft(il)))/(3.0*mubar)
     eta=(2.0+emis_v(pft(il)))/(3.0*mubar)
     exk=sqrt(eta**2-zeta**2)
     exki=1.0/exk
     zetai=1.0/zeta
     source(il) = emis_v(pft(il))*stefan*T_veg(il)**4
     forcing(il) = -(zeta+eta)*source(il)
     explai(il)=exp(exk*LAI(il))
     exmlai(il)=exp(-exk*LAI(il))

     ! coefficient of lambda1 (and minus the coefficient of lambda2) for the bottom of a layer, downwelling radiation.
     A_dw(il) = 0.5 * exki

     ! coefficient of lambda1, top of layer, downwelling radiation
     B_dw(il) = 0.5*exki*explai(il)

     ! coefficient of lambda2, top of layer, downwelling radiation
     C_dw(il) = -0.5*exki*exmlai(il)

     ! term of downwelling radiation not multiplying a lambda
     D_dw(il) = 0.5*(exki**2)*forcing(il) * (explai(il) + exmlai(il) - 2.0)

     A_uw(il) = 0.5*zetai*(eta*exki-1.0)
     B_uw(il) = -0.5*zetai*(eta*exki+1.0)
     C_uw(il) =  -source(il)*zetai
     D_uw(il) = A_uw(il) * explai(il)
     E_uw(il) = B_uw(il) * exmlai(il)
     F_uw(il) = -source(il)*zetai  &
          +0.5*zetai*(eta*exki-1.0)*explai(il)  &
          *(forcing(il)*exki*(1.0           &
          -exmlai(il)))  &
          -0.5*zetai*(eta*exki+1.0)*exmlai(il)  &
          *(forcing(il)*exki*(explai(il)-1.0))
  enddo

  do j=1,ncoh2
     do i=1,ncoh2
        masmatcp(i,j)=0.0
     enddo
     mastervec_surf(j)=0.0
     mastervec_incid(j)=0.0
  enddo

  ! Vector is of the form: (lambda_N, lambda_{N-1},...,lambda_1)

  masmatcp(1,1)=B_dw(ncoh)
  masmatcp(1,2)=C_dw(ncoh)
  mastervec_surf(1)=-D_dw(ncoh)
  mastervec_incid(1)=1.0

  do i=2,ncoh2-2,2
     ind = nint((ncoh2-i)*0.5)
     masmatcp(i,i-1)=-A_dw(ind+1)
     masmatcp(i,i)=A_dw(ind+1)
     masmatcp(i,i+1)=B_dw(ind)
     masmatcp(i,i+2)=C_dw(ind)
     mastervec_surf(i)=-D_dw(ind)
     mastervec_incid(i)=0.0
  enddo
  do i=3,ncoh2-1,2
     ind = nint((ncoh2-i+1)*0.5)
     masmatcp(i,i-2)=-A_uw(ind+1)
     masmatcp(i,i-1)=-B_uw(ind+1)
     masmatcp(i,i)=D_uw(ind)
     masmatcp(i,i+1)=E_uw(ind)
     mastervec_surf(i)=C_uw(ind+1)-F_uw(ind)
     mastervec_incid(i)=0.0
  enddo
  masmatcp(ncoh2,ncoh2-1)=A_uw(1)-(1.0-emgs)*A_dw(1)
  masmatcp(ncoh2,ncoh2)=B_uw(1)+(1.0-emgs)*A_dw(1)
  mastervec_surf(ncoh2)=emgs*stefan*T_grnd**4 - C_uw(1)
  mastervec_incid(ncoh2)=0.0

  mastermat(1,1)=0.
  mastermat(1,2)=0.
  mastermat(1,3)=masmatcp(1,1)
  mastermat(1,4)=masmatcp(1,2)
  mastermat(1,5)=0.
  do i=2,ncoh2-2,2
     mastermat(i,1)=0.
     mastermat(i,2)=masmatcp(i,i-1)
     mastermat(i,3)=masmatcp(i,i)
     mastermat(i,4)=masmatcp(i,i+1)
     mastermat(i,5)=masmatcp(i,i+2)
  enddo
  do i=3,ncoh2-1,2
     mastermat(i,1)=masmatcp(i,i-2)
     mastermat(i,2)=masmatcp(i,i-1)
     mastermat(i,3)=masmatcp(i,i)
     mastermat(i,4)=masmatcp(i,i+1)
     mastermat(i,5)=0.
  enddo
  mastermat(ncoh2,1)=0.
  mastermat(ncoh2,2)=masmatcp(ncoh2,ncoh2-1)
  mastermat(ncoh2,3)=masmatcp(ncoh2,ncoh2)
  mastermat(ncoh2,4)=0.
  mastermat(ncoh2,5)=0.
  
  ! Invert matrix
  call bandec(mastermat,ncoh2,2,2,matal,indx,d)

  ! Backsubstitute for contributions of ground and vegetation
  call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_surf)

  ! Backsubstitute for contribution of incident longwave at canopy top
  call banbks(mastermat,ncoh2,2,2,matal,indx,mastervec_incid)

  do i=3,ncoh2-1,2
     ind = nint((ncoh2-i+1)*0.5)
     UW_surf(ind+1)=masmatcp(i,i)*mastervec_surf(i)   &
          +masmatcp(i,i+1)*mastervec_surf(i+1)+F_uw(ind)
     UW_incid(ind+1)=masmatcp(i,i)*mastervec_incid(i)   &
          +masmatcp(i,i+1)*mastervec_incid(i+1)
  enddo
  do i=2,ncoh2-2,2
     ind = nint((ncoh2-i)*0.5)
     DW_surf(ind+1)=masmatcp(i,i+1)*mastervec_surf(i+1)   &
          +masmatcp(i,i+2)*mastervec_surf(i+2)+D_dw(ind)
     DW_incid(ind+1)=masmatcp(i,i+1)*mastervec_incid(i+1)   &
          +masmatcp(i,i+2)*mastervec_incid(i+2)
  enddo

  UW_surf(ncoh+1)=D_uw(ncoh)*mastervec_surf(1)+E_uw(ncoh)*mastervec_surf(2)  &
       + F_uw(ncoh)
  UW_incid(ncoh+1)=D_uw(ncoh)*mastervec_incid(1)+E_uw(ncoh)*mastervec_incid(2)
  DW_surf(ncoh+1)=0.0
  DW_incid(ncoh+1)=1.0
  DW_surf(1)=A_dw(1)*(mastervec_surf(ncoh2-1)-mastervec_surf(ncoh2))
  DW_incid(1)=A_dw(1)*(mastervec_incid(ncoh2-1)-mastervec_incid(ncoh2))
  UW_surf(1)=(1.0-emgs)*DW_surf(1)+emgs*stefan*T_grnd**4
  UW_incid(1)=(1.0-emgs)*DW_incid(1)

  do il = 1,ncoh
     lw_v_surf(il) = sngl(DW_surf(il+1) - DW_surf(il) + UW_surf(il) -   &
          UW_surf(il+1))
     lw_v_incid(il) = sngl(DW_incid(il+1) - DW_incid(il) + UW_incid(il) -   &
          UW_incid(il+1))
  enddo

  downward_lw_below_surf = sngl(DW_surf(1))
  downward_lw_below_incid = sngl(DW_incid(1))
  upward_lw_below_surf = sngl(UW_surf(1))
  upward_lw_below_incid = sngl(UW_incid(1))
  upward_lw_above_surf = sngl(UW_surf(ncoh+1))
  upward_lw_above_incid = sngl(UW_incid(ncoh+1))

  return
end subroutine lw_twostream

!===================================================================

subroutine ed2land_radiation(cs)
  use ed_structure_defs
  use mem_leaf, only: land

  implicit none

  type(site)           :: cs
  type(patch), pointer :: cp

  ! initialize arrays to zero
  land%albedo_beam(cs%iland) = 0.0
  land%albedo_diffuse(cs%iland) = 0.0
  land%rlongup(cs%iland) = 0.0
  land%rlong_albedo(cs%iland) = 0.0

  ! loop over patches
  cp => cs%oldest_patch
  do while(associated(cp))

     ! compute cell-level albedo and upward longwave, weighting patches by fractional area
!     land%albedo(cs%iland) = land%albedo(cs%iland) + cp%albedo * cp% area
     land%albedo_beam(cs%iland) = land%albedo_beam(cs%iland) + cp%albedo_beam * cp%area
     land%albedo_diffuse(cs%iland) = land%albedo_diffuse(cs%iland) + cp%albedo_diffuse * cp%area
     land%rlongup(cs%iland) = land%rlongup(cs%iland) + cp%rlongup * cp% area
     land%rlong_albedo(cs%iland) = land%rlong_albedo(cs%iland) + cp%rlong_albedo * cp% area

     cp => cp%younger
  enddo

  return
end subroutine ed2land_radiation

!==========================================================================

subroutine scale_ed_radiation(cp)
  use ed_structure_defs
  use mem_leaf, only: land
  use canopy_radiation_coms, only: lai_min

  implicit none

  type(patch)           :: cp
  type(cohort), pointer :: cc
  real :: beam_radiation
  integer :: k

  beam_radiation = land%rshort(cp%siteptr%iland) -   &
       land%rshort_diffuse(cp%siteptr%iland)
  cc => cp%shortest
  do while(associated(cc))
     if(cc%lai > lai_min .and. cc%hite > cp%total_snow_depth)then
        cc%rshort_v_beam = cc%rshort_v_beam * beam_radiation
        cc%rshort_v_diffuse = cc%rshort_v_diffuse *   &
             land%rshort_diffuse(cp%siteptr%iland)
        cc%rshort_v = cc%rshort_v_beam + cc%rshort_v_diffuse

        cc%par_v_beam = cc%par_v_beam * beam_radiation
        cc%par_v_diffuse = cc%par_v_diffuse *   &
             land%rshort_diffuse(cp%siteptr%iland)
        cc%par_v = cc%par_v_beam + cc%par_v_diffuse

        cc%rlong_v_incid = cc%rlong_v_incid * land%rlong(cp%siteptr%iland)
        cc%rlong_v = cc%rlong_v_incid + cc%rlong_v_surf
     endif
     cc => cc%taller
  enddo

  cp%rshort_g_beam = cp%rshort_g_beam * beam_radiation
  cp%rshort_g_diffuse = cp%rshort_g_diffuse *   &
       land%rshort_diffuse(cp%siteptr%iland)
  cp%rshort_g = cp%rshort_g_beam + cp%rshort_g_diffuse

  ! Absorption rate of short wave by the surface water
  do k=1,cp%nlev_sfcwater
     cp%rshort_s_beam(k) = cp%rshort_s_beam(k) * beam_radiation
     cp%rshort_s_diffuse(k) = cp%rshort_s_diffuse(k) *   &
          land%rshort_diffuse(cp%siteptr%iland)
     cp%rshort_s(k) = cp%rshort_s_beam(k) + cp%rshort_s_diffuse(k)
  enddo

  cp%rlong_s_incid = cp%rlong_s_incid * land%rlong(cp%siteptr%iland)
  cp%rlong_g_incid = cp%rlong_g_incid * land%rlong(cp%siteptr%iland)
  cp%rlong_s = cp%rlong_s_surf + cp%rlong_s_incid
  cp%rlong_g = cp%rlong_g_surf + cp%rlong_g_incid

  return
end subroutine scale_ed_radiation
