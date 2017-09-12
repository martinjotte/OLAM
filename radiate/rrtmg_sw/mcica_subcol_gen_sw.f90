!     path:      $Source$
!     author:    $Author: mike $
!     revision:  $Revision: 11661 $
!     created:   $Date: 2009-05-22 18:22:22 -0400 (Fri, 22 May 2009) $
!

      module mcica_subcol_gen_sw

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2006-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! Purpose: Create McICA stochastic arrays for cloud physical or optical properties.   
! Two options are possible:
! 1) Input cloud physical properties: cloud fraction, ice and liquid water
!    paths, ice fraction, and particle sizes.  Output will be stochastic
!    arrays of these variables.  (inflag = 1)
! 2) Input cloud optical properties directly: cloud optical depth, single
!    scattering albedo and asymmetry parameter.  Output will be stochastic
!    arrays of these variables.  (inflag = 0)

! --------- Modules ----------

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrsw, only : nbndsw, ngptsw
      use rrsw_con, only: grav
      use rrsw_wvn, only: ngb
      use rrsw_vsn

      implicit none

      integer, parameter :: ncol = 1

! public interfaces/functions/subroutines
      public :: mcica_subcol_sw, generate_stochastic_clouds_sw

      contains

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine mcica_subcol_sw(nlay, icld, permuteseed, seeds, irng, &
                       cldfrac, ciwp, clwp, rei, rel, tauc, ssac, asmc, fsfc, &
                       cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                       taucmcl, ssacmcl, asmcmcl, fsfcmcl, inflg)

! ----- Input -----
! Control

      integer(kind=im), intent(in) :: nlay            ! number of model layers
      integer(kind=im), intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer(kind=im), intent(in) :: permuteseed     ! if the cloud generator is called multiple times,
                                                      ! permute the seed between each call;
                                                      ! between calls for LW and SW, recommended
                                                      ! permuteseed differs by 'ngpt'
      integer(kind=im), intent(in) :: inflg           ! randomize cloud ice/liquid (1)
                                                      ! or optical properties (0)
      integer(kind=im), intent(inout) :: seeds(4)     ! seeds for the kissvec random number generator
      integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                      !  0 = kissvec
                                                      !  1 = Mersenne Twister
        
! Atmosphere/clouds - cldprop
      real(kind=rb), intent(in) :: cldfrac(ncol,nlay)       ! layer cloud fraction
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: tauc(nbndsw,ncol,nlay)! in-cloud optical depth
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=rb), intent(in) :: ssac(nbndsw,ncol,nlay)! in-cloud single scattering albedo (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=rb), intent(in) :: asmc(nbndsw,ncol,nlay)! in-cloud asymmetry parameter (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=rb), intent(in) :: fsfc(nbndsw,ncol,nlay)! in-cloud forward scattering fraction (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=rb), intent(in) :: ciwp(ncol,nlay)          ! in-cloud ice water path
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: clwp(ncol,nlay)          ! in-cloud liquid water path
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: rei(ncol,nlay)           ! cloud ice particle size
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: rel(ncol,nlay)           ! cloud liquid particle size
                                                      !    Dimensions: (ncol,nlay)

! ----- Output -----
! Atmosphere/clouds - cldprmc [mcica]
      real(kind=rb), intent(out) :: cldfmcl(ngptsw,ncol,nlay)    ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: ciwpmcl(ngptsw,ncol,nlay)    ! in-cloud ice water path [mcica]
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: clwpmcl(ngptsw,ncol,nlay)    ! in-cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: relqmcl(ncol,nlay)      ! liquid particle size (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(out) :: reicmcl(ncol,nlay)      ! ice partcle size (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(out) :: taucmcl(ngptsw,ncol,nlay)    ! in-cloud optical depth [mcica]
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: ssacmcl(ngptsw,ncol,nlay)    ! in-cloud single scattering albedo [mcica]
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: asmcmcl(ngptsw,ncol,nlay)    ! in-cloud asymmetry parameter [mcica]
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: fsfcmcl(ngptsw,ncol,nlay)    ! in-cloud forward scattering fraction [mcica]
                                                      !    Dimensions: (ngptsw,ncol,nlay)

! ----- Local -----

! Stochastic cloud generator variables [mcica]
      integer(kind=im), parameter :: nsubcsw = ngptsw ! number of sub-columns (g-point intervals)
!     integer(kind=im) :: ilev                        ! loop index

!      real(kind=rb) :: pdel(ncol,nlay)               ! layer pressure thickness (Pa) 
!      real(kind=rb) :: qi(ncol,nlay)                 ! ice water (specific humidity)
!      real(kind=rb) :: ql(ncol,nlay)                 ! liq water (specific humidity)


! Return if clear sky; or stop if icld out of range
      if (icld.eq.0) return
      if (icld.lt.0.or.icld.gt.3) then 
         stop 'MCICA_SUBCOL: INVALID ICLD'
      endif 

! NOTE: For GCM mode, permuteseed must be offset between LW and SW by at least number of subcolumns


! Pass particle sizes to new arrays, no subcolumns for these properties yet
! Convert pressures from mb to Pa

      if (inflg > 0) then
         reicmcl(1:ncol,1:nlay) = rei(1:ncol,1:nlay)
         relqmcl(1:ncol,1:nlay) = rel(1:ncol,1:nlay)
      endif

! Convert input ice and liquid cloud water paths to specific humidity ice and liquid components 

!      cwp =  (q * pdel * 1000.) / gravit)
!           = (kg/kg * kg m-1 s-2 *1000.) / m s-2
!           = (g m-2)
!
!      q  = (cwp * gravit) / (pdel *1000.)
!         = (g m-2 * m s-2) / (kg m-1 s-2 * 1000.)
!         =  kg/kg

!      do ilev = 1, nlay
!         qi(ilev) = (ciwp(ilev) * grav) / (pdel(ilev) * 1000._rb)
!         ql(ilev) = (clwp(ilev) * grav) / (pdel(ilev) * 1000._rb)
!      enddo

!  Generate the stochastic subcolumns of cloud optical properties for the shortwave;
      call generate_stochastic_clouds_sw (nlay, nsubcsw, icld, irng, cldfrac, clwp, ciwp, &
                              tauc, ssac, asmc, fsfc, cldfmcl, clwpmcl, ciwpmcl, &
                              taucmcl, ssacmcl, asmcmcl, fsfcmcl, permuteseed, seeds, inflg)

      end subroutine mcica_subcol_sw


!-------------------------------------------------------------------------------------------------
      subroutine generate_stochastic_clouds_sw(nlay, nsubcol, icld, irng, cld, clwp, ciwp, &
                             tauc, ssac, asmc, fsfc, cld_stoch, clwp_stoch, ciwp_stoch, &
                             tauc_stoch, ssac_stoch, asmc_stoch, fsfc_stoch, changeSeed, seeds, inflg)
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
  ! ---------------------
  ! Contact: Cecile Hannay (hannay@ucar.edu)
  ! 
  ! Original code: Based on Raisanen et al., QJRMS, 2004.
  !
  ! Modifications: Generalized for use with RRTMG and added Mersenne Twister as the default
  !   random number generator, which can be changed to the optional kissvec random number generator
  !   with flag 'irng'. Some extra functionality has been commented or removed.  
  !   Michael J. Iacono, AER, Inc., February 2007
  !
  ! Given a profile of cloud fraction, cloud water and cloud ice, we produce a set of subcolumns.
  ! Each layer within each subcolumn is homogeneous, with cloud fraction equal to zero or one 
  ! and uniform cloud liquid and cloud ice concentration.
  ! The ensemble as a whole reproduces the probability function of cloud liquid and ice within each layer 
  ! and obeys an overlap assumption in the vertical.   
  ! 
  ! Overlap assumption:
  !  The cloud are consistent with 4 overlap assumptions: random, maximum, maximum-random and exponential. 
  !  The default option is maximum-random (option 3)
  !  The options are: 1=random overlap, 2=max/random, 3=maximum overlap, 4=exponential overlap
  !  This is set with the variable "overlap" 
  !mji - Exponential overlap option (overlap=4) has been deactivated in this version
  !  The exponential overlap uses also a length scale, Zo. (real,    parameter  :: Zo = 2500. ) 
  ! 
  ! Seed:
  !  If the stochastic cloud generator is called several times during the same timestep, 
  !  one should change the seed between the call to insure that the subcolumns are ifferent.
  !  This is done by changing the argument 'changeSeed'
  !  For example, if one wants to create a set of columns for the shortwave and another set for the longwave ,
  !  use 'changeSeed = 1' for the first call and'changeSeed = 2' for the second call 
  !
  ! PDF assumption:
  !  We can use arbitrary complicated PDFS. 
  !  In the present version, we produce homogeneuous clouds (the simplest case).  
  !  Future developments include using the PDF scheme of Ben Johnson. 
  !
  ! History file:
  !  Option to add diagnostics variables in the history file. (using FINCL in the namelist)
  !  nsubcol = number of subcolumns
  !  overlap = overlap type (1-3)
  !  Zo = length scale 
  !  CLOUD_S = mean of the subcolumn cloud fraction ('_S" means Stochastic)
  !  CLDLIQ_S = mean of the subcolumn cloud water
  !  CLDICE_S = mean of the subcolumn cloud ice 
  !
  ! Note:
  !   Here: we force that the cloud condensate to be consistent with the cloud fraction 
  !   i.e we only have cloud condensate when the cell is cloudy. 
  !   In CAM: The cloud condensate and the cloud fraction are obtained from 2 different equations 
  !   and the 2 quantities can be inconsistent (i.e. CAM can produce cloud fraction 
  !   without cloud condensate or the opposite).
  !---------------------------------------------------------------------------------------------------------------

      use mcica_random_numbers
! The Mersenne Twister random number engine
      use MersenneTwister, only: randomNumberSequence, &   
                                 new_RandomNumberSequence, getRandomReal

      type(randomNumberSequence) :: randomNumbers

! -- Arguments

      integer(kind=im), intent(in) :: inflg
      integer(kind=im), intent(in) :: nlay            ! number of layers
      integer(kind=im), intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer(kind=im), intent(inout) :: irng         ! flag for random number generator
                                                      !  0 = kissvec
                                                      !  1 = Mersenne Twister
      integer(kind=im), intent(in) :: nsubcol         ! number of sub-columns (g-point intervals)
      integer(kind=im), intent(inout) :: seeds(4)     ! seeds for kissvec
      integer(kind=im), optional, intent(in) :: changeSeed     ! allows permuting seed

! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state 
      real(kind=rb), intent(in) :: cld(ncol,nlay)           ! cloud fraction 
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: clwp(ncol,nlay)          ! in-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: ciwp(ncol,nlay)          ! in-cloud ice water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: tauc(nbndsw,ncol,nlay)        ! in-cloud optical depth (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=rb), intent(in) :: ssac(nbndsw,ncol,nlay)        ! in-cloud single scattering albedo (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=rb), intent(in) :: asmc(nbndsw,ncol,nlay)        ! in-cloud asymmetry parameter (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=rb), intent(in) :: fsfc(nbndsw,ncol,nlay)        ! in-cloud forward scattering fraction (non-delta scaled)
                                                      !    Dimensions: (nbndsw,ncol,nlay)

      real(kind=rb), intent(out) :: cld_stoch(ngptsw,ncol,nlay)  ! subcolumn cloud fraction 
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: clwp_stoch(ngptsw,ncol,nlay) ! subcolumn in-cloud liquid water path
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: ciwp_stoch(ngptsw,ncol,nlay) ! subcolumn in-cloud ice water path
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: tauc_stoch(ngptsw,ncol,nlay) ! subcolumn in-cloud optical depth
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: ssac_stoch(ngptsw,ncol,nlay) ! subcolumn in-cloud single scattering albedo
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: asmc_stoch(ngptsw,ncol,nlay) ! subcolumn in-cloud asymmetry parameter
                                                      !    Dimensions: (ngptsw,ncol,nlay)
      real(kind=rb), intent(out) :: fsfc_stoch(ngptsw,ncol,nlay) ! subcolumn in-cloud forward scattering fraction
                                                      !    Dimensions: (ngptsw,ncol,nlay)

! -- Local variables
      real(kind=rb) :: cldf(ncol,nlay)                ! cloud fraction 
                                                      !    Dimensions: (ncol,nlay)

! Mean over the subcolumns (cloud fraction, cloud water , cloud ice) - inactive
!      real(kind=rb) :: mean_cld_stoch(ncol,nlay)     ! cloud fraction 
!      real(kind=rb) :: mean_clwp_stoch(ncol,nlay)    ! cloud water
!      real(kind=rb) :: mean_ciwp_stoch(ncol,nlay)    ! cloud ice
!      real(kind=rb) :: mean_tauc_stoch(ncol,nlay)    ! cloud optical depth
!      real(kind=rb) :: mean_ssac_stoch(ncol,nlay)    ! cloud single scattering albedo
!      real(kind=rb) :: mean_asmc_stoch(ncol,nlay)    ! cloud asymmetry parameter
!      real(kind=rb) :: mean_fsfc_stoch(ncol,nlay)    ! cloud forward scattering fraction

! Set overlap
      integer(kind=im) :: overlap                     ! 1 = random overlap, 2 = maximum/random,
                                                      ! 3 = maximum overlap, 
!      real(kind=rb), parameter  :: Zo = 2500._rb        ! length scale (m) 
!      real(kind=rb) :: zm(ncon,nlay)                    ! Height of midpoints (above surface)
!      real(kind=rb), dimension(nlay) :: alpha=0.0_rb    ! overlap parameter  

! Constants (min value for cloud fraction and cloud water and ice)
      real(kind=rb), parameter :: cldmin = 1.0e-7_rb ! min cloud fraction
!      real(kind=rb), parameter :: qmin   = 1.0e-10_rb   ! min cloud water and cloud ice (not used)

! Variables related to random number and seed 
      real(kind=rb), dimension(nsubcol, ncol, nlay) :: CDF !, CDF2       ! random numbers
      real(kind=rb)                                 :: rand_num          ! random number (kissvec)

! Flag to identify cloud fraction in subcolumns
!     logical,  dimension(nsubcol, ncol, nlay) :: isCloudy   ! flag that says whether a gridbox is cloudy

! Indices
      integer(kind=im) :: ilev, ilevp, isubcol, i, n, ngbm    ! indices

!------------------------------------------------------------------------------------------ 

! Check that irng is in bounds; if not, set to default
      if (irng .ne. 0) irng = 1

! Pass input cloud overlap setting to local variable
      overlap = icld

! Ensure that cloud fractions are in bounds 
      do ilev = 1, nlay
         do i = 1, ncol
            cldf(i,ilev) = cld(i,ilev)
            if (cldf(i,ilev) < cldmin) then
               cldf(i,ilev) = 0._rb
            endif
         enddo
      enddo

! ----- Create seed  --------
   
! Advance randum number generator by changeseed values
!      if (irng.eq.0) then   
! For kissvec, create a seed that depends on the state of the columns. Maybe not the best way, but it works.  
! Must use pmid from bottom four layers. 
!        do i=1,ncol
!            if (pmid(i,1).lt.pmid(i,2)) then
!               stop 'MCICA_SUBCOL: KISSVEC SEED GENERATOR REQUIRES PMID FROM BOTTOM FOUR LAYERS.'
!            endif
!            seed1(i) = (pmid(i,1) - int(pmid(i,1)))  * 1000000000_im
!            seed2(i) = (pmid(i,2) - int(pmid(i,2)))  * 1000000000_im
!            seed3(i) = (pmid(i,3) - int(pmid(i,3)))  * 1000000000_im
!            seed4(i) = (pmid(i,4) - int(pmid(i,4)))  * 1000000000_im
!         do i=1,changeSeed
!            call kissvec(seed1, seed2, seed3, seed4, rand_num)
!         enddo
      if (irng.eq.1) then
         randomNumbers = new_RandomNumberSequence(seed = changeSeed)
      endif 

! ------ Apply overlap assumption --------

! generate the random numbers  

      select case (overlap)

      case(1) 
! Random overlap
! i) pick a random value at every level

         do ilev = 1, nlay
            do i = 1, ncol

               if (cldf(i,ilev) < cldmin) then
                  CDF(1:nsubcol,i,ilev) = 0.0_rb
               elseif(cldf(i,ilev) > 0.999) then
                  CDF(1:nsubcol,i,ilev) = 1.0_rb
               elseif (irng == 0) then
                  do isubcol = 1,nsubcol
                     call kissvec(seeds, rand_num)
                     CDF(isubcol,i,ilev) = rand_num
                  enddo
               else
                  do isubcol = 1, nsubcol
                     CDF(isubcol,i,ilev) = getRandomReal(randomNumbers)
                  enddo
               endif

            enddo
         enddo

      case(2)
! Maximum-Random overlap
! i) pick  a random number for top layer.
! ii) walk down the column: 
!    - if the layer above is cloudy, we use the same random number than in the layer above
!    - if the layer above is clear, we use a new random number 

         ilev  = 1
         ilevp = 2
         do i = 1, ncol
            if (cldf(i,ilev) > 0.999 .and. cldf(i,ilevp) < cldmin) then
               CDF(1:nsubcol,i,ilev) = 1._rb
            elseif (cldf(i,ilev) < cldmin) then
               CDF(1:nsubcol,i,ilev) = 0.0_rb
            elseif (irng == 0) then
               do isubcol = 1,nsubcol
                  call kissvec(seeds, rand_num)
                  CDF(isubcol,i,ilev) = rand_num
               enddo
            else
               do isubcol = 1, nsubcol
                  CDF(isubcol,i,ilev) = getRandomReal(randomNumbers)
               enddo
            endif
         enddo

         do ilev = 2,nlay
            ilevp = merge( ilev+1, ilev, ilev<nlay )
            do i = 1, ncol
               if (cldf(i,ilev) > 0.999 .and. cldf(i,ilevp) < cldmin) then
                  CDF(1:nsubcol,i,ilev) = 1._rb
               elseif (cldf(i,ilev) < cldmin) then
                  CDF(1:nsubcol,i,ilev) = 0.0_rb
               elseif (cldf(i,ilev) == cldf(i,ilev-1)) then
                  CDF(1:nsubcol,i,ilev) = CDF(1:nsubcol,i,ilev-1) 
               else 
                  do isubcol = 1, nsubcol
                     if (CDF(isubcol, i, ilev-1) > 1._rb - cldf(i,ilev-1) ) then
                        CDF(isubcol,i,ilev) = CDF(isubcol,i,ilev-1) 
                     elseif (irng == 0) then
                        call kissvec(seeds, rand_num)
                        CDF(isubcol,i,ilev) = (1._rb - cldf(i,ilev-1)) * rand_num
                     else
                        CDF(isubcol,i,ilev) = (1._rb - cldf(i,ilev-1)) &
                                            * getRandomReal(randomNumbers)
                     endif
                  enddo
               endif
            enddo
         enddo

      case(3) 
! Maximum overlap
! i) pick same random numebr at every level  

         if (irng.eq.0) then 
            do isubcol = 1,nsubcol
               do i = 1, ncol
                  call kissvec(seeds, rand_num)
                  do ilev = 1,nlay
                     CDF(isubcol,i,ilev) = rand_num
                  enddo
               enddo
            enddo
         elseif (irng.eq.1) then
            do isubcol = 1, nsubcol
               do i = 1, ncol
                  rand_num = getRandomReal(randomNumbers)
                  do ilev = 1, nlay
                     CDF(isubcol,i,ilev) = rand_num
                  enddo
               enddo
             enddo
         endif

!    case(4)  - inactive
!       ! Exponential overlap: weighting between maximum and random overlap increases with the distance. 
!       ! The random numbers for exponential overlap verify:
!       ! j=1   RAN(j)=RND1
!       ! j>1   if RND1 < alpha(j,j-1) => RAN(j) = RAN(j-1)
!       !                                 RAN(j) = RND2
!       ! alpha is obtained from the equation
!       ! alpha = exp(- (Zi-Zj-1)/Zo) where Zo is a characteristic length scale    


!       ! compute alpha
!       zm    = state%zm     
!       alpha(:, 1) = 0._rb
!       do ilev = 2,nlay
!          alpha(:, ilev) = exp( -( zm (:, ilev-1) -  zm (:, ilev)) / Zo)
!       end do
       
!       ! generate 2 streams of random numbers
!       do isubcol = 1,nsubcol
!          do ilev = 1,nlay
!             call kissvec(seeds, rand_num)
!             CDF(isubcol, :, ilevvv) = rand_num
!             call kissvec(seeds, rand_num)
!             CDF2(isubcol, :, ilev) = rand_num
!          end do
!       end do

!       ! generate random numbers
!       do ilev = 2,nlay
!          where (CDF2(:, :, ilev) < spread(alpha (:,ilev), dim=1, nCopies=nsubcol) )
!             CDF(:,:,ilev) = CDF(:,:,ilev-1) 
!          end where
!       end do

      end select

 
! -- generate subcolumns for homogeneous clouds -----
!      do ilev = 1, nlay
!         isCloudy(:,:,ilev) = (CDF(:,:,ilev) >= 1._rb - spread(cldf(:,ilev), dim=1, nCopies=nsubcol) )
!      enddo

! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
! where there is a cloud, define the subcolumn cloud properties,
! otherwise set these to zero

      ngbm = ngb(1) - 1

      if (inflg == 0) then

         do ilev = 1,nlay
            do i = 1, ncol
               do isubcol = 1, nsubcol
                  if (cdf(isubcol,i,ilev) > 1._rb - cldf(i,ilev)) then
                     n = ngb(isubcol) - ngbm
                     cld_stoch (isubcol,i,ilev) = 1._rb
                     tauc_stoch(isubcol,i,ilev) = tauc(n,i,ilev)
                     ssac_stoch(isubcol,i,ilev) = ssac(n,i,ilev)
                     asmc_stoch(isubcol,i,ilev) = asmc(n,i,ilev)
                     fsfc_stoch(isubcol,i,ilev) = fsfc(n,i,ilev)
                  else
                     cld_stoch(isubcol,i,ilev) = 0._rb
                     tauc_stoch(isubcol,i,ilev) = 0._rb
                     ssac_stoch(isubcol,i,ilev) = 1._rb
                     asmc_stoch(isubcol,i,ilev) = 0._rb
                     fsfc_stoch(isubcol,i,ilev) = 0._rb
                  endif
               enddo
            enddo
         enddo

      else

         do ilev = 1,nlay
            do i = 1, ncol
               do isubcol = 1, nsubcol
                  if (cdf(isubcol,i,ilev) > 1._rb - cldf(i,ilev)) then
                     cld_stoch (isubcol,i,ilev) = 1._rb
                     clwp_stoch(isubcol,i,ilev) = clwp(i,ilev)
                     ciwp_stoch(isubcol,i,ilev) = ciwp(i,ilev)
                  else
                     cld_stoch (isubcol,i,ilev) = 0._rb
                     clwp_stoch(isubcol,i,ilev) = 0._rb
                     ciwp_stoch(isubcol,i,ilev) = 0._rb
                  endif
               enddo
            enddo
         enddo
         
      endif

! -- compute the means of the subcolumns ---
!      mean_cld_stoch(:,:) = 0._rb
!      mean_clwp_stoch(:,:) = 0._rb
!      mean_ciwp_stoch(:,:) = 0._rb
!      mean_tauc_stoch(:,:) = 0._rb
!      mean_ssac_stoch(:,:) = 0._rb
!      mean_asmc_stoch(:,:) = 0._rb
!      mean_fsfc_stoch(:,:) = 0._rb
!      do i = 1, nsubcol
!         mean_cld_stoch(:,:) =  cld_stoch(i,:,:) + mean_cld_stoch(:,:) 
!         mean_clwp_stoch(:,:) =  clwp_stoch( i,:,:) + mean_clwp_stoch(:,:) 
!         mean_ciwp_stoch(:,:) =  ciwp_stoch( i,:,:) + mean_ciwp_stoch(:,:) 
!         mean_tauc_stoch(:,:) =  tauc_stoch( i,:,:) + mean_tauc_stoch(:,:) 
!         mean_ssac_stoch(:,:) =  ssac_stoch( i,:,:) + mean_ssac_stoch(:,:) 
!         mean_asmc_stoch(:,:) =  asmc_stoch( i,:,:) + mean_asmc_stoch(:,:) 
!         mean_fsfc_stoch(:,:) =  fsfc_stoch( i,:,:) + mean_fsfc_stoch(:,:) 
!      end do
!      mean_cld_stoch(:,:) = mean_cld_stoch(:,:) / nsubcol
!      mean_clwp_stoch(:,:) = mean_clwp_stoch(:,:) / nsubcol
!      mean_ciwp_stoch(:,:) = mean_ciwp_stoch(:,:) / nsubcol
!      mean_tauc_stoch(:,:) = mean_tauc_stoch(:,:) / nsubcol
!      mean_ssac_stoch(:,:) = mean_ssac_stoch(:,:) / nsubcol
!      mean_asmc_stoch(:,:) = mean_asmc_stoch(:,:) / nsubcol
!      mean_fsfc_stoch(:,:) = mean_fsfc_stoch(:,:) / nsubcol

      end subroutine generate_stochastic_clouds_sw


!-------------------------------------------------------------------------------------------------- 
      subroutine kissvec(seed,rand)
!-------------------------------------------------------------------------------------------------- 

! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04 for NCAR CAM
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004

! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123; 
!
      real   (kind=rb), intent(out)    :: rand
      integer(kind=im), intent(inout)  :: seed(4)
      integer(kind=im) :: kiss
      integer(kind=im) :: m, k, n

! inline function 
      m(k, n) = ieor (k, ishft (k, n) )

      seed(1) = 69069_im * seed(1) + 1327217885_im
      seed(2) = m (m (m (seed(2), 13_im), - 17_im), 5_im)
      seed(3) = 18000_im * iand (seed(3), 65535_im) + ishft (seed(3), - 16_im)
      seed(4) = 30903_im * iand (seed(4), 65535_im) + ishft (seed(4), - 16_im)
      kiss  = seed(1) + seed(2) + ishft (seed(3), 16_im) + seed(4)
      rand  = kiss*2.328306e-10_rb + 0.5_rb
      
      end subroutine kissvec

      end module mcica_subcol_gen_sw


