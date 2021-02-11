!     path:      $Source$
!     author:    $Author: mike $
!     revision:  $Revision: 11661 $
!     created:   $Date: 2009-05-22 18:22:22 -0400 (Fri, 22 May 2009) $

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

      use parkind,  only : im => kind_im, rb => kind_rb, cldmin, cldmax
      use parrrsw,  only : nbndsw, ngptsw
      use rrsw_wvn, only: ngb

      implicit none

! public interfaces/functions/subroutines
      private
      public :: mcica_subcol_sw, generate_stochastic_clouds_sw

      contains

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine mcica_subcol_sw(nlay, icld, seeds, &
                       cldfrac, ciwp, clwp, tauc, ssac, asmc, &
                       ciwpmcl, clwpmcl, &
                       taucmcl, ssacmcl, asmcmcl, inflg)

! ----- Input -----
! Control

      integer(kind=im), intent(in) :: nlay            ! number of model layers
      integer(kind=im), intent(in) :: icld            ! clear/cloud, cloud overlap flag
      integer(kind=im), intent(in) :: inflg           ! randomize cloud ice/liquid (1)
                                                      ! or optical properties (0)
      integer(kind=im), intent(inout) :: seeds(4)     ! seeds for the kissvec random number generator

! Atmosphere/clouds - cldprop
      real(kind=rb), intent(in) :: cldfrac(nlay)       ! layer cloud fraction
                                                      !    Dimensions: (nlay)
      real(kind=rb), intent(in) :: tauc(nbndsw,nlay)! in-cloud optical depth
                                                      !    Dimensions: (nlay,nbndsw)
      real(kind=rb), intent(in) :: ssac(nbndsw,nlay)! in-cloud single scattering albedo (non-delta scaled)
                                                      !    Dimensions: (nlay,nbndsw)
      real(kind=rb), intent(in) :: asmc(nbndsw,nlay)! in-cloud asymmetry parameter (non-delta scaled)
                                                      !    Dimensions: (nlay,nbndsw)
      real(kind=rb), intent(in) :: ciwp(nlay)          ! in-cloud ice water path
                                                      !    Dimensions: (nlay)
      real(kind=rb), intent(in) :: clwp(nlay)          ! in-cloud liquid water path
                                                      !    Dimensions: (nlay)

! ----- Output -----
! Atmosphere/clouds - cldprmc [mcica]
      real(kind=rb), intent(out) :: ciwpmcl(ngptsw,nlay)    ! in-cloud ice water path [mcica]
                                                      !    Dimensions: (nlay,ngptsw)
      real(kind=rb), intent(out) :: clwpmcl(ngptsw,nlay)    ! in-cloud liquid water path [mcica]
                                                      !    Dimensions: (nlay,ngptsw)
      real(kind=rb), intent(out) :: taucmcl(ngptsw,nlay)    ! in-cloud optical depth [mcica]
                                                      !    Dimensions: (nlay,ngptsw)
      real(kind=rb), intent(out) :: ssacmcl(ngptsw,nlay)    ! in-cloud single scattering albedo [mcica]
                                                      !    Dimensions: (nlay,ngptsw)
      real(kind=rb), intent(out) :: asmcmcl(ngptsw,nlay)    ! in-cloud asymmetry parameter [mcica]
                                                      !    Dimensions: (nlay,ngptsw)

! Return if clear sky; or stop if icld out of range
      if (icld.eq.0) return
      if (icld.lt.0.or.icld.gt.3) then
         stop 'MCICA_SUBCOL: INVALID ICLD'
      endif

! NOTE: For GCM mode, permuteseed must be offset between LW and SW by at least number of subcolumns

!  Generate the stochastic subcolumns of cloud optical properties for the shortwave;
      call generate_stochastic_clouds_sw (nlay, icld, cldfrac, clwp, ciwp, &
                              tauc, ssac, asmc, clwpmcl, ciwpmcl, &
                              taucmcl, ssacmcl, asmcmcl, seeds, inflg)

      end subroutine mcica_subcol_sw


!-------------------------------------------------------------------------------------------------
      subroutine generate_stochastic_clouds_sw(nlay, overlap, cld, clwp, ciwp, &
                             tauc, ssac, asmc, clwp_stoch, ciwp_stoch, &
                             tauc_stoch, ssac_stoch, asmc_stoch, seeds, inflg)
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
  !----------------------------------------------------------------------------------------------

! -- Arguments

      integer(kind=im), intent(in) :: inflg
      integer(kind=im), intent(in) :: nlay            ! number of layers
      integer(kind=im), intent(in) :: overlap         ! clear/cloud, cloud overlap flag
      integer(kind=im), intent(inout) :: seeds(4)     ! seeds for kissvec

! Column state (cloud fraction, cloud water, cloud ice) + variables needed to read physics state
      real(kind=rb), intent(in) :: cld(nlay)           ! cloud fraction
                                                      !    Dimensions: (nlay)
      real(kind=rb), intent(in) :: clwp(nlay)          ! in-cloud liquid water path (g/m2)
                                                      !    Dimensions: (nlay)
      real(kind=rb), intent(in) :: ciwp(nlay)          ! in-cloud ice water path (g/m2)
                                                      !    Dimensions: (nlay)
      real(kind=rb), intent(in) :: tauc(nbndsw,nlay)        ! in-cloud optical depth (non-delta scaled)
                                                      !    Dimensions: (nlay,nbndsw)
      real(kind=rb), intent(in) :: ssac(nbndsw,nlay)        ! in-cloud single scattering albedo (non-delta scaled)
                                                      !    Dimensions: (nlay,nbndsw)
      real(kind=rb), intent(in) :: asmc(nbndsw,nlay)        ! in-cloud asymmetry parameter (non-delta scaled)
                                                      !    Dimensions: (nlay,nbndsw)
      real(kind=rb), intent(out) :: clwp_stoch(ngptsw,nlay) ! subcolumn in-cloud liquid water path
                                                      !    Dimensions: (nlay,ngptsw)
      real(kind=rb), intent(out) :: ciwp_stoch(ngptsw,nlay) ! subcolumn in-cloud ice water path
                                                      !    Dimensions: (nlay,ngptsw)
      real(kind=rb), intent(out) :: tauc_stoch(ngptsw,nlay) ! subcolumn in-cloud optical depth
                                                      !    Dimensions: (nlay,ngptsw)
      real(kind=rb), intent(out) :: ssac_stoch(ngptsw,nlay) ! subcolumn in-cloud single scattering albedo
                                                      !    Dimensions: (nlay,ngptsw)
      real(kind=rb), intent(out) :: asmc_stoch(ngptsw,nlay) ! subcolumn in-cloud asymmetry parameter
                                                      !    Dimensions: (nlay,ngptsw)
! -- Local variables

! Variables related to random number and seed
      real(kind=rb), dimension(ngptsw, nlay) :: CDF          ! random numbers
      real(kind=rb)                          :: rand_num     ! random number (kissvec)

! Indices
      integer(kind=im) :: ilev, ilevp, ilevm, isubcol        ! indices

!------------------------------------------------------------------------------------------

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

! ------ Apply overlap assumption --------

! generate the random numbers

      select case (overlap)

      case(1)

! Random overlap
! i) pick a random value at every level

         do ilev = 1, nlay

            if (cld(ilev) >= cldmax) then
               CDF(1:ngptsw,ilev) = 1.0_rb
            elseif (cld(ilev) >= cldmin) then
               do isubcol = 1,ngptsw
                  call kissvec(seeds, rand_num)
                  CDF(isubcol,ilev) = rand_num
               enddo
            endif

         enddo

      case(2)

! Maximum-Random overlap
! i) pick  a random number for top layer.
! ii) walk down the column:
!    - if the layer above is cloudy, we use the same random number than in the layer above
!    - if the layer above is clear, we use a new random number

         if (cld(1) >= cldmax .and. cld(2) < cldmin) then
            CDF(1:ngptsw,1) = 1.
         elseif (cld(1) >= cldmin) then
            do isubcol = 1, ngptsw
               call kissvec(seeds, rand_num)
               CDF(isubcol,1) = rand_num
            enddo
         endif

         do ilev  = 2, nlay
            ilevp = merge( ilev+1, ilev, ilev<nlay )
            ilevm = ilev - 1

            if (cld(ilev) >= cldmax) then
               if (cld(ilevp) < cldmin .or. ilev == nlay) then
                  CDF(1:ngptsw,ilev) = 1.
               elseif (cld(ilevm) < cldmin) then
                  do isubcol = 1, ngptsw
                     call kissvec(seeds, rand_num)
                     CDF(isubcol,ilev) = rand_num
                  enddo
               else
                  CDF(1:ngptsw,ilev) = CDF(1:ngptsw,ilevm)
               endif
            elseif (cld(ilev) >= cldmin) then
               if (cld(ilevm) < cldmin) then
                  do isubcol = 1, ngptsw
                     call kissvec(seeds, rand_num)
                     CDF(isubcol,ilev) = rand_num
                  enddo
               else
                  do isubcol = 1, ngptsw
                     if (CDF(isubcol, ilevm) >= 1. - cld(ilevm) ) then
                        CDF(isubcol,ilev) = CDF(isubcol,ilevm)
                     else
                        call kissvec(seeds, rand_num)
                        CDF(isubcol,ilev) = (1. - cld(ilevm)) * rand_num
                     endif
                  enddo
               endif
            endif
         enddo

      case(3)
! Maximum overlap
! i) pick same random number at every level

         do isubcol = 1,ngptsw
            call kissvec(seeds, rand_num)
            do ilev = 1,nlay
               CDF(ilev,isubcol) = rand_num
            enddo
         enddo

      end select

! where the subcolumn is cloudy, the subcolumn cloud fraction is 1;
! where the subcolumn is not cloudy, the subcolumn cloud fraction is 0;
! where there is a cloud, define the subcolumn cloud properties,
! otherwise set these to zero

      if (inflg == 0) then

         do ilev = 1,nlay
            if (cld(ilev) >= cldmax) then
               do isubcol = 1, ngptsw
                  tauc_stoch(isubcol,ilev) = tauc(ngb(isubcol),ilev)
                  ssac_stoch(isubcol,ilev) = ssac(ngb(isubcol),ilev)
                  asmc_stoch(isubcol,ilev) = asmc(ngb(isubcol),ilev)
               enddo
            elseif (cld(ilev) >= cldmin) then
               do isubcol = 1, ngptsw
                  if (cdf(isubcol,ilev) >= 1._rb - cld(ilev)) then
                     tauc_stoch(isubcol,ilev) = tauc(ngb(isubcol),ilev)
                     ssac_stoch(isubcol,ilev) = ssac(ngb(isubcol),ilev)
                     asmc_stoch(isubcol,ilev) = asmc(ngb(isubcol),ilev)
                  else
                     tauc_stoch(isubcol,ilev) = 0.
                     ssac_stoch(isubcol,ilev) = 0.
                     asmc_stoch(isubcol,ilev) = 0.
                  endif
               enddo
            endif
         enddo

      else

         do ilev = 1,nlay
            if (cld(ilev) >= cldmax) then
               clwp_stoch(:,ilev) = clwp(ilev)
               ciwp_stoch(:,ilev) = ciwp(ilev)
            elseif (cld(ilev) >= cldmin) then
               do isubcol = 1, ngptsw
                  if (cdf(ilev,isubcol) >= 1._rb - cld(ilev)) then
                     clwp_stoch(isubcol,ilev) = clwp(ilev)
                     ciwp_stoch(isubcol,ilev) = ciwp(ilev)
                  else
                     clwp_stoch(isubcol,ilev) = 0.
                     ciwp_stoch(isubcol,ilev) = 0.
                  endif
               enddo
            endif
         enddo

      endif

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
      integer(kind=im)                 :: kiss

      seed(1) = 69069_im * seed(1) + 1327217885_im
      seed(2) = m (m (m (seed(2), 13_im), - 17_im), 5_im)
      seed(3) = 18000_im * iand (seed(3), 65535_im) + ishft (seed(3), - 16_im)
      seed(4) = 30903_im * iand (seed(4), 65535_im) + ishft (seed(4), - 16_im)
      kiss  = seed(1) + seed(2) + ishft (seed(3), 16_im) + seed(4)
      rand  = kiss*2.328306e-10_rb + 0.5_rb

      contains

        integer function m(k,n)
          implicit none
          integer, intent(in) :: k, n
          m = ieor(k, ishft(k, n))
        end function m

      end subroutine kissvec

      end module mcica_subcol_gen_sw
