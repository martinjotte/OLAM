subroutine olam2ed_state_init(cgrid)

  use ed_state_vars,  only: edtype, polygontype
  use ed_consts_coms, only: cp, p00i, rocp, rdry, p00
  use mem_leaf,       only: land, itab_wl
  use sp_therm_lib,   only: thetaeiv
  use misc_coms,      only: isubdomain
  use mem_basic,      only: rho
  use mem_ijtabs,     only: itabg_w

  implicit none

  type(edtype), target :: cgrid
  type(polygontype), pointer :: cpoly
  integer :: ipy, iwl, isi, iw, kw
  real :: rvaux
  real, parameter :: initial_co2=370.0
!  integer, save :: first_call=1


  do ipy = 1, cgrid%npolygons
     iwl = cgrid%iwl(ipy)

     iw = itab_wl(iwl)%iw         ! global index
     kw = itab_wl(iwl)%kw

     ! If run is parallel, get local rank indices
     if (isubdomain == 1) then
        iw = itabg_w(iw)%iw_myrank
     endif

!     if(first_call == 1)then
!        first_call = 0
!        cpoly => cgrid%polygon(ipy)
!        do isi=1,cpoly%nsites
!           cpoly%met(isi)%vels = 2.
!        enddo
!     endif
     
     cgrid%met(ipy)%prss = rho(kw,iw) * rdry * land%cantemp(iwl) *   &
          (1.0 + 0.61 * land%canshv(iwl))
     cgrid%met(ipy)%exner = cp * (p00i * cgrid%met(ipy)%prss)**rocp
     cgrid%met(ipy)%atm_shv = land%canshv(iwl)
     cgrid%met(ipy)%atm_theta = land%cantemp(iwl) * (p00/cgrid%met(ipy)%prss)**rocp
     cgrid%met(ipy)%atm_co2 = initial_co2
     cgrid%met(ipy)%atm_tmp = land%cantemp(iwl)
     cgrid%met(ipy)%rhos = rho(kw,iw)

     !------------------------------------------------------------------------------------!     
     !    We now find the equivalent potential temperature.                               !     
     !------------------------------------------------------------------------------------!     
     rvaux  = cgrid%met(ipy)%atm_shv / (1. - cgrid%met(ipy)%atm_shv)
     cgrid%met(ipy)%atm_theiv = thetaeiv(cgrid%met(ipy)%atm_theta,cgrid%met(ipy)%prss     &
          ,cgrid%met(ipy)%atm_tmp,rvaux,rvaux,1)
     
     !------ Apply met to sites, and adjust met variables for topography. ----------------!     
     call calc_met_lapse(cgrid,ipy)
     
     cpoly => cgrid%polygon(ipy)
     siteloop: do isi = 1,cpoly%nsites
        !----- CO2.  In case we used the namelist, use that value. -----------------------!     
        cpoly%met(isi)%atm_co2 = initial_co2


        !---------------------------------------------------------------------------------!     
        !     We now find some derived properties.  In case several sites exist, the      !     
        ! lapse rate was applied to pressure, temperature, and mixing ratio.  Then we     !     
        ! calculate the Exner function, potential temperature and equivalent potential    !     
        ! temperature, so it will respect the ideal gas law and first law of thermo-      !     
        ! dynamics.                                                                       !     
        !---------------------------------------------------------------------------------!    
        cpoly%met(isi)%exner        = cp * (p00i * cpoly%met(isi)%prss) **rocp
        cpoly%met(isi)%atm_theta    = cp * cpoly%met(isi)%atm_tmp / cpoly%met(isi)%exner

        !----- Find the atmospheric equivalent potential temperature. --------------------!     
        rvaux  = cpoly%met(isi)%atm_shv / (1. - cpoly%met(isi)%atm_shv)
        cpoly%met(isi)%atm_theiv = thetaeiv(cpoly%met(isi)%atm_theta,cpoly%met(isi)%prss  &
             ,cpoly%met(isi)%atm_tmp,rvaux,rvaux,2)
        
     end do siteloop
     
  end do

  return
end subroutine olam2ed_state_init

!==========================================================================================
subroutine olam2ed_state(cgrid)

  use ed_state_vars,  only: edtype, polygontype
  use ed_consts_coms, only: cp, p00i, rocp!, rdry, p00
  use mem_leaf,       only: land, itab_wl
  use misc_coms,      only: isubdomain
  use mem_basic,      only: press, rho, vxe, vye, vze
  use mem_ijtabs,     only: itabg_w
! use sp_therm_lib,   only: thetaeiv

  implicit none

  type(edtype), target :: cgrid
  type(polygontype), pointer :: cpoly
  integer :: ipy, iwl, isi, iw, kw

!  real :: rvaux


  do ipy = 1, cgrid%npolygons
     iwl = cgrid%iwl(ipy)

     iw = itab_wl(iwl)%iw         ! global index
     kw = itab_wl(iwl)%kw

     ! If run is parallel, get local rank indices
     if (isubdomain == 1) then
        iw = itabg_w(iw)%iw_myrank
     endif

     cgrid%met(ipy)%prss = press(kw,iw)
     cgrid%met(ipy)%exner = cp * (p00i * cgrid%met(ipy)%prss)**rocp
     cgrid%met(ipy)%vels = sqrt(vxe(kw,iw)**2 + vye(kw,iw)**2 + vze(kw,iw)**2)
     cgrid%met(ipy)%rhos = rho(kw,iw)

!     cgrid%met(ipy)%atm_shv = land%canshv(iwl)
!     cgrid%met(ipy)%atm_theta = land%cantemp(iwl) * (p00/cgrid%met(ipy)%prss)**rocp
!     cgrid%met(ipy)%atm_co2 = initial_co2
!     cgrid%met(ipy)%atm_tmp = land%cantemp(iwl)

     !------------------------------------------------------------------------------------!     
     !    We now find the equivalent potential temperature.                               !     
     !------------------------------------------------------------------------------------!     
!     rvaux  = cgrid%met(ipy)%atm_shv / (1. - cgrid%met(ipy)%atm_shv)
!     cgrid%met(ipy)%atm_theiv = thetaeiv(cgrid%met(ipy)%atm_theta,cgrid%met(ipy)%prss     &
!          ,cgrid%met(ipy)%atm_tmp,rvaux,rvaux,1)
     
     !------ Apply met to sites, and adjust met variables for topography. ----------------!     
     call calc_met_lapse(cgrid,ipy)
     
     cpoly => cgrid%polygon(ipy)
     siteloop: do isi = 1,cpoly%nsites
  
        !---------------------------------------------------------------------------------!     
        !     We now find some derived properties.  In case several sites exist, the      !     
        ! lapse rate was applied to pressure, temperature, and mixing ratio.  Then we     !     
        ! calculate the Exner function, potential temperature and equivalent potential    !     
        ! temperature, so it will respect the ideal gas law and first law of thermo-      !     
        ! dynamics.                                                                       !     
        !---------------------------------------------------------------------------------!    
        cpoly%met(isi)%exner        = cp * (p00i * cpoly%met(isi)%prss) **rocp
!        cpoly%met(isi)%atm_theta    = cp * cpoly%met(isi)%atm_tmp / cpoly%met(isi)%exner

        !----- Find the atmospheric equivalent potential temperature. --------------------!     
!        rvaux  = cpoly%met(isi)%atm_shv / (1. - cpoly%met(isi)%atm_shv)
!        cpoly%met(isi)%atm_theiv = thetaeiv(cpoly%met(isi)%atm_theta,cpoly%met(isi)%prss  &
!             ,cpoly%met(isi)%atm_tmp,rvaux,rvaux,2)
        
     end do siteloop
     
  end do

  return
end subroutine olam2ed_state

!==========================================================================================

subroutine olam2ed_soils(cgrid)
  use ed_state_vars, only: edtype, polygontype, sitetype
  use mem_leaf, only: land
  use leaf_coms, only: slcpd, nzg, nzs

  implicit none

  type(edtype), target :: cgrid
  type(polygontype), pointer :: cpoly
  type(sitetype), pointer :: csite
  integer :: ipy, iwl, isi, ipa, k

  do ipy = 1, cgrid%npolygons
     cpoly => cgrid%polygon(ipy)
     iwl = cgrid%iwl(ipy)
     do isi = 1, cpoly%nsites
        csite => cpoly%site(isi)
        do ipa = 1, csite%npatches
           do k = 1, nzg
              csite%soil_water(k,ipa) = land%soil_water(k,iwl)
              csite%soil_energy(k,ipa) = land%soil_energy(k,iwl)
              call qwtk(csite%soil_energy(k,ipa),csite%soil_water(k,ipa)*1000.,  &
                   slcpd(cpoly%ntext_soil(k,isi)), csite%soil_tempk(k,ipa),  &
                   csite%soil_fracliq(k,ipa))
           enddo
           csite%nlev_sfcwater(ipa) = land%nlev_sfcwater(iwl)
           do k = 1, nzs
              csite%sfcwater_energy(k,ipa) = land%sfcwater_energy(k,iwl)
              csite%sfcwater_depth(k,ipa) = land%sfcwater_depth(k,iwl)
              csite%sfcwater_mass(k,ipa) = land%sfcwater_mass(k,iwl)
              call qtk(csite%sfcwater_energy(k,ipa), csite%sfcwater_tempk(k,ipa),  &
                   csite%sfcwater_fracliq(k,ipa))
           enddo
        enddo
     enddo
  enddo

  return
end subroutine olam2ed_soils
