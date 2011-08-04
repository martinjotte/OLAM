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
subroutine cuparth_shal(iw, dtime)

  use mem_tend,    only: thilt, sh_wt
  use mem_cuparm,  only: thsrcsh, rtsrcsh
  use mem_grid, only: mwa, mza, lpw, zm, zt, arw0
  use consts_coms, only: grav, rdry, p00, rocp, tkmin
  use mem_basic, only: theta, press, rho, sh_v
  use mem_micro,  only: sh_c
  use grell_coms, only: iact_gr
  use mem_turb, only: tkep
  use misc_coms, only: io6, nqparm_sh

  implicit none

  integer, intent(in) :: iw

  real, intent(in) :: dtime

  integer :: k2
  integer :: mrl
  integer :: j
  integer :: kr
  integer :: k

  real, dimension(mza) :: tkeg
  real, dimension(mza) :: rcpg
  real, dimension(mza) :: t
  real, dimension(mza) :: tn
  real, dimension(mza) :: q
  real, dimension(mza) :: qo
  real, dimension(mza) :: p
  real, dimension(mza) :: po
  real, dimension(mza) :: outt
  real, dimension(mza) :: outq
  real :: z1
  real :: psur
  integer :: tke_calc
  real :: tscl_KF

! Initialize

  thsrcsh(:,iw) = 0.0
  rtsrcsh(:,iw) = 0.0

! A. Betts: suggestion for the KF timescale < DELTAX  / 25 m/s

  tscl_KF = sqrt(arw0(iw)) / 25.0

!------- transfere valores do rams para o eschema

  z1 = zm(lpw(iw)-1)

! surface pressure in mb

  psur = exp(  &
       grav * (zt(lpw(iw)) - z1) / ( rdry * theta(lpw(iw),iw) *   &
       (sngl(press(lpw(iw),iw)) / p00)**rocp *   &
       (1.0+0.608*sh_v(lpw(iw),iw))) + log(sngl(press(lpw(iw),iw))*0.01))

  k2 = mza - lpw(iw)

  do kr = lpw(iw), mza-1
     k = kr - lpw(iw) + 1
           
     po(k) = sngl(press(kr,iw)) * 1.0e-2 ! pressure in mbar

     if (allocated(tkep)) then
        tkeg(k) = tkep(kr,iw)
        tke_calc = 1
     else
        tkeg(k) = 0.0
        tke_calc = 0
     endif
     
     rcpg(k) = 0.

     if (allocated(sh_c)) then
        rcpg(k) = sh_c(kr,iw)
     endif
      
     t(k) = theta(kr,iw) * (sngl(press(kr,iw)) / p00)**rocp
     q(k) = sh_v(kr,iw)

     tn(k) = t(k) + thilt(kr,iw) * (sngl(press(kr,iw)) / p00)**rocp *   &
          dtime / sngl(rho(kr,iw))
     qo(k) = q(k) + sh_wt(kr,iw) * dtime / sngl(rho(kr,iw))

     if (tn(k) < 200.0) tn(k) = t(k)
     qo(k) = max(qo(k),1.0e-8)

     p(k) = po(k)

     outt(k) = 0.0
     outq(k) = 0.0

  enddo

  call CUP_enss_shal(k2, iact_gr(iw), z1, psur, tkmin, dtime, t,   &
       tn, p, po, q, qo, rcpg, tkeg, outt, outq, tke_calc, tscl_KF)

  do k = 1, k2
     kr = k + lpw(iw) - 1
     thsrcsh(kr,iw) = theta(kr,iw) / t(k) * outt(k)
     rtsrcsh(kr,iw) = outq(k)
  enddo

  return
end subroutine cuparth_shal

!=====================================================================

subroutine CUP_enss_shal(k2, iact_gr, z1, psur, tkmin, dtime, t,   &
     tn, p, po, q, qo, rcpg, tkeg, outt, outq, tke_calc, tscl_KF)

  use grell_shallow_coms, only: entr_rate, cap_maxs, maxens, ensdim,  &
       maxens2, maxens3, cap_inc, zkbmax, icoic
  use consts_coms, only: cpi, alvl, t00
  use misc_coms,   only: io6

  implicit none

  integer, intent(in) :: k2
  integer, intent(in) :: iact_gr
  real, intent(in) :: z1
  real, intent(in) :: psur
  real, intent(in) :: tkmin
  real, intent(in) :: dtime
  real, intent(in) :: tscl_KF
  real, dimension(k2), intent(in) :: t
  real, dimension(k2), intent(in) :: tn
  real, dimension(k2), intent(in) :: p
  real, dimension(k2), intent(in) :: po
  real, dimension(k2), intent(in) :: rcpg
  real, dimension(k2), intent(in) :: tkeg
  real, dimension(k2), intent(inout) :: q
  real, dimension(k2), intent(inout) :: qo
  real, dimension(k2), intent(inout) :: outt
  real, dimension(k2), intent(inout) :: outq

  real :: mentr_rate ! entrainment of mass
  real, dimension(k2) :: cd
  real :: cap_max
  integer :: ierr
  integer :: nens

  integer :: kzi
  integer :: iedt
  integer :: k22
  integer :: kbcon
  integer :: k
  integer :: ktop
  integer :: kbmax
  real :: hkb
  real :: xhkb
  real :: hkbo
  real :: mbdt
  real :: xmb
  real :: aa0
  real :: xaa0
  real :: aa1
  real :: pwav
  real :: xpwav
  real :: pwavo
  real, dimension(k2) :: t_cup
  real, dimension(k2) :: tn_cup
  real, dimension(k2) :: p_cup
  real, dimension(k2) :: po_cup
  real, dimension(k2) :: q_cup
  real, dimension(k2) :: qo_cup
  real, dimension(k2) :: z
  real, dimension(k2) :: zo
  real, dimension(k2) :: z_cup
  real, dimension(k2) :: zo_cup
  real, dimension(k2) :: he
  real, dimension(k2) :: xhe
  real, dimension(k2) :: xq
  real, dimension(k2) :: xhes
  real, dimension(k2) :: xqes
  real, dimension(k2) :: xz
  real, dimension(k2) :: xt
  real, dimension(k2) :: heo
  real, dimension(k2) :: he_cup
  real, dimension(k2) :: heo_cup
  real, dimension(k2) :: qes
  real, dimension(k2) :: qeso
  real, dimension(k2) :: qes_cup
  real, dimension(k2) :: qeso_cup
  real, dimension(k2) :: hes
  real, dimension(k2) :: heso
  real, dimension(k2) :: hc
  real, dimension(k2) :: xhc
  real, dimension(k2) :: zu
  real, dimension(k2) :: xzu
  real, dimension(k2) :: zuo
  real, dimension(k2) :: dby
  real, dimension(k2) :: xdby
  real, dimension(k2) :: hco
  real, dimension(k2) :: dbyo
  real, dimension(k2) :: hes_cup
  real, dimension(k2) :: heso_cup
  real, dimension(k2) :: gamma_cup
  real, dimension(k2) :: gammao_cup
  real, dimension(k2) :: qc
  real, dimension(k2) :: xqc
  real, dimension(k2) :: qrc
  real, dimension(k2) :: xqrc
  real, dimension(k2) :: qco
  real, dimension(k2) :: qrco
  real, dimension(k2) :: pw
  real, dimension(k2) :: xpw
  real, dimension(k2) :: dellaq
  real, dimension(k2) :: dellat
  real, dimension(k2) :: dellah
  real, dimension(k2) :: xq_cup
  real, dimension(k2) :: xz_cup
  real, dimension(k2) :: xt_cup
  real, dimension(k2) :: xqes_cup
  real, dimension(k2) :: xhes_cup
  real, dimension(k2) :: xhe_cup
  real, dimension(k2) :: pwo
  real, dimension(k2, maxens2) :: dellat_ens
  real, dimension(k2, maxens2) :: dellaq_ens
  real, dimension(maxens) :: mbdt_ens
  real, dimension(maxens) :: xaa0_ens
  real, dimension(ensdim) :: xf_ens
  integer, intent(in) :: tke_calc

  !--- initial detrainment rates
  ! strong lateral mixing, entrtainment=dedtrainment,
  ! constant mass flux with height
  mentr_rate = entr_rate
  cd(1:k2) = entr_rate
  cap_max = cap_maxs
     
  !don't permit shallow if there is deep convection
  if(iact_gr > 10000 ) then
     ierr = 20
  else
     ierr = 0
  endif

  do nens=1,maxens
     mbdt_ens(nens) = dtime * ( (nens - 3) * 1.0e-3 + 5.0e-3 )
  enddo

!--- environmental conditions, FIRST HEIGHTS
  if(ierr /= 20)then
     xf_ens(1:ensdim)= 0.0
  endif

!--- calculate moist static energy, heights, qes
  if(ierr == 0)then
     call cup_env(k2, t, p, t00, 0, z1, psur, q, he, hes, qes, z) 
     call cup_env(k2, tn, po, t00, 0, z1, psur, qo, heo, heso, qeso, zo) 

!--- environmental values on cloud levels
     call cup_env_clev(k2, psur, z1, q, q_cup, z, z_cup, p, p_cup, t,   &
          t_cup, qes, qes_cup, hes, hes_cup, he, he_cup, gamma_cup)
     call cup_env_clev(k2, psur, z1, qo, qo_cup, zo, zo_cup,  &
       po, po_cup, tn, tn_cup, qeso, qeso_cup, heso, heso_cup, heo, heo_cup,  &
       gammao_cup)

     find_kbmax: do k=1,k2
        if(zo_cup(k) > (zkbmax + z1))then
           kbmax = k
           exit find_kbmax
        endif
     enddo find_kbmax


!--- DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!
!Gr-dec2002
!   here, ,should just take top of PBL for shallow clouds, also maybe for
!   deep clouds in tropics: k22=level(pbltop)
!
!srf-fev2003

!       New way to define K22
!----   Determine PBL top using TKE (TKEG) and liquid water mixing ratio (RCPG)
     if(tke_calc == 1)then
        call get_zi(k2, zo, rcpg, tkeg, z1, tkmin, kzi)
        k22 = max(2, kzi - 1)
     else
        call maximi(k2, heo_cup, 2, kbmax, ierr, k22)
     endif

!srf-14-fev-2003
!A segunda forma produz uma cobertura de shallow mais representativa
!Em caso de alta resolucao vertical tente a versao original (k22=kzi)

     if(k22 >= kbmax)ierr = 2

!--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON

     call cup_kbcon(k22, kbcon, heo_cup, heso_cup, ierr,   &
          kbmax, po_cup, cap_max, k2)

!--- Calculate incloud moist static energy
     call cup_up_he(k2, k22, hkb, z_cup, cd, mentr_rate,   &
          he_cup, hc, kbcon, dby, he, hes_cup)
     call cup_up_he(k2, k22, hkbo, zo_cup, cd, mentr_rate,   &
          heo_cup, hco, kbcon, dbyo, heo, heso_cup)
    
  else ! ierr /= 0.

     kbcon = 1

  endif

!--- Determine cloud top - KTOP
  call cup_ktop(1, dbyo, kbcon, ktop, k2, ierr)

!--- Normalized updraft mass flux profile
  call cup_up_nms(zu, z_cup, mentr_rate, cd, kbcon, ktop, k2, ierr, k22)
  call cup_up_nms(zuo, zo_cup, mentr_rate, cd, kbcon, ktop, k2, ierr, k22)

!--- Calculate moisture properties of updraft
  call cup_up_moisture(k2, ierr, z_cup, qc, qrc, pw, pwav, &
     kbcon, ktop, cd, dby, mentr_rate, q, gamma_cup, zu, qes_cup,   &
     k22, q_cup)
  call cup_up_moisture(k2, ierr, zo_cup, qco, qrco, pwo, pwavo, &
     kbcon, ktop, cd, dbyo, mentr_rate, qo, gammao_cup, zuo, qeso_cup,   &
     k22, qo_cup)

!--- Calculate workfunctions for updrafts
  call cup_up_aa0(k2, aa0, z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr) 
  call cup_up_aa0(k2, aa1, zo, zuo, dbyo, gammao_cup, tn_cup, kbcon,   &
       ktop, ierr) 

  if(ierr == 0 .and. aa1 == 0.0)ierr = 17

!srf - Big loop starts here!

  iedt = 1
  dellat_ens(1:k2,iedt)=0.
  dellaq_ens(1:k2,iedt)=0.

!--- changes per unit mass from shallow convection

  call cup_dellas_shallow(k2, ierr, zo_cup, po_cup, heo, dellah,   &
       zuo, cd, hco, ktop, k22, kbcon, mentr_rate, heo_cup)

  call cup_dellas_shallow(k2, ierr, zo_cup, po_cup, qo, dellaq,   &
       zuo, cd, qco, ktop, k22, kbcon, mentr_rate, qo_cup)

!--- Using dellas, calculate changed environmental profiles
  
  do nens = 1, maxens

     mbdt = mbdt_ens(nens)

     xaa0_ens(nens)=0.0

     do k = 1, k2 - 1
        dellat(k) = 0.0
        if(ierr == 0)then
           xhe(k) = dellah(k) * mbdt + heo(k)
           xq(k) = dellaq(k) * mbdt + qo(k)
           dellat(k) = cpi * (dellah(k) - alvl * dellaq(k))
           xt(k) = dellat(k) * mbdt + tn(k)
           if(xq(k) <= 0.) xq(k) = 1.e-08
        endif
     enddo

     dellat(k2) = 0.0
     if(ierr == 0)then
        xhe(k2) = heo(k2)
        xq(k2) = qo(k2)
        xt(k2) = tn(k2)
        if(xq(k2) <= 0.0) xq(k2) = 1.e-08
     endif
     
     if(ierr == 0)then

!--- calculate moist static energy, heights, qes
        call cup_env(k2, xt, po, t00, 2, z1, psur, xq, xhe, xhes, xqes, xz) 

!--- Environmental values on cloud levels
        call cup_env_clev(k2, psur, z1, xq, xq_cup, xz, xz_cup, po, po_cup,   &
             xt, xt_cup, xqes, xqes_cup, xhes, xhes_cup, xhe, xhe_cup,   &
             gamma_cup)

!**************************** Static Control ************

!--- Moist static energy inside cloud

        xhkb = xhe(k22)

        call cup_up_he(k2, k22, xhkb, xz_cup, cd, mentr_rate,   &
             xhe_cup, xhc, kbcon, xdby, xhe, xhes_cup)
     endif

!--- Normalized mass flux profile

     call cup_up_nms(xzu, xz_cup, mentr_rate, cd, kbcon, ktop, k2, ierr, k22)

!--- Moisture updraft
     
     call cup_up_moisture(k2, ierr, xz_cup, xqc, xqrc, xpw, xpwav, &
          kbcon, ktop, cd, xdby, mentr_rate, xq, gamma_cup, xzu,   &
          xqes_cup, k22, xq_cup)

!--- Workfunctions for updraft
     call cup_up_aa0(k2, xaa0, xz, xzu, xdby, gamma_cup, xt_cup,   &
          kbcon, ktop, ierr) 

!srf-feb-2003
     if(ierr == 0) xaa0_ens(nens) = xaa0 

  enddo

!--------- LARGE SCALE FORCING  -----------------------------------

  call cup_forcing_ens_shal(ensdim, k2, aa0, aa1, xaa0_ens, mbdt_ens,   &
       dtime, xmb, ierr, xf_ens, maxens, maxens2, maxens3, p_cup,  &
       ktop, icoic, iedt, tscl_KF)
  
  if(ierr == 0)then
     do k = 1, k2
        dellat_ens(k,iedt) =  dellat(k)
        dellaq_ens(k,iedt) =  dellaq(k)
     enddo
  else 
     do k = 1, k2
        dellat_ens(k,iedt) = 0.0
        dellaq_ens(k,iedt) = 0.0
     enddo
  endif

!------------------------------ FEEDBACK ------------------------------------

  call cup_output_ens_shal(k2, xf_ens, ierr, dellat_ens, dellaq_ens, &
       outt, outq, xmb, ktop, 'shallow', maxens2, maxens, maxens3,  &
       ensdim)

  return
end subroutine CUP_enss_shal

!======================================================================

subroutine cup_dellas_shallow(k2, ierr, z_cup, p_cup, he, della,   &
     zu, cd, hc, ktop, k22, kbcon, mentr_rate, he_cup)

  use consts_coms, only: grav

  implicit none

  integer, intent(in) :: k2
  integer, intent(in) :: ierr
  real, dimension(k2), intent(in) :: z_cup
  real, dimension(k2), intent(in) :: p_cup
  real, dimension(k2), intent(in) :: he
  real, dimension(k2), intent(out) :: della
  real, dimension(k2), intent(in) :: zu
  real, dimension(k2), intent(in) :: cd
  real, dimension(k2), intent(in) :: hc
  integer, intent(in) :: ktop
  integer, intent(in) :: k22
  integer, intent(in) :: kbcon
  real, intent(in) :: mentr_rate
  real, dimension(k2), intent(in) :: he_cup
  integer :: k
  real :: dz
  real :: subin
  real :: entup
  real :: detup
  real :: subdown
  real :: entupk
  real :: detupk
  real :: totmas
  real :: dp

  della(1:k2) = 0.0

  do k = 2, k2 - 1

     if(ierr /= 0 .or. k > ktop)cycle

!--- Specify detrainment of downdraft, has to be consistent
!--- with zd calculations in soundd.

     dz = z_cup(k+1) - z_cup(k)
     subin = zu(k+1)
     entup = 0.
     detup = 0.
     if(k >= kbcon .and. k < ktop)then
        entup = mentr_rate * dz * zu(k)
        detup = cd(k+1) * dz * zu(k)
     endif
     subdown = zu(k)
     entupk = 0.
     detupk = 0.
     
     if(k == (k22 - 1)) entupk = zu(k22)
     
     if(k == ktop)then
        detupk = zu(ktop)
        subin = 0.
     endif
     
     if(k < kbcon) detup = 0.

!--- Changed due to subsidence and entrainment

     totmas = subin - subdown + detup - entup - entupk + detupk
     
     dp = 100. * (p_cup(k) - p_cup(k+1) )
     della(k)=(subin * he_cup(k+1) - subdown * he_cup(k) + detup * .5 *   &
          (hc(k+1)+ hc(k)) - entup * he(k) - entupk * he_cup(k22) +  &
          detupk * hc(ktop)) * grav / dp
     
  enddo

  return
end subroutine cup_dellas_shallow

!=====================================================================

subroutine cup_forcing_ens_shal(ensdim, k2, aa0, aa1, xaa0, mbdt,   &
     dtime, xmb, ierr, xf, maxens, maxens2, maxens3, p_cup,  &
     ktop, icoic, iedt, tscl_KF)
  
  implicit none

  integer, intent(in) :: maxens
  integer, intent(in) :: maxens2
  integer, intent(in) :: maxens3
  integer, intent(in) :: ensdim
  integer, intent(in) :: k2
  integer, intent(in) :: ktop
  integer, intent(in) :: iedt
  integer, intent(in) :: icoic
  integer, intent(in) :: ierr

  real, intent(in) :: aa0
  real, intent(in) :: aa1
  real, intent(in) :: xaa0(maxens)
  real, intent(in) :: mbdt(maxens)
  real, intent(in) :: dtime
  real, intent(in) :: tscl_KF
  real, intent(in) :: p_cup(k2)

  real, intent(out) :: xmb
  real, intent(out) :: xf(ensdim)

  integer, parameter :: mkxcrt = 25
  real, dimension(mkxcrt), parameter :: pcrit = (/ &
       850., 837.5, 825., 812.5, 800., 787.5, &
       775., 762.5, 750., 737.5, 725., 712.5, &
       700., 687.5, 675., 662.5, 650., 637.5, &
       625., 612.5, 600., 550.,  500., 450.,  &
       400. /)
  real, dimension(mkxcrt), parameter :: acrit = (/   &
       6.323E-02,  5.795E-02,  5.390E-02, 5.236E-02, &
       4.450E-02,  4.965E-02, 5.000E-02,  4.983E-02, &
       5.530E-02,  5.289E-02, 6.080E-02,  5.883E-02, &
       6.640E-02,  6.766E-02, 7.070E-02,  7.937E-02, &
       7.500E-02,  9.396E-02, 0.108,      0.111,     &
       0.130,      0.152,     0.221,      0.315,     &
       0.368 /)
  real, dimension(mkxcrt), parameter :: acritt = (/ &
       0.203,  0.299, 0.359,  0.403, 0.515,  0.478, &
       0.518,  0.530, 0.521,  0.565, 0.543,  0.588, &
       0.566,  0.602, 0.596,  0.611, 0.625,  0.619, &
       0.645,  0.627, 0.665,  0.659, 0.688,  0.743, &
       0.813 /)

  real    :: aclim1
  real    :: aclim2
  real    :: aclim3
  real    :: aclim4
  integer :: kclim
  integer :: k
  real    :: xff0
  real    :: xff_ens3(maxens3)
  integer :: nens
  real    :: xk(maxens)
  integer :: iresultd
  integer :: iresulte
  integer :: ne
  integer :: nall
  integer :: nens3

!--- LARGE SCALE FORCING

  xmb = 0.
  if(ierr == 0)then

     kclim = 0
     check_pcrit_shal: do
        do k = mkxcrt, 1, -1
           if(p_cup(ktop) < pcrit(k))then
              kclim = k
              exit check_pcrit_shal
           endif
        enddo
        if(p_cup(ktop) > pcrit(1)) kclim = 1
        exit check_pcrit_shal
     enddo check_pcrit_shal

     k = max(kclim - 1,1)
     aclim1 = acrit(kclim) * 1.e3
     aclim2 = acrit(k) * 1.e3
     aclim3 = acritt(kclim) * 1.e3
     aclim4 = acritt(k) * 1.e3

!---- Grell's closure
     
     xff0 =  (aa1 - aa0) / dtime
     xff_ens3(1) = xff0
     xff_ens3(2) = .9 * xff_ens3(1)
     xff_ens3(3) = 1.1 * xff_ens3(1)
        
!--- More original Arakawa-Schubert 
        
     xff_ens3(4) = max(0., (aa1 - aclim1) / dtime)
     xff_ens3(5) = max(0., (aa1 - aclim2) / dtime)
     xff_ens3(6) = max(0., (aa1 - aclim3) / dtime)
     xff_ens3(7) = max(0., (aa1 - aclim4) / dtime)

!--- More like Fritsch Chappel or Kain Fritsch (plus triggers)
        
     xff_ens3(8) = aa1 / tscl_KF
     xff_ens3(9) = aa1 / tscl_KF
     xff_ens3(10) = aa1 / tscl_KF
                
     do nens=1,maxens
        xk(nens) = (xaa0(nens) - aa1) / mbdt(nens)
        if(xk(nens) <= 0 .and. xk(nens) > -1.e-9)   &
             xk(nens) = -1.e-9
        if(xk(nens) > 0 .and. xk(nens) < 1.e-9)    &
             xk(nens) = 1.e-9
     enddo

!--- Add up all ensembles
        
     do ne = 1, maxens

!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!

!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive

!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....

        nall = (iedt - 1) * maxens3 * maxens + (ne - 1) * maxens3

!--- Special treatment for stability closures

        if(xff0 > 0 .and. xk(ne) < 0.)then
           xf(nall+1) = max(0., -xff_ens3(1) / xk(ne))
           xf(nall+2) = max(0., -xff_ens3(2) / xk(ne))
           xf(nall+3) = max(0., -xff_ens3(3) / xk(ne))
        endif
           
        if(xk(ne) < 0.)then
           xf(nall+4) = max(0., -xff_ens3(4) / xk(ne))
           xf(nall+5) = max(0., -xff_ens3(5) / xk(ne))
           xf(nall+6) = max(0., -xff_ens3(6) / xk(ne))
           xf(nall+7) = max(0., -xff_ens3(7) / xk(ne))
           xf(nall+8) = max(0., -xff_ens3(8) / xk(ne))
           xf(nall+9) = max(0., -xff_ens3(9) / xk(ne))
           xf(nall+10)= max(0., -xff_ens3(10)/ xk(ne))
        endif

        if(icoic >= 1)then
           do nens3 = 1, maxens3
              xf(nall+nens3) = xf(nall+icoic)
           enddo
        endif
     enddo

  elseif(ierr /= 20 .and. ierr /= 0)then
     xf(1:ensdim) = 0.0
  endif

  return
end subroutine cup_forcing_ens_shal

!================================================================

subroutine cup_output_ens_shal(k2, xf, ierr, dellat, dellaq,     &
     outt, outq, xmb, ktop, name, maxens2, maxens, maxens3,   &
     ensdim)

  implicit none

  integer, intent(in) :: ktop
  integer, intent(in) :: k2
  integer, intent(in) :: ensdim
  integer, intent(in) :: maxens
  integer, intent(in) :: maxens2
  integer, intent(in) :: maxens3

  real, intent(in) :: xf(ensdim)
  real, intent(in) :: dellat(k2,maxens2)
  real, intent(in) :: dellaq(k2,maxens2)

  character(*), intent(in)    :: name
  integer,      intent(inout) :: ierr

  real, intent(out) :: outt(k2)
  real, intent(out) :: outq(k2)
  real, intent(out) :: xmb

  real    :: xmb_ave
  real    :: xmb_std
  real    :: ddtes
  integer :: k
  real    :: dtt
  real    :: dtq
  integer :: n
  real    :: outtes
  integer :: ncount

  outt(1:k2)=0.
  outq(1:k2)=0.
  xmb = 0.
  
!--- calculate ensemble average mass fluxes

  ncount = 0
  if(ierr == 0)then
     do n = 1, ensdim
        if(xf(n) > 0.0)then
           xmb = xmb + xf(n)
           ncount = ncount + 1
        endif
     enddo
     if(ncount > 0)then
        xmb = xmb / real(ncount)
        xmb = min(xmb, 1.0)
     else
        xmb = 0.0
        ierr = 13
     endif
  endif

!-- now do feedback

  ddtes = 250.0
  if(trim(name) == 'shallow')ddtes = 500.0

  do k = 1, k2

     dtt = 0.0
     dtq = 0.0

     if(ierr == 0 .and. k <= ktop)then
        do n = 1, maxens2
           dtt = dtt + dellat(k,n)
           dtq = dtq + dellaq(k,n)
        enddo
        outtes = dtt * xmb * 86400.0 / real(maxens2)

        if(outtes > 2. * ddtes .and. k > 2) then
           xmb = 2. * ddtes / outtes * xmb
           outtes = ddtes
        endif
        if(outtes < -ddtes) then
           xmb = -ddtes / outtes * xmb
           outtes = -ddtes
        endif
        if(outtes > .5 * ddtes .and. k <= 2) then
           xmb = ddtes / outtes * xmb
           outtes = .5 * ddtes
        endif
        
        outt(k)= xmb *dtt / real(maxens2)
        outq(k)= xmb *dtq / real(maxens2)
        
     endif
        
  enddo

  return
end subroutine cup_output_ens_shal
