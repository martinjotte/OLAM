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
subroutine cuparth(iw,dtime)

  use mem_tend,    only: thilt, sh_wt
  use mem_cuparm,  only: thsrc, rtsrc, conprr
  use mem_ijtabs, only: itab_w
  use mem_grid, only: mza, zm, lpw, lpu, lpv, aru, arv, zt,  &
                      unx, uny, unz, vnx, vny, vnz, xew, yew, zew
  use grell_coms, only: iact_gr
  use consts_coms, only: grav, p00, rocp, cp, rdry, erad
  use mem_basic, only: umc, vmc, wmc, uc, vc, wc, theta, press, rho, sh_v
  use grell_deep_coms, only: ensdim
  use misc_coms,  only: io6, meshtype

  implicit none

  integer, intent(in) :: iw
  real, intent(in) :: dtime

  integer :: npoly
  integer :: k2
  integer :: mrl
  integer :: j
  integer :: jv,iv
  real :: mconv
  real :: flux(7),flux0
  real :: z1
  real :: psur
  real :: polyi
  real :: vx, vy, vz
  real :: raxis

  integer :: k
  integer :: kr
  integer :: kv

  real :: exner
  real :: cpdTdt

  real, dimension(mza) :: outt
  real, dimension(mza) :: outq
  real, dimension(mza) :: outqc
  real, dimension(mza) :: po
  real, dimension(mza) :: us
  real, dimension(mza) :: vs
  real, dimension(mza) :: omeg
  real, dimension(mza) :: t
  real, dimension(mza) :: q
  real, dimension(mza) :: tn
  real, dimension(mza) :: qo
  real, dimension(mza) :: p

  integer :: jupwind, iw_upwind
  integer :: iconv_upwind

  real :: dq
  real :: pret

!--- prepare input, erase output

  mconv = 0.0
  flux = 0.0
  thsrc(:, iw) = 0.0
  rtsrc(:, iw) = 0.0
  conprr(iw) = 0.0

! Number of polygon edges for current IW

  npoly = itab_w(iw)%npoly
  polyi = 1. / real(npoly)

!------- transfere valores do rams para o eschema

  z1 = zm(lpw(iw)-1)

! surface pressure in mb

  psur = 0.01 * sngl(press(lpw(iw),iw))
        
! number of vertical levels used in Grell scheme

  k2 = mza - lpw(iw)

  do kr = lpw(iw), mza - 1

! nivel k da grade do grell corresponde ao nivel k + 1 do rams

     k = kr - lpw(iw) + 1
           
     po(k) = sngl(press(kr,iw)) * 1.0e-2 ! pressure in mbar
           
     omeg(k) = - grav * 0.5 * (wmc(kr,iw) + wmc(kr-1,iw))

     t(k) = theta(kr,iw) * (sngl(press(kr,iw)) / p00)**rocp
     q(k) = sh_v(kr,iw)

     tn(k) = t(k) + thilt(kr,iw) * (sngl(press(kr,iw)) / p00)**rocp *   &
             dtime / sngl(rho(kr,iw))
     qo(k) = q(k) + sh_wt(kr,iw) * dtime / sngl(rho(kr,iw))
     p(k) = po(k)

     if (tn(k) < 200.0)tn(k) = t(k)
     qo(k) = max(qo(k), 1.0e-8)
     outt(k) = 0.0
     outq(k) = 0.0
     outqc(k) = 0.0

! Zero out EARTH components of wind and flux array

     vx = 0.
     vy = 0.
     vz = 0.
     
     flux(1:7) = 0.

! Loop over V neighbors of W

     do jv = 1,npoly

! Avoid underground levels for IU or IV point   
   
        if (meshtype == 1) then
           iv = itab_w(iw)%iu(jv)
           kv = max(kr,lpu(iv))
        else
           iv = itab_w(iw)%iv(jv)
           kv = max(kr,lpv(iv))
        endif
           
! Sum EARTH wind components

        vx = vx + uc(kv,iv) * unx(iv) + vc(kv,iv) * vnx(iv)
        vy = vy + uc(kv,iv) * uny(iv) + vc(kv,iv) * vny(iv)
        vz = vz + uc(kv,iv) * unz(iv) + vc(kv,iv) * vnz(iv)
                 
! If p(k) is at least 150 mb less than surface pressure and is also 
! greater than 300 mb, add this k-layer contribution to flux(jv)

        if ((psur - p(k)) > 150.0 .and. p(k) > 300.0) then

           if (meshtype == 1) then
              flux(jv) = flux(jv) + itab_w(iw)%diru(jv) * umc(k,iv) * aru(k,iv)
           else
              flux(jv) = flux(jv) + itab_w(iw)%dirv(jv) * vmc(k,iv) * arv(k,iv)
           endif

        endif

     enddo
     
! Get average EARTH wind components     
                        
     vx = vx / real(npoly)
     vy = vy / real(npoly)
     vz = vz / real(npoly)

! Compute zonal and meridional wind components

     raxis = sqrt(xew(iw) ** 2 + yew(iw) ** 2)  ! dist from earth axis

     if (raxis > 1.e3) then
        us(k) = (vy * xew(iw) - vx * yew(iw)) / raxis

        vs(k) = vz * raxis / erad  &
              - (vx * xew(iw) + vy * yew(iw)) * zew(iw) / (raxis * erad) 
     else
        us(k) = vx
        vs(k) = vy
     endif

  enddo

! Determine which is the upwind cell.

  flux0 = -1.e12
  jupwind = 1

  do jv = 1,npoly

     if (flux0 < flux(jv)) then
        flux0 = flux(jv)
        jupwind = jv
     endif

  enddo
    
  iw_upwind = itab_w(iw)%iw(jupwind)
     
  iconv_upwind = iact_gr(iw_upwind)

  do k = 2, k2 - 1
     dq = 0.5 * (q(k+1) - q(k-1))

!convergencia de umidade da coluna (omega em pa/s)

     mconv = mconv + omeg(k) * dq / grav
  enddo
  if (mconv < 0.0) mconv = 0.0

!---  cumulus parameterization

  call cup_enss(iact_gr(iw), iconv_upwind,  &
     mconv, mza, k2, dtime, outt,   &
     outq, outqc, pret, t, p, z1, psur, q, tn, po, qo, us, vs, omeg)

!----------------------------- output-----------------------
        
  if (pret > 0.0)then
     do k = 1, k2
        kr = k + lpw(iw) - 1

! converte tendencia da temperatura (outt) em tendencia de 
! theta (outtem)
! cp*t=pi*theta => cp dt/dt = theta*dpi/dt + pi*dtheta/dt,
! assumindo dpi/dt (=pt(kr,i,j)) << (exner/theta)*dtheta/dt:

        exner = cp * t(k) / theta(kr,iw)

! tendencia do theta  devida aos cumulus

        thsrc(kr,iw) = cp / exner * outt(k)

! tendencia do rtotal devida aos cumulus

        rtsrc(kr,iw) = outq(k) + outqc(k)

     enddo
     conprr(iw) = pret
     iact_gr(iw) = 1
  else
     iact_gr(iw) = 0
  endif

  return
end subroutine cuparth

!====================================================================

subroutine cup_enss(iact_gr, iconv_upwind,  &
     mconv, mza, k2, dtime, outt, outq, outqc, pre, t, p,  &
     z1, psur, q, tn, po, qo, us, vs, omeg)

  use grell_deep_coms, only: entr_rate, cap_maxs, maxens, maxens2, ensdim,  &
       zkbmax, z_detr, cap_max_increment, zcutdown, depth_min, maxens3,  &
       edtmin, edtmax
  use consts_coms, only: grav, cpi, alvl, t00

  implicit none

  integer, intent(in) :: k2
  integer, intent(in) :: iact_gr
  integer, intent(in) :: iconv_upwind
  real, intent(in) :: dtime
  integer, intent(in) :: mza
  real, dimension(k2), intent(out) :: outt
  real, dimension(k2), intent(out) :: outq
  real, dimension(k2), intent(out) :: outqc
  real, intent(out) :: pre
  real, intent(in) :: mconv
  real, intent(in) :: psur
  real, intent(in) :: z1
  real, dimension(mza), intent(in) :: po
  real, dimension(mza), intent(in) :: us
  real, dimension(mza), intent(in) :: vs
  real, dimension(mza), intent(in) :: omeg
  real, dimension(mza), intent(in) :: t
  real, dimension(mza), intent(in) :: tn
  real, dimension(mza), intent(in) :: p
  real, dimension(mza), intent(inout) :: q
  real, dimension(mza), intent(inout) :: qo
  
  integer :: kstabm
  integer :: ierr
  integer :: iresult
  integer :: nens
  integer :: k
  integer :: kbmax
  integer :: kdet
  integer :: k22
  integer :: kstart
  integer :: kbcon
  integer :: kstabi
  integer :: ktop
  integer :: kzdown
  integer :: jmin
  integer :: ki

  integer :: nens3
  integer :: k22x
  integer :: ierr2
  integer :: ierr3
  integer :: kbconx
  integer :: iedt
  integer :: nall

  real :: mentrd_rate
  real :: mentr_rate
  real :: cap_max
  real :: massfld
  real :: hkb
  real :: hkbo
  real :: zktop
  real :: dz
  real :: dh

  real :: xpwav
  real :: aa0
  real :: aa1
  real :: bu
  real :: buo
  real :: pwev
  real :: pwevo
  real :: pwav
  real :: pwavo
  real :: vshear
  real :: sdp
  real :: vws
  real :: edt
  real :: edto
  real :: edtx
  real :: mbdt
  real :: xhkb
  real :: xpwev
  real :: xaa0

  real, dimension(k2) :: cd
  real, dimension(k2) :: cdd
  real, dimension(maxens) :: mbdt_ens
  real, dimension(maxens2) :: edt_ens
  real, dimension(k2,maxens2) :: dellat_ens
  real, dimension(k2,maxens2) :: dellaq_ens
  real, dimension(k2,maxens2) :: dellaqc_ens
  real, dimension(k2,maxens2) :: pwo_ens
  real, dimension(maxens2) :: edtc
  real, dimension(ensdim) :: xf_ens
  real, dimension(ensdim) :: pr_ens
  real, dimension(ensdim) :: outt_ens
  real, dimension(k2) :: he
  real, dimension(k2) :: hes
  real, dimension(k2) :: qes
  real, dimension(k2) :: z
  real, dimension(k2) :: heo
  real, dimension(k2) :: heso
  real, dimension(k2) :: qeso
  real, dimension(k2) :: zo
  real, dimension(k2) :: q_cup
  real, dimension(k2) :: qo_cup
  real, dimension(k2) :: z_cup
  real, dimension(k2) :: zo_cup
  real, dimension(k2) :: p_cup
  real, dimension(k2) :: po_cup
  real, dimension(k2) :: t_cup
  real, dimension(k2) :: tn_cup
  real, dimension(k2) :: qes_cup
  real, dimension(k2) :: qeso_cup
  real, dimension(k2) :: hes_cup
  real, dimension(k2) :: heso_cup
  real, dimension(k2) :: he_cup
  real, dimension(k2) :: heo_cup
  real, dimension(k2) :: gamma_cup
  real, dimension(k2) :: gammao_cup
  real, dimension(k2) :: hc
  real, dimension(k2) :: hco
  real, dimension(k2) :: dby
  real, dimension(k2) :: dbyo
  real, dimension(k2) :: hcdo
  real, dimension(k2) :: zu
  real, dimension(k2) :: zuo
  real, dimension(k2) :: zd
  real, dimension(k2) :: zdo
  real, dimension(k2) :: hcd
  real, dimension(k2) :: dbyd
  real, dimension(k2) :: dbydo

  real, dimension(maxens) :: xaa0_ens
  real, dimension(maxens3) :: xff_ens3
  real, dimension(k2) :: xhes_cup
  real, dimension(k2) :: xhcd
  real, dimension(k2) :: xhc
  real, dimension(k2) :: xpw
  real, dimension(k2) :: xqc
  real, dimension(k2) :: xqrc
  real, dimension(k2) :: xqrcd
  real, dimension(k2) :: xdby
  real, dimension(k2) :: xqes_cup
  real, dimension(k2) :: xhe_cup
  real, dimension(k2) :: xq_cup
  real, dimension(k2) :: xz_cup
  real, dimension(k2) :: xt_cup
  real, dimension(k2) :: xhes
  real, dimension(k2) :: xqes
  real, dimension(k2) :: xz
  real, dimension(k2) :: xhe
  real, dimension(k2) :: xt
  real, dimension(k2) :: xq
  real, dimension(k2) :: qcd
  real, dimension(k2) :: xqcd
  real, dimension(k2) :: scr1
  real, dimension(k2) :: dellaqc
  real, dimension(k2) :: qcdo
  real, dimension(k2) :: pw
  real, dimension(k2) :: pwo
  real, dimension(k2) :: qrcd
  real, dimension(k2) :: qrcdo
  real, dimension(k2) :: pwd
  real, dimension(k2) :: pwdo
  real, dimension(k2) :: dellah
  real, dimension(k2) :: qco
  real, dimension(k2) :: qrco
  real, dimension(k2) :: qc
  real, dimension(k2) :: qrc
  real, dimension(k2) :: xzd
  real, dimension(k2) :: xzu
  real, dimension(k2) :: dellat
  real, dimension(k2) :: dellaq
  real, dimension(k2) :: xpwd

!--- entrainment of mass
  mentrd_rate = entr_rate
  mentr_rate = entr_rate

!--- initial detrainmentrates
  cd(1:k2) = 0.1 * entr_rate
  cdd(1:k2) = 0.0
  kstabm = k2 - 2
  ierr = 0
      
!---  initialize cap_max 
  cap_max = cap_maxs
      
!--- first check for upstream convection 
!  call cup_direction2(0, massflx, massflx_upwind, iact_gr,   &
!       iconv_upwind, iresult, massfld)
  iresult = 0
  massfld = 0.0


!     (increase cap_max if is there upstream convection)
  if(iresult == 1) cap_max = cap_max + 20.0

!--- MBDT parameter 
  do nens = 1, maxens
     mbdt_ens(nens) = dtime * ( (real(nens) - 3.0) * 1.0e-3 + 5.0e-3)
  enddo

  do nens=1,maxens2
     edt_ens(nens) = 0.95 - real(nens) * 0.01
  enddo

!--- environmental conditions, FIRST HEIGHTS
  if(ierr /= 20)then
     do k = 1, ensdim
        xf_ens(k) = 0.0
        pr_ens(k) = 0.0
        outt_ens(k) = 0.0
     enddo
  endif

!--- calculate moist static energy, heights, qes
  call cup_env(k2, t, p, t00, 0, z1, psur, q, he, hes, qes, z) 
  call cup_env(k2, tn, po, t00, 0, z1, psur, qo, heo, heso, qeso, zo) 

!---  environmental values on cloud levels
  call cup_env_clev(k2, psur, z1, q, q_cup, z, z_cup, p, p_cup, t,   &
       t_cup, qes, qes_cup, hes, hes_cup, he, he_cup, gamma_cup)
  call cup_env_clev(k2, psur, z1, qo, qo_cup, zo, zo_cup, po,   &
       po_cup, tn, tn_cup, qeso, qeso_cup, heso, heso_cup, heo, heo_cup,  &
       gammao_cup)
     
  find_kbmax: do k = 1, k2
     if(zo_cup(k) > (zkbmax + z1))then
        kbmax = k
        exit find_kbmax
     endif
  enddo find_kbmax
     
!---  level where detrainment for downdraft starts
  find_kdet: do k = 1, k2
     if(zo_cup(k) > (z_detr + z1))then
        kdet = k
        exit find_kdet
     endif
  enddo find_kdet
      
!---  DETERMINE LEVEL WITH HIGHEST MOIST STATIC ENERGY CONTENT - K22
!     Gr-dec2002
!     here, ,should just take top of PBL for shallow clouds, also maybe for
!     deep clouds in tropics: kstart=level(pbltop)
  kstart = 3
  call maximi(k2, heo_cup, kstart, kbmax, ierr, k22)

  if(ierr == 0 .and. k22 >= kbmax) ierr = 2

  !--- DETERMINE THE LEVEL OF CONVECTIVE CLOUD BASE  - KBCON
  call cup_kbcon(k22, kbcon, heo_cup, heso_cup,   &
       ierr, kbmax, po_cup, cap_max, k2)

!--- Increase detrainment in stable layers
  call minimi(k2, heso_cup, kbcon, kstabm, kstabi, ierr)

  if(ierr == 0)then
     if((kstabm - 1) > kstabi)then
        do k = kstabi, kstabm - 1
           cd(k) = cd(k-1) + 1.5 * entr_rate
           if(cd(k) > 10.0*entr_rate) cd(k) = 10.0 * entr_rate
        enddo
     endif
  endif

!--- Calculate incloud moist static energy
  if(ierr == 0)then
     call cup_up_he(k2, k22, hkb, z_cup, cd, mentr_rate,   &
          he_cup, hc, kbcon, dby, he, hes_cup)
     call cup_up_he(k2, k22, hkbo, zo_cup, cd, mentr_rate,   &
          heo_cup, hco, kbcon, dbyo, heo, heso_cup)
  endif

  call cup_ktop(1, dbyo, kbcon, ktop, k2, ierr)

  kzdown = 0
  if(ierr == 0)then
     zktop = (zo_cup(ktop) - z1) * 0.6
     zktop = min(zktop + z1, zcutdown + z1)
     find_kzdown: do k = 1, k2
        if(zo_cup(k) > zktop)then
           kzdown = k
           exit find_kzdown
        endif
     enddo find_kzdown
  endif

!--- Downdraft originating level - JMIN
  call minimi(k2, heso_cup, k22, kzdown, jmin, ierr)

  if(ierr == 0)then

!--- Check whether it would have buoyancy, if there where
!--- no entrainment/detrainment
     
     jmin_loop: do
        
        if((jmin - 1) < kdet) kdet = jmin - 1
        if(jmin >= (ktop-1)) jmin = ktop - 2
        ki = jmin
        hcdo(ki) = heso_cup(ki)
        dh = 0.0
        do k = ki - 1, 1, -1
           hcdo(k) = heso_cup(jmin)
           dz = zo_cup(k+1) - zo_cup(k)
           dh = dh + dz * (hcdo(k) - heso_cup(k))
           if(dh > 0.0)then
              jmin = jmin - 1
              if(jmin > 3)then
                 cycle jmin_loop
              else
                 ierr = 9
                 exit jmin_loop
              endif
           endif
        enddo
        
        if(jmin <= 3) ierr = 4
        exit jmin_loop
        
     enddo jmin_loop
     
  endif

!--- Must have at least depth_min m between cloud convective base
!    and cloud top.
  if(ierr == 0)then
     if((zo_cup(ktop) - zo_cup(kbcon)) < depth_min) ierr = 6
  endif

!--- Normalized updraft mass flux profile
  call cup_up_nms(zu, z_cup, mentr_rate, cd, kbcon, ktop, k2, ierr, k22)
  call cup_up_nms(zuo, zo_cup, mentr_rate, cd, kbcon, ktop, k2, ierr, k22)

!--- Normalized downdraft mass flux profile,also work on bottom detrainment
!--- in this routine
  call cup_dd_nms(k2, zd, z_cup, cdd, mentrd_rate, jmin, ierr, 0, kdet, z1)
  call cup_dd_nms(k2, zdo, zo_cup, cdd, mentrd_rate, jmin, ierr, 1, kdet, z1)

!--- Downdraft moist static energy
  call cup_dd_he(k2, hes_cup, hcd, z_cup, cdd, mentrd_rate, jmin,   &
       ierr, he, dbyd)
  call cup_dd_he(k2, heso_cup, hcdo, zo_cup, cdd, mentrd_rate, jmin,   &
       ierr, heo, dbydo)

!--- Calculate moisture properties of downdraft
  call cup_dd_moisture(k2, zd, hcd, hes_cup, qcd, qes_cup, &
       pwd, q_cup, z_cup, cdd, mentrd_rate, jmin, ierr, gamma_cup,   &
       pwev, bu, qrcd, q, 2)
  call cup_dd_moisture(k2, zdo, hcdo, heso_cup, qcdo, qeso_cup, &
       pwdo, qo_cup, zo_cup, cdd, mentrd_rate, jmin, ierr, gammao_cup,   &
       pwevo, bu, qrcdo, qo, 1)

!--- Calculate moisture properties of updraft
  call cup_up_moisture(k2, ierr, z_cup, qc, qrc, pw, pwav, &
       kbcon, ktop, cd, dby, mentr_rate, q, gamma_cup, zu, qes_cup,   &
       k22, q_cup)
  call cup_up_moisture(k2, ierr, zo_cup, qco, qrco, pwo, pwavo, &
       kbcon, ktop, cd, dbyo, mentr_rate, qo, gammao_cup, zuo,   &
       qeso_cup, k22, qo_cup)

  !---  Calculate workfunctions for updrafts
  call cup_up_aa0(k2, aa0, z, zu, dby, gamma_cup, t_cup, kbcon, ktop, ierr) 
  call cup_up_aa0(k2, aa1, zo, zuo, dbyo, gammao_cup, tn_cup, kbcon,   &
       ktop, ierr) 

  if(ierr == 0 .and. aa1 == 0.0) ierr = 17
        
!--- Determine downdraft strength in terms of windshear
  edtc(1:maxens2) = 0.0
  call cup_dd_edt(k2, ierr, us, vs, zo, ktop, kbcon, edt, po, pwavo,   &
       pwevo, edtmax, edtmin, edtc, vshear, sdp, vws)

  do iedt = 1, maxens2

     if(ierr == 0)then
        edt = edtc(iedt)
        edto = edtc(iedt)
        edtx = edtc(iedt)
     endif

     do k = 1, k2
        dellat_ens(k,iedt) = 0.0
        dellaq_ens(k,iedt) = 0.0
        dellaqc_ens(k,iedt) = 0.0
        pwo_ens(k,iedt) = 0.0
     enddo
         
!--- Downdraft workfunctions
! aad is not used.
!     call cup_dd_aa0(mza, k2, edto, ierr, aad, jmin, gammao_cup, tn_cup, &
!          hcdo, heso_cup, zo, zdo)

!---  Change per unit mass that a model cloud would modify the environment

!--- 1. in bottom layer ::: But the lowest della is overwritten in cup_dellas. 
     call cup_dellabot(k2, heo_cup, ierr, zo_cup, po_cup, hcdo, edto, zdo,   &
          cdd, heo, dellah(1), mentrd_rate)
     call cup_dellabot(k2, qo_cup, ierr, zo_cup, po_cup, qrcdo, edto, zdo,   &
          cdd, qo, dellaq(1), mentrd_rate)

!--- 2. everywhere else
     call cup_dellas(k2, ierr, zo_cup, po_cup, hcdo, edto, zdo, cdd,   &
          heo, dellah, mentrd_rate, zuo, cd, hco, ktop, k22, kbcon,   &
          mentr_rate, jmin, heo_cup, kdet, k22)

!--   Take out cloud liquid water for detrainment
     scr1(1:k2) = 0.0
     dellaqc(1:k2) = 0.0
     if(ierr == 0)then
        do k = 1, k2

           scr1(k) = qco(k) - qrco(k)
           if(k == ktop) dellaqc(k) = 0.01 * zuo(ktop) * qrco(ktop) *    &
                grav / (po_cup(k) - po_cup(k+1))
           
           if(k < ktop .and. k > kbcon)then
              dz = zo_cup(k+1) - zo_cup(k)
              dellaqc(k) = 0.01 * grav * cd(k) * dz * zuo(k) * 0.5 *   &
                   (qrco(k) + qrco(k+1)) / (po_cup(k) - po_cup(k+1))
           endif

        enddo
     endif
     
     call cup_dellas(k2, ierr, zo_cup, po_cup, qrcdo, edto, zdo, cdd,   &
          qo, dellaq, mentrd_rate, zuo, cd, scr1, ktop, k22, kbcon,   &
          mentr_rate, jmin, qo_cup, kdet, k22)

!--- Using dellas, calculate changed environmental profiles
     do nens = 1, maxens
        mbdt = mbdt_ens(nens)
        xaa0_ens(nens) = 0.0
        dellat(1:(k2-1)) = 0.0
        if(ierr == 0)then
           do k = 1, k2 - 1
              xhe(k) = dellah(k) * mbdt + heo(k)
              xq(k) = dellaq(k) * mbdt + qo(k)
              dellat(k) = cpi * (dellah(k) - alvl * dellaq(k))
              xt(k) = dellat(k) * mbdt + tn(k)
              if(xq(k) <= 0.0) xq(k) = 1.0e-8
           enddo
           xhe(k2) = heo(k2)
           xq(k2) = qo(k2)
           xt(k2) = tn(k2)
           if(xq(k2) <= 0.0) xq(k2) = 1.0e-8

!--- Calculate moist static energy, heights, qes
           call cup_env(k2, xt, po, t00, 2, z1, psur, xq, xhe, xhes,   &
                xqes, xz)
         
!--- Environmental values on cloud levels
           call cup_env_clev(k2, psur, z1, xq, xq_cup, xz, xz_cup, po,   &
                po_cup, xt, xt_cup, xqes, xqes_cup, xhes, xhes_cup, xhe,   &
                xhe_cup, gamma_cup)

!**************************** Static Control

!--- Moist static energy inside cloud
           xhkb = xhe(k22)
           call cup_up_he(k2, k22, xhkb, xz_cup, cd, mentr_rate,   &
                xhe_cup, xhc, kbcon, xdby, xhe, xhes_cup)
        endif

!--- Normalized mass flux profile
        call cup_up_nms(xzu, xz_cup, mentr_rate, cd, kbcon, ktop, k2,   &
             ierr, k22)
        call cup_dd_nms(k2, xzd, xz_cup, cdd, mentrd_rate, jmin, ierr, &
             1, kdet, z1)
     
!--- Moisture downdraft
        call cup_dd_he(k2, xhes_cup, xhcd, xz_cup, cdd, mentrd_rate,   &
             jmin, ierr, xhe, dbyd)
        call cup_dd_moisture(k2, xzd, xhcd, xhes_cup, xqcd, xqes_cup, &
             xpwd, xq_cup, xz_cup, cdd, mentrd_rate, jmin, ierr, gamma_cup,   &
             xpwev, bu, xqrcd, xq, 3)

!--- Moisture updraft
        call cup_up_moisture(k2, ierr, xz_cup, xqc, xqrc, xpw, xpwav, &
             kbcon, ktop, cd, xdby, mentr_rate, xq, gamma_cup, xzu,   &
             xqes_cup, k22, xq_cup)
        
!--- Workfunctions for updraft
        call cup_up_aa0(k2, xaa0, xz, xzu, xdby, gamma_cup, xt_cup,   &
             kbcon, ktop, ierr) 

!------------------------- LOOP at  CAP_MAX ENSEMBLE ----------
        if(ierr == 0)then

           xaa0_ens(nens) = xaa0
           nall = (iedt-1) * maxens3 * maxens + (nens-1) * maxens3

           do k = 1, k2
              if(k <= ktop)then
                 do nens3 = 1, maxens3
                    if(nens3 == 7)then
                       pr_ens(nall+nens3) = pr_ens(nall+nens3) + pwo(k) +   &
                            edto * pwdo(k)
                    elseif(nens3 == 8)then
                       pr_ens(nall+nens3) = pr_ens(nall+nens3) + pwo(k)
                    elseif(nens3 == 9)then
                       pr_ens(nall+nens3) = pr_ens(nall+nens3) + pwo(k) +  &
                            0.5 * edto * pwdo(k)
                    else
                       pr_ens(nall+nens3) = pr_ens(nall+nens3) + pwo(k) +  &
                            edto * pwdo(k)
                    endif
                 enddo
              endif
           enddo

           
           do nens3=1,maxens3
              outt_ens(nall+nens3) = dellat(1)
           enddo

        endif

     enddo

!--- LARGE SCALE FORCING

!------- CHECK wether aa0 should have been zero
!Gr-dec2002
!   here, ,should just take top of PBL for shallow clouds, also maybe for
!   deep clouds in tropics: kstart=level(pbltop)
     kstart = 3
     call maximi(k2, he_cup, kstart, kbmax, ierr, k22x)
     if(ierr == 0)then
        if(k22x >= kbmax) ierr = 998
     endif

     call cup_forcing_ens_16(k2, aa0, aa1, xaa0_ens, mbdt_ens, dtime,   &
          ierr, xf_ens, 'deeps', iedt, mconv, omeg, kbcon, pr_ens,   &
          edto, kdet, iact_gr, iconv_upwind,   &
          massfld, iresult, xff_ens3, p_cup, ktop, ierr2, ierr3)

     if(ierr == 0)then
        do k = 1, k2
           dellat_ens(k,iedt) = dellat(k)
           dellaq_ens(k,iedt) = dellaq(k)
           dellaqc_ens(k,iedt) = dellaqc(k)
           pwo_ens(k,iedt) = pwo(k) + edt * pwdo(k)
        enddo
     else
        do k = 1, k2
           dellat_ens(k,iedt) = 0.0
           dellaq_ens(k,iedt) = 0.0
           dellaqc_ens(k,iedt) = 0.0
           pwo_ens(k,iedt) = 0.0
        enddo
     endif
  enddo

!---------------------- FEEDBACK   ---------------------------
  call cup_output_ens(k2, xf_ens,ierr, dellat_ens, dellaq_ens,   &
       dellaqc_ens, outt, outq, outqc, pre, pwo_ens, ktop, pr_ens,   &
       outt_ens)

  pre = max(pre, 0.0)

!---------------------------done Cumulus scheme -------------------------

!-------Salva parametros nas analises do RAMS
            
  return
end subroutine cup_enss

!=======================================================================

subroutine cup_dellabot(k2, he_cup, ierr, z_cup, p_cup, hcd, edt, zd,   &
     cdd, he, della, mentrd_rate)

  use consts_coms, only: grav, cp

  implicit none

  integer, intent(in) :: k2
  real, dimension(k2), intent(in) :: he_cup
  integer, intent(in) :: ierr
  real, dimension(k2), intent(in) :: z_cup
  real, dimension(k2), intent(in) :: p_cup
  real, dimension(k2), intent(in) :: hcd
  real, intent(in) :: edt
  real, dimension(k2), intent(in) :: zd
  real, dimension(k2), intent(in) :: cdd
  real, dimension(k2), intent(in) :: he
  real, intent(out) :: della
  real, intent(in) :: mentrd_rate
  real :: dz
  real :: dp
  real :: detdo1
  real :: detdo2
  real :: entdo
  real :: subin
  real :: detdo

  della =0.0
  if(ierr == 0)then
     dz = z_cup(2) - z_cup(1)
     dp = 100.0 * (p_cup(1) - p_cup(2))
     detdo1 = edt * zd(2) * cdd(1) * dz
     detdo2 = edt * zd(1)
     entdo = edt * zd(2) * mentrd_rate * dz
     subin = -edt * zd(2)
     detdo = detdo1 + detdo2 - entdo + subin
     della = (detdo1 * 0.5 * (hcd(1) + hcd(2)) + detdo2 * hcd(1) +  & 
          subin * he_cup(2) - entdo * he(1)) * grav / dp
  endif

  return
end subroutine cup_dellabot

!=======================================================================

subroutine cup_dellas(k2, ierr, z_cup, p_cup, hcd, edt, zd, cdd,   &
     he, della, mentrd_rate, zu, cd, hc, ktop, k22, kbcon, mentr_rate,  &
     jmin, he_cup, kdet, kpbl)

  use consts_coms, only: grav

  implicit none

  integer, intent(in) :: k2
  integer, intent(in) :: ierr
  real, dimension(k2), intent(in) :: z_cup
  real, dimension(k2), intent(in) :: p_cup
  real, dimension(k2), intent(in) :: hcd
  real, intent(in) :: edt
  real, dimension(k2), intent(in) :: zd
  real, dimension(k2), intent(in) :: cdd
  real, dimension(k2), intent(in) :: he
  real, dimension(k2), intent(out) :: della
  real, intent(in) :: mentrd_rate
  real, dimension(k2), intent(in) :: zu
  real, dimension(k2), intent(in) :: cd
  real, dimension(k2), intent(in) :: hc
  integer, intent(in) :: ktop
  integer, intent(in) :: k22
  integer, intent(in) :: kbcon
  real, intent(in) :: mentr_rate
  integer, intent(in) :: jmin
  real, dimension(k2), intent(in) :: he_cup
  integer, intent(in) :: kdet
  integer, intent(in) :: kpbl

  integer :: k
  real :: dz
  real :: detdo
  real :: entdo
  real :: subin
  real :: entup
  real :: detup
  real :: subdown
  real :: entdoj
  real :: entupk
  real :: detupk
  real :: dp

  della(2:k2) = 0.0
  if (ierr /= 0)return

  check_ktop: do k = 2, k2 - 1
     if(k > ktop)exit check_ktop

!--- specify detrainment of downdraft, has to be consistent
!--- with zd calculations in soundd.
     dz    = z_cup(k+1) - z_cup(k)
     detdo = edt * cdd(k) * dz * zd(k+1)
     entdo = edt * mentrd_rate * dz * zd(k+1)
     subin = zu(k+1) - zd(k+1) * edt

     if(k >= kbcon .and. k < ktop)then
        entup = mentr_rate * dz *zu(k)
        detup = cd(k+1) * dz * zu(k)
     else
        entup = 0.0
        detup = 0.0
     endif

     subdown = (zu(k) - zd(k) * edt)
     entdoj  = 0.0
     entupk  = 0.0
     detupk  = 0.0
     
     if(k == jmin) entdoj = zd(k) * edt
     if(k == (k22-1)) entupk = zu(kpbl)
     if(k > kdet) detdo = 0.0
     
     if(k == ktop)then
        detupk = zu(ktop)
        subin = 0.0
     endif

     if(k < kbcon) detup = 0.0
     
!--- changed due to subsidence and entrainment
     dp = 100.0 * (p_cup(k) - p_cup(k+1))
     della(k) = (subin * he_cup(k+1) - subdown * he_cup(k) + detup * 0.5 *  &
          (hc(k+1) + hc(k)) + detdo * 0.5 * (hcd(k+1) + hcd(k)) - entup *   &
          he(k) - entdo * he(k) - entupk * he_cup(k22) - entdoj *   &
          he_cup(jmin) + detupk * hc(ktop)) * grav / dp
     
  enddo check_ktop

  return
end subroutine cup_dellas

!==========================================================================

subroutine cup_forcing_ens_16(k2, aa0, aa1, xaa0, mbdt, dtime, ierr,  &
     xf, name, iedt, mconv, omeg, k22, pr_ens, edt, kbcon,   &
     iact_gr, iconv_upwind, massfld, iresult,   &
     xff_ens3, p_cup, ktop, ierr2, ierr3)

  use consts_coms, only: grav
  use grell_deep_coms, only: maxens, maxens3, ensdim

  implicit none

  integer, intent(in) :: k2
  real, intent(inout) :: aa0
  real, intent(in) :: aa1
  real, dimension(maxens), intent(in) :: xaa0
  real, dimension(maxens), intent(in) :: mbdt
  real, intent(in) :: dtime
  integer, intent(inout) :: ierr
  real, dimension(ensdim), intent(inout) :: xf
  character(len=*), intent(in) :: name
  integer, intent(in) :: iedt
  real, intent(in) :: mconv
  real, dimension(k2), intent(in) :: omeg
  integer, intent(in) :: k22
  real, dimension(ensdim), intent(in) :: pr_ens
  real, intent(in) :: edt
  integer, intent(in) :: kbcon
  integer, intent(in) :: iact_gr
  integer, intent(in) :: iconv_upwind
  real, intent(inout) :: massfld
  integer, intent(inout) :: iresult
  real, dimension(maxens3), intent(out) :: xff_ens3
  real, dimension(maxens) :: xk
  real, dimension(k2), intent(in) :: p_cup
  integer, intent(in) :: ktop
  integer, intent(in) :: ierr2
  integer, intent(in) :: ierr3
  

  integer, parameter :: mkxcrt = 15
  integer :: n
  integer :: nens
  integer :: nens3
  integer :: kclim
  integer :: k
  integer :: ne
  integer :: nall
  integer :: iresulte
  integer :: iresultd
  real, dimension(mkxcrt), parameter :: pcrit = (/ 850., 800., 750., 700.,  &
       650., 600., 550., 500., 450., 400., 350., 300., 250., 200., 150. /)
  real, dimension(mkxcrt), parameter :: acrit = (/ .0633, .0445, .0553,   &
       .0664, .075, .1082, .1521, .2216, .3151, .3677, .41, .5255, .7663,  &
       1.1686, 1.6851 /)
  real, dimension(mkxcrt), parameter :: acritt = (/ .203, .515, .521, .566,  &
       .625, .665, .659, .688, .743, .813, .886, .947, 1.138, 1.377, 1.896 /)
  real :: aclim1
  real :: a1
  real :: aclim2
  real :: aclim3
  real :: aclim4
  real :: xff0
  real :: xomg
  
!--- large scale forcing
  if(trim(name) == 'deeps' .and. ierr > 995)then
     aa0 = 0.0
     ierr = 0
  endif
  if(ierr == 0)then

!added for ensemble 3 with dimension = 16
     kclim = 0
     check_pcrit: do
        do k = mkxcrt, 1, -1
           if(p_cup(ktop) < pcrit(k))then
              kclim = k
              exit check_pcrit
           endif
        enddo
        if(p_cup(ktop) > pcrit(1)) kclim = 1
        exit check_pcrit
     enddo check_pcrit
     
     k = max(kclim-1, 1)
     aclim1 = acrit(kclim) * 1.e3
     aclim2 = acrit(k) * 1.e3
     aclim3 = acritt(kclim) * 1.e3
     aclim4 = acritt(k) * 1.e3
     
!--- treatment different for this closure
     if(trim(name) == 'deeps')then
        
        xff0 = (aa1 - aa0) / dtime
        xff_ens3(1) = (aa1 - aa0) / dtime
        xff_ens3(2) = 0.9 * xff_ens3(1)
        xff_ens3(3) = 1.1 * xff_ens3(1)

!     more like brown (1979), or frank-cohen (199?)

!---  omeg is in bar/s, mconv done with omeg in pa/s
        xff_ens3(4) = -omeg(k22) / grav
        xff_ens3(5) = -omeg(kbcon) / grav
        xff_ens3(6) = -omeg(1) / grav
        do k = 2, kbcon - 1
           xomg = -omeg(k) / grav
           if(xomg > xff_ens3(6)) xff_ens3(6) = xomg
        enddo

!--- more like krishnamurti et al.

        xff_ens3(7) = mconv
        xff_ens3(8) = 0.9 * mconv
        xff_ens3(9) = 1.1 * mconv

!--- more like fritsch chappel or kain fritsch (plus triggers)

!srf - changed at dec/2002 - greater timescale instab. removal
!               xff_ens3(10)=aa1/(60.*20.)
!               xff_ens3(11)=aa1/(60.*30.)
!               xff_ens3(12)=aa1/(60.*40.)
        xff_ens3(10) = aa1 / (60.0 * 50.0)
        xff_ens3(11) = aa1 / (60.0 * 60.0)
        xff_ens3(12) = aa1 / (60.0 * 70.0)
        
   
!--- more original arakawa-schubert (climatologic value of aa0)

        xff_ens3(13) = max(0.0, (aa1 - aclim1) / dtime)
!srf- dec/2002 - later xff_ens3(14) will be equal to xff_ens3(13)
        xff_ens3(14) = max(0.0, (aa1 - aclim2) / dtime)
        xff_ens3(15) = max(0.0, (aa1 - aclim3) / dtime)
        xff_ens3(16) = max(0.0, (aa1 - aclim4) / dtime)
                
!srf - dec/2002
! loop at cap_max ensemble, observe that we are using mbdt = mbdt(1)
! change later
! observe that xk(nens) is the same for any nens.
        do nens = 1, maxens
           xk(nens) = (xaa0(nens) - aa1) / mbdt(nens)
           if(xk(nens) <= 0 .and. xk(nens) > -1.0e-9) xk(nens) = -1.0e-9
           if(xk(nens) > 0 .and. xk(nens) < 1.e-9) xk(nens) = 1.0e-9
        enddo

!--- add up all ensembles

! loop at cap_max ensemble, using mbdt = mbdt(1)
        do ne = 1, maxens
           nall = (iedt-1) * maxens3 * maxens + (ne - 1) * maxens3
!     
!          observe the mass flux calculation:
!--------------------------------------------------------------!
!           ne   |     ierr     | mass flux                    !
!           1    |     ierr =0  |  mf1 = xff_ens3 / xk (ne)    !
!           1    |     ierr >0  |  mf1 =  0                    !
!           2    |     ierr2=0  |  mf2 = mf1                   !
!           2    |     ierr2>0  |  mf2 =  0                    !
!           3    |     ierr3=0  |  mf3 = mf1                   !
!           3    |     ierr3>0  |  mf3 =  0                    !
!                                                              !
! xk(ne) is the same for any 'ne'.                             !
!--------------------------------------------------------------!
! if ierr2 > 0 (convection was not permited for that cap_max)
! then equal to zero the mass flux for the second member of the ensemble (maxens)
!           if(ne == 2 .and. ierr2 > 0)then
!              do nens3 = 1, maxens3
!                 xf(nall+nens3) = 0.0
!                 massfln(nall+nens3) = 0.0
!              enddo
!              cycle maxens_loop
!           endif

! if ierr3 > 0 (convection was not permited for that cap_max)
! then equal to zero the mass flux for the third member of the ensemble (maxens)
!           if(ne == 3 .and. ierr3 > 0)then
!              do nens3 = 1, maxens3
!                 xf(nall+nens3) = 0.0
!                 massfln(nall+nens3) = 0.0
!              enddo
!              cycle maxens_loop
!           endif
!
!--- for every xk, we have maxens3 xffs
!--- iens is from outermost ensemble (most expensive!
!
!--- iedt (maxens2 belongs to it)
!--- is from second, next outermost, not so expensive
!
!--- so, for every outermost loop, we have maxens*maxens2*3
!--- ensembles!!! nall would be 0, if everything is on first
!--- loop index, then ne would start counting, then iedt, then iens....
!
!           iresultd = 0
!           iresulte = 0

!--- check for upwind convection
!           call cup_direction2(1, massflx, massflx_upwind, iact_gr,   &
!                iconv_upwind, iresult, massfld)
!           if(xk(ne) < 0.0 .and. xff0 > 0.0) iresultd = 1
!           iresulte = max(iresult, iresultd)

!           if(iresulte == 1)then

!--- special treatment for stability closures
           if(xff0 > 0.0)then
              xf(nall+1) = max(0.0, -xff_ens3(1) / xk(ne)) + massfld
              xf(nall+2) = max(0.0, -xff_ens3(2) / xk(ne)) + massfld
              xf(nall+3) = max(0.0, -xff_ens3(3) / xk(ne)) + massfld
              xf(nall+13) = max(0.0, -xff_ens3(13) / xk(ne)) + massfld 
              xf(nall+14) = max(0.0, -xff_ens3(14) / xk(ne)) + massfld
              xf(nall+15) = max(0.0, -xff_ens3(15) / xk(ne)) + massfld
              xf(nall+16) = max(0.0, -xff_ens3(16) / xk(ne)) + massfld
           else
              xf(nall+1) = massfld
              xf(nall+2) = massfld
              xf(nall+3) = massfld
              xf(nall+13) = massfld
              xf(nall+14) = massfld
              xf(nall+15) = massfld
              xf(nall+16) = massfld
           endif

!--- if iresult.eq.1, following independent of xff0
           xf(nall+4) = max(0.0, xff_ens3(4) + massfld)
           xf(nall+5) = max(0.0, xff_ens3(5) + massfld)
           xf(nall+6) = max(0.0, xff_ens3(6) + massfld)
           xf(nall+7) = max(0.0, xff_ens3(7) / pr_ens(nall+7))
           xf(nall+8) = max(0.0, xff_ens3(8) / pr_ens(nall+8))
           xf(nall+9) = max(0.0, xff_ens3(9) / pr_ens(nall+9))

!for physical initialization================     
!xff_ens3(7) = pcp rate observed
!icoic = 7 
!only for the initial time
!check units for moist convergence
!at first 20 minutes
!============================================
!              a1 = max(1.0e-9, pr_ens(nall+7))
!              xf(nall+7) = max(0.0, xff_ens3(7) / a1)
!              a1 = max(1.0e-9, pr_ens(nall+8))
!              xf(nall+8) = max(0.0, xff_ens3(8) / a1)
!              a1 = max(1.0e-9, pr_ens(nall+9))
!              xf(nall+9) = max(0.0, xff_ens3(9) / a1)
              
           if(xk(ne) < 0.0)then
              xf(nall+10) = max(0.0, -xff_ens3(10) / xk(ne)) + massfld
              xf(nall+11) = max(0.0, -xff_ens3(11) / xk(ne)) + massfld
              xf(nall+12) = max(0.0, -xff_ens3(12) / xk(ne)) + massfld
           else
              xf(nall+10) = massfld
              xf(nall+11) = massfld
              xf(nall+12) = massfld
           endif
              

!==============
!05-12-2002
!srf - forcing 14 is too bad, use the same for 13:
!a&s (14) = a&s (13)
           xf(nall+14) = xf(nall+13)
!==============


        enddo

     endif

  elseif(ierr /= 20 .and. ierr /= 0)then

     do n = 1, ensdim
        xf(n) = 0.0
     enddo

  endif

  return
end subroutine cup_forcing_ens_16

!======================================================================

subroutine cup_output_ens(k2, xf_ens,ierr, dellat, dellaq, dellaqc,   &
     outt, outq, outqc, pre, pw, ktop, pr_ens, outt_ens)
  
  use grell_deep_coms, only: ensdim, maxens, maxens2, maxens3

  implicit none

  integer, intent(in) :: k2
  real, dimension(ensdim), intent(in) :: xf_ens
  integer, intent(inout) :: ierr
  real, dimension(k2,maxens2), intent(in) :: dellat
  real, dimension(k2,maxens2), intent(in) :: dellaq
  real, dimension(k2,maxens2), intent(in) :: dellaqc
  real, dimension(k2), intent(out) :: outt
  real, dimension(k2), intent(out) :: outq
  real, dimension(k2), intent(out) :: outqc
  real, intent(out) :: pre
  real, dimension(k2,maxens2), intent(in) :: pw
  real :: xmb
  integer, intent(in) :: ktop
  real, dimension(ensdim), intent(inout) :: pr_ens
  real, dimension(ensdim), intent(inout) :: outt_ens

  integer :: k
  real :: xfac_for_dn
  integer, parameter :: i_use_stat_prop = 1
  real :: xmb_ave
  real :: xmb_std
  real :: pr_ave
  real :: pr_std
  integer :: ncount
  integer :: n
  real, parameter :: tunning = 0.7
  real, parameter :: ddtes = 250.0
  real :: dtt
  real :: dtq
  real :: dtqc
  real :: dtpw
  real :: outtes

  do k = 1, k2
     outt(k) = 0.0
     outq(k) = 0.0
     outqc(k) = 0.0
  enddo
  pre = 0.0
  xmb = 0.0
  ncount = 0

!--- calculate ensemble average mass fluxes
  if(ierr == 0)then
           
     do n = 1, ensdim
        pr_ens(n) = pr_ens(n) * xf_ens(n)
        outt_ens(n) = outt_ens(n) * xf_ens(n)
        if(xf_ens(n) > 0.0)then
           xmb = xmb + xf_ens(n)
           ncount = ncount + 1
        endif
     enddo

     if(ncount > 0)then
        xmb = xmb / real(ncount)
     else
        xmb = 0.0
        ierr = 13
     endif

  endif

!------------ now do feedback -----------------------
  do k = 1, k2
     dtt = 0.0
     dtq = 0.0
     dtqc = 0.0
     dtpw = 0.0
     if(ierr == 0.and. k <= ktop)then
        do n = 1, maxens2
           dtt = dtt + dellat(k,n)
           dtq = dtq + dellaq(k,n)
           dtqc = dtqc + dellaqc(k,n)
           dtpw = dtpw + pw(k,n)
        enddo
        outtes = dtt * xmb * 86400.0 / real(maxens2)
        
        if((outtes > (2.0 * ddtes) .and. k > 2))then
           xmb = 2.0 * ddtes / outtes * xmb
           outtes = ddtes
        endif
        
        if (outtes < -ddtes) then
           xmb = -ddtes / outtes * xmb
           outtes = -ddtes
        endif
        
        if (outtes > (0.5 * ddtes) .and. k <= 2) then
           xmb = ddtes / outtes * xmb
           outtes = 0.5 * ddtes
        endif
        outt(k) = outt(k) + xmb * dtt  / real(maxens2)
        outq(k) = outq(k) + xmb * dtq  / real(maxens2)
        outqc(k) = outqc(k) + xmb * dtqc / real(maxens2)
        pre = pre + xmb * dtpw / real(maxens2) ! unit : kg[liq water]/(m^2 s)
     endif
  enddo

  return
end subroutine cup_output_ens

!==========================================================================

subroutine massflx_stats(xf_ens, ensdim, maxens, maxens2, maxens3,    &
     xt_ave, xt_std, ierr)

  implicit none

  integer, intent(in) :: ensdim
  real, dimension(ensdim), intent(in) :: xf_ens
  integer, intent(in) :: maxens
  integer, intent(in) :: maxens2
  integer, intent(in) :: maxens3
  real, intent(out) :: xt_ave
  real, intent(out) :: xt_std
  integer, intent(in) :: ierr
  real, dimension(maxens3) :: x_ave
  integer :: num
  integer :: num2
  integer :: k
  integer :: kk
  integer :: iedt

  num = ensdim / maxens3
  num2 = ensdim / maxens
  
  xt_ave = 0.
  xt_std = 0.
  x_ave(1:maxens3) = 0.
  
  if(ierr /= 0)return

  do kk=1,num
     do k=1,maxens3
        !srf- average in maxens and maxens2 ensemble elements for each
        ! maxens3 (closure) ensemble element (the mean for each closure)
        x_ave(k) = x_ave(k) + xf_ens(maxens3 * (kk - 1) + k)
     enddo
  enddo
     
  x_ave(1:maxens3) = x_ave(1:maxens3) / real(num)
  
  !srf- total average in maxens, maxen2 and maxens3 ensemble elements 
  xt_ave = sum(x_ave(1:maxens3)) / real(maxens3)
  
!--- now do std, skewness,curtosis

  do kk = 1, num
     do k = 1, maxens3
        if(x_ave(k) > 0.)then
!srf- stats properties in maxens and maxens2 ensemble elements for each
! maxens3 (closure) ensemble element (for each closure)
           xt_std  = xt_std  &
                +(xf_ens(maxens3 * (kk - 1) + k) - xt_ave)**2
        endif
     enddo
  enddo

  if(xt_std > 0.)xt_std = sqrt(xt_std / real(num * maxens3))

  return
end subroutine massflx_stats

