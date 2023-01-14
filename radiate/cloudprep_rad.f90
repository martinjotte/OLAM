subroutine cloudprep_rad(iw,ka,mcat,jhcat,rhov,rx,emb,ktop)

! This subroutine was developed from parts of subroutine MIC_COPY in
! omic_driv.f90 and subroutines EACH_COLUMN and ENEMB in omic_misc.f90.

! Arrays rx and cx are the bulk mass and number of hydrometeors PER KG OF AIR,
! NOT PER M^3 AS IN THE ORIGINAL MICROPHYSICS SUBROUTINES.

  use micro_coms, only: ncat, jnmb, rxmin, jhabtab, miclevel, iccn, ccnparm, &
                        emb0, emb1, emb2, zfactor_ccn, iccn, pwmas, dmncofi
  use mem_micro,  only: rr_c, rr_r, rr_p, rr_s, rr_a, rr_g, rr_h, rr_d,  &
                        con_c, con_r, con_p, con_s, con_a, con_g, con_h, con_d, &
                        cldnum, ccntyp
  use ccnbin_coms, only: nccntyp
  use mem_basic,  only: tair, rho
  use mem_grid,   only: mza
  use therm_lib,  only: rhovsl_inv
  use mem_cuparm, only: kcutop, kcubot, iactcu, qwcon

  implicit none

  integer, intent(in) :: iw
  integer, intent(in) :: ka

  integer, intent(out) :: mcat              ! # of active hydrom categories
  integer, intent(out) :: jhcat(mza,ncat+1) ! hydrom category table with ice habits
  real,    intent( in) :: rhov(mza)         ! water vapor density [kg_vap/m^3]
  real,    intent(out) :: rx(mza,ncat+1)    ! hydrom bulk spec dens [kg_hyd/kg_air]
  real,    intent(out) :: emb(mza,ncat+1)   ! hydrom mean particle mass [kg/particle]
  integer, intent(out) :: ktop(ncat+1)

  integer :: k, kb
  integer :: ih
  integer :: ns
  integer :: nt
  integer :: iccntyp

  real :: relhum
  real :: tairc
  real :: dmean
  real :: cx

! If miclevel <= 1, there is no condensate of any type in this simulation.

  if (miclevel <= 1) then

     mcat = 0

  elseif (miclevel == 2) then

! If miclevel = 2, cloud water is the only form of condensate that may exist.
! Set mcat to 1

     mcat = 1

! Zero out microphysics scratch arrays for the present iw column

     rx (:,1) = 0.
     emb(:,1) = 0.
     ktop (1) = 1

! In OLAM, with miclevel = 2, cloud number concentration is specified in cldnum.
! Diagnose cloud droplet mean mass.

     do k = ka,mza
        jhcat(k,1) = 1
        if (rr_c(k,iw) * real(rho(k,iw)) >= rxmin(1)) then
           rx   (k,1) = rr_c(k,iw)
           cx         = cldnum(iw) * zfactor_ccn(k)
           emb  (k,1) = max(emb0(1), min(emb1(1), rx(k,1) / cx))
           ktop   (1) = k
        endif
     enddo

  else

! If miclevel = 3, up to 8 forms of condensate may exist.

     mcat = 8

! Zero out microphysics scratch arrays for the present iw column

     rx (:,1:8) = 0.
     emb(:,1:8) = 0.
     ktop (1:8) = 1

! Copy hydrometeor bulk mass and number concentration from main model arrays
! to microphysics column arrays rx and cx

! Cloud water

     if (jnmb(1) >= 1) then
        do k = ka, mza

           ! If cloud bulk density is sufficiently abundant, copy to rx.
           if (rr_c(k,iw) * real(rho(k,iw)) > rxmin(1)) then

              if (jnmb(1) == 5) then
                 cx = con_c(k,iw)
              elseif (iccn == 1) then
                 if (ccnparm > 1.e6) then
                    cx = ccnparm * zfactor_ccn(k)
                 else
                    cx = cldnum(iw) * zfactor_ccn(k)
                 endif
              elseif (nccntyp == 1) then
                 cx = ccntyp(1)%con_ccn(k,iw)
              else
!Bob: this wouldn't compile:   cx = sum(ccntyp(1:nccntyp)%con_ccn(k,iw))
!Bob: next 4 lines are replacement
                 cx = 0.
                 do iccntyp = 1,nccntyp
                    cx = cx + ccntyp(iccntyp)%con_ccn(k,iw)
                 enddo
              endif

              rx (k,1) = rr_c(k,iw)
              emb(k,1) = max(emb0(1), min(emb1(1), rx(k,1) / max(1.e-12,cx)))

              ktop(1) = k
           endif

        enddo
     endif

! Rain

     if (jnmb(2) >= 1) then
        do k = ka, mza

           ! If rain bulk density is sufficiently abundant, copy to rx,
           if (rr_r(k,iw) * real(rho(k,iw)) > rxmin(2)) then
              rx(k,2) = rr_r(k,iw)
              if (jnmb(2) == 2) then
                 emb(k,2) = emb2(2)
              else
                 emb(k,2) = max(emb0(2), min(emb1(2), rx(k,2) / con_r(k,iw)))
              endif
              ktop(2) = k
           endif

        enddo
     endif

! Pristine ice

     if (jnmb(3) >= 1) then
        do k = ka, mza

           ! If pristine ice bulk density is sufficiently abundant, copy to rx.
           if (rr_p(k,iw) * real(rho(k,iw)) > rxmin(3)) then
              rx (k,3) = rr_p(k,iw)
              emb(k,3) = max(emb0(3), min(emb1(3), rx(k,3) / con_p(k,iw)))
              ktop(3) = k
           endif

        enddo
     endif

! Snow

     if (jnmb(4) >= 1) then

        do k = ka, mza

           ! If snow bulk density is sufficiently abundant, copy to rx.
           if (rr_s(k,iw) * real(rho(k,iw)) > rxmin(4)) then
              rx(k,4) = rr_s(k,iw)
              if (jnmb(4) == 2) then
                 emb(k,4) = emb2(4)
              else
                 emb(k,4) = max(emb0(4), min(emb1(4), rx(k,4) / con_s(k,iw)))
              endif
              ktop(4) = k
           endif

        enddo
     endif

! Aggregates

     if (jnmb(5) >= 1) then

        do k = ka, mza

           ! If aggregates bulk density is sufficiently abundant, copy to rx.
           if (rr_a(k,iw) * real(rho(k,iw)) > rxmin(5)) then
              rx(k,5) = rr_a(k,iw)
              if (jnmb(5) == 2) then
                 emb(k,5) = emb2(5)
              else
                 emb(k,5) = max(emb0(5), min(emb1(5), rx(k,5) / con_a(k,iw)))
              endif
              ktop(5) = k
           endif

        enddo
     endif

! Graupel

     if (jnmb(6) >= 1) then

        do k = ka, mza

           ! If graupel bulk density is sufficiently abundant, copy to rx,
           if (rr_g(k,iw) * real(rho(k,iw)) > rxmin(6)) then
              rx(k,6) = rr_g(k,iw)
              if (jnmb(6) == 2) then
                 emb(k,6) = emb2(6)
              else
                 emb(k,6) = max(emb0(6), min(emb1(6), rx(k,6) / con_g(k,iw)))
              endif
              ktop(6) = k
           endif

        enddo
     endif

! Hail

     if (jnmb(7) >= 1) then

        do k = ka, mza

           ! If hail bulk density is sufficiently abundant, copy to rx,
           if (rr_h(k,iw) * real(rho(k,iw)) > rxmin(7)) then
              rx(k,7) = rr_h(k,iw)
              if (jnmb(7) == 2) then
                 emb(k,7) = emb2(7)
              else
                 emb(k,7) = max(emb0(7), min(emb1(7), rx(k,7) / con_h(k,iw)))
              endif
              ktop(7) = k
           endif
        enddo
     endif

! Drizzle

     if (jnmb(8) >= 1) then

        do k = ka, mza

           ! If drizzle bulk density is sufficiently abundant, copy to rx,
           if (rr_d(k,iw) * real(rho(k,iw)) > rxmin(8)) then
              rx (k,8) = rr_d(k,iw)
              emb(k,8) = max(emb0(8), min(emb1(8), rx(k,8) / con_d(k,iw)))
              ktop(8) = k
           endif

        enddo
     endif

! Get microphysics habits

     if (iactcu(iw) > 0) then
        kb = max( maxval(ktop(1:8)), kcutop(iw) )
     else
        kb = maxval(ktop(1:8))
     endif

     do k = ka, kb
        tairc  = tair(k,iw) - 273.15
        relhum = min(1.,rhov(k) * rhovsl_inv(tairc))

        ns = max(1,nint(100. * relhum))
        nt = max(1,min(31,-nint(tairc)))

        jhcat(k,1) = 1
        jhcat(k,2) = 2
        jhcat(k,3) = jhabtab(nt,ns,1)
        jhcat(k,4) = jhabtab(nt,ns,2)
        jhcat(k,5) = 5
        jhcat(k,6) = 6
        jhcat(k,7) = 7
        jhcat(k,8) = 8
     enddo

  endif

! Convective (subgrid) clouds

  mcat = mcat + 1

  rx (:,mcat) = 0.
  emb(:,mcat) = 0.
  ktop (mcat) = 1

  if (iactcu(iw) > 0) then

     ktop(mcat) = kcutop(iw)

     do k = kcubot(iw), kcutop(iw)

        rx(k,mcat) = qwcon(k,iw)
        tairc = tair(k,iw) - 273.15

        if (miclevel /= 3 .or. tairc > -12.) then
           jhcat(k,mcat) = 1
           cx          = cldnum(iw) * zfactor_ccn(k)
           emb(k,mcat) = max(emb0(1), min(emb1(1), rx(k,mcat) / cx))
        else
           ih = jhcat(k,3)
           jhcat(k,mcat) = ih

           ! Mean maximum dimension of ice crystals as a function of T
           ! (see Kristjansson et al., 2000, JGR)
           dmean = 1.0307e-3 * exp(0.05522*(tair(k,iw)-279.5))

           emb(k,mcat) = max(emb0(3), min(emb1(3), (dmean*dmncofi(ih))**pwmas(ih)))
        endif
     enddo

  endif

end subroutine cloudprep_rad
