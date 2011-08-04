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
subroutine harr_swrad(nrad,iw,albedt_beam,albedt_diffuse,cosz,time,   &
                      u,pl,tl,dzl,vp,tp,omgp,gp,fu,fd,flxus,flxds,ngass,  &
                      flx_diff)

! Bob's interface subroutine

use mem_harr, only: mg, ng, mb, nb, nsolb, npsb, xp, alpha, beta, wght,  &
                    prf,trf,ralcs, solar1, ulim

implicit none

integer, intent(in) :: nrad
integer, intent(in) :: iw
integer, intent(in) :: ngass(mg)

real, intent(in) :: albedt_beam
real, intent(in) :: albedt_diffuse
real, intent(in) :: cosz
real, intent(in) :: time
real, intent(in) :: u   (nrad,3)
real, intent(in) :: pl  (nrad)
real, intent(in) :: tl  (nrad)
real, intent(in) :: dzl (nrad)
real, intent(in) :: vp  (nrad)
real, intent(in) :: tp  (nrad,mb)
real, intent(in) :: omgp(nrad,mb)
real, intent(in) :: gp  (nrad,mb)
real, intent(in) :: fu  (nrad,6)
real, intent(in) :: fd  (nrad,6)

real, intent(inout) :: flxus(nrad)
real, intent(inout) :: flxds(nrad)
real :: flx_diff

call swrad(nrad,ng,nb,nsolb,npsb,     &
   u,pl,tl,dzl,vp,                    &
   xp,alpha,beta,wght,prf,trf,ralcs,  &
   solar1,ngass,                      &
   albedt_beam,albedt_diffuse,cosz,                       &
   tp,omgp,gp,fu,fd,flxus,flxds,ulim,iw,time,flx_diff)

return
end subroutine harr_swrad

!===============================================================================

subroutine harr_lwrad(nrad,u,pl,tl,dzl,vp,tp,omgp,gp,fu,fd,flxul,flxdl,ngast)

! Bob's interface subroutine

use mem_harr, only: mg, ng, mb, nb, nsolb, npsb, nuum, xp, alpha, beta, wght,  &
                    prf, trf, ralcs, a0, a1, a2, a3, ulim

implicit none

integer, intent(in) :: nrad
integer, intent(in) :: ngast(mg)

real, intent(in) :: u   (nrad,3)
real, intent(in) :: pl  (nrad)
real, intent(in) :: tl  (nrad)
real, intent(in) :: dzl (nrad)
real, intent(in) :: vp  (nrad)
real, intent(in) :: tp  (nrad,mb)
real, intent(in) :: omgp(nrad,mb)
real, intent(in) :: gp  (nrad,mb)
real, intent(in) :: fu  (nrad,6)
real, intent(in) :: fd  (nrad,6)

real, intent(inout) :: flxul(nrad)
real, intent(inout) :: flxdl(nrad)

call lwrad(nrad,ng,nb,nsolb,npsb,nuum,   &
   u,pl,tl,dzl,vp,                       &
   xp,alpha,beta,wght,prf,trf,ralcs,     &
   a0,a1,a2,a3,                          &
   ngast,                                &
   tp,omgp,gp,fu,fd,flxul,flxdl,ulim)
return
end subroutine harr_lwrad

!===============================================================================

! * New swrad (Jerry, March 8, 1996)

subroutine swrad(nz,ng,nb,ns,npsb,               &
           u,pl,tl,dz,vp,                        &
           xp,alpha,beta,wght,prf,trf,ral,       &
           solar,ngas,                           &
           alb_b,alb_d,amu0,                             &
           tp,omgp,asym,fu,fd,flxsu,flxsd,ulim,i,time,flx_diff)

implicit none
integer :: nz,ng,nb,ns,i
real :: time
real, intent(out) :: flx_diff

!
!     two-stream radiative transfer code.
!     This version will attempt the use of FESFT
!     as outlined in R&G (MWR 1991)
!      JH 5-95
!
integer, parameter :: mb=8,mg=3,mk=7
real, parameter :: top=1800.,tm=1800./296.,gma=0.002,cp=1004.

!     input

real    :: u(nz,mg),pl(nz),tl(nz),dz(nz),vp(nz)
real    :: ulim(mg,mb)
integer :: npsb(mg,mb),na(mg),ngas(mg)

!     input parameters

real    :: ral(mb),wlenlo(mb),wlenhi(mb)
real    :: xp(mg,mk,mb),alpha(mg,mk,mb),beta(mg,mk,mb)
real    :: wght(mg,mk,mb),prf(mg,mb),trf(mg,mb)
real    :: solar(mb)

!     arrays used in fluxes

real    :: tg(nz),tp(nz,nb),tcr(nz),omgp(nz,nb)
real    :: alb_b,alb_d,amu0,asym(nz,nb)
real    :: t(nz),r(nz),tc(nz),sigu(nz),sigd(nz)
real    :: re(nz),vd(nz),td(nz),vu(nz)
real    :: fu(nz,6),fd(nz,6),flxsu(nz),flxsd(nz)

integer ib,iz,mxa,ig,ik,ig1,ig2,ik1,ik2,ig3,ik3
real :: uu,fact,ufac,dfac,tu1,tu2,td1,td2,fbu,fbd,tu3,td3
real, dimension(6) :: vdflx
real :: diffuse_fac,tdiffuse1,tdiffuse2,f_diffuse_d,tdiffuse3

!     output
flx_diff = 0.0
!     remember, for this code it is assumed that u, the gaseous
!     absorber amounts, are specified in Pa
!----------------------------------------------------------------------------
!     zero out flxsu and flxsd

!     do iz = 1,nz
!      flxsu(iz) = 0.0
!      flxsd(iz) = 0.0
!     enddo

!     loop through each band

do ib=1,ns

!        calculate here the properties that are considered grey,
!        i.e. averaged values are used across the band...
!        rayleigh scatter...cloud is now done outside of this routine

!           get rayleigh scattering AND
!           zero out the local flux arrays

   do iz=1,nz
      tcr(iz) = ral(ib)*pl(iz)*dz(iz)/  &
           tl(iz)
      fu(iz,2) = 0.
      fu(iz,3) = 0.
      fu(iz,4) = 0.
      fd(iz,2) = 0.
      fd(iz,3) = 0.
      fd(iz,4) = 0.
   enddo
   vdflx(:) = 0.

!        determine if, and how many overlaps...also check
!        to see if the gas is used

   mxa = 0
   do ig=1,mg
      if (npsb(ig,ib).gt.0.and.ngas(ig).eq.1)then
         mxa = mxa+1
         na(mxa) = ig
      endif
   enddo

   if (mxa.eq.0) go to 111
   if (mxa.eq.1) then

!-------------------------------------------------------------------------
!           no overlapping gasses, single aborber

      ig = na(1)

      do ik=1,npsb(ig,ib)

          do iz=nz,2,-1
            uu=min(ulim(ig,ib),u(iz,ig))
            tg(iz) = xp(ig,ik,ib)*uu*  &
               (pl(iz)/prf(ig,ib))**alpha(ig,ik,ib)*  &
               (trf(ig,ib)/tl(iz))**beta(ig,ik,ib)
          enddo

!               now do rest of stuff in subroutine

          call flxsw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb_b,alb_d,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib),i,time,vdflx(6))

!               add pseudo-band fluxes to the total flux

          do iz=1,nz
            flxsu(iz)=flxsu(iz)+wght(ig,ik,ib)*fu(iz,6)
            flxsd(iz)=flxsd(iz)+wght(ig,ik,ib)*fd(iz,6)
          enddo
          flx_diff = flx_diff + wght(ig,ik,ib) * vdflx(6)
       enddo
   else if (mxa.eq.2) then

!--------------------------------------------------------------------------
!           overlap of two gasses using the FESFT

      ig1 = na(1)
      ig2 = na(2)

!           do the gray fluxes first

          tg(1:nz) = 0.

          call flxsw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb_b,alb_d,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,1),fd(1,1),asym(1,ib),i,time,vdflx(1))

!           do the 1st gas


      do ik1=1,npsb(ig1,ib)

          do iz=2,nz
            uu=min(ulim(ig1,ib),u(iz,ig1))
            tg(iz) = xp(ig1,ik1,ib)*uu*  &
               (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
               (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
          enddo

          call flxsw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb_b,alb_d,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib),i,time,vdflx(1))

          fact = wght(ig1,ik1,ib)
          do iz=1,nz
            fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
            fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
          enddo
          vdflx(2) = vdflx(2) + fact * vdflx(6)

       enddo

!            do the 2nd gas

       do ik2=1,npsb(ig2,ib)

          do iz=1,nz
            uu=min(ulim(ig2,ib),u(iz,ig2))
            tg(iz) = xp(ig2,ik2,ib)*uu*  &
               (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
               (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
          enddo

          call flxsw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb_b,alb_d,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib),i,time,vdflx(6))

          fact = wght(ig2,ik2,ib)
          do iz=1,nz
            fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
            fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
          enddo
          vdflx(3) = vdflx(3) + fact * vdflx(6)
      enddo

!           add together the fluxes: FESFT method

!           do iz=1,nz
!             ufac=max(1.e-14,fu(iz,1))
!             dfac=max(1.e-14,fd(iz,1))
!             fbu=fu(iz,2)/ufac*fu(iz,3)
!    +         /ufac*fu(iz,1)
!             fbd=fd(iz,2)/dfac*fd(iz,3)
!    +         /dfac*fd(iz,1)
!             flxsu(iz)=flxsu(iz)+fbu
!             flxsd(iz)=flxsd(iz)+fbd
!           enddo

! this is the stuff that is fixed (old stuff above)

      do iz=1,nz
        ufac=max(1.e-14,fu(iz,1))
        dfac=max(1.e-14,fd(iz,1))
        tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
        tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
        td1 = max(0.0,min(1.1,fd(iz,2)/dfac))
        td2 = max(0.0,min(1.1,fd(iz,3)/dfac))
        fbu = tu1*tu2*fu(iz,1)
        fbd = td1*td2*fd(iz,1)

        flxsu(iz)=flxsu(iz)+fbu
        flxsd(iz)=flxsd(iz)+fbd
      enddo

      diffuse_fac = max(1.e-14,vdflx(1))
      tdiffuse1 = max(0.0,min(1.1,vdflx(2)/diffuse_fac))
      tdiffuse2 = max(0.0,min(1.1,vdflx(3)/diffuse_fac))
      f_diffuse_d = tdiffuse1 * tdiffuse2 * vdflx(1)
      flx_diff = flx_diff + f_diffuse_d

      else if (mxa.eq.3) then

!--------------------------------------------------------------------------
!           overlap of three gasses

      ig1 = na(1)
      ig2 = na(2)
      ig3 = na(3)

!           do the gray fluxes first

          tg(1:nz) = 0.

          call flxsw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb_b,alb_d,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,1),fd(1,1),asym(1,ib),i,time,vdflx(1))

!           do the 1st gas


      do ik1=1,npsb(ig1,ib)

          do iz=2,nz
            uu=min(ulim(ig1,ib),u(iz,ig1))
            tg(iz) = xp(ig1,ik1,ib)*uu*  &
               (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
               (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
          enddo

          call flxsw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb_b,alb_d,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib),i,time,vdflx(6))

          fact = wght(ig1,ik1,ib)
          do iz=1,nz
            fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
            fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
          enddo
          vdflx(2) = vdflx(2) + fact * vdflx(6)
        enddo

!           do the 2nd gas

      do ik2=1,npsb(ig2,ib)

          do iz=2,nz
            uu=min(ulim(ig2,ib),u(iz,ig2))
            tg(iz) = xp(ig2,ik2,ib)*uu*  &
               (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
               (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
          enddo

          call flxsw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb_b,alb_d,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib),i,time,vdflx(6))

          fact = wght(ig2,ik2,ib)
          do iz=1,nz
            fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
            fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
          enddo
          vdflx(3) = vdflx(3) + fact * vdflx(6)
        enddo

!          do the 3rd gas

      do ik3=1,npsb(ig3,ib)

          do iz=2,nz
            uu=min(ulim(ig3,ib),u(iz,ig3))
            tg(iz) = xp(ig3,ik3,ib)*uu*  &
               (pl(iz)/prf(ig3,ib))**alpha(ig3,ik3,ib)*  &
               (trf(ig3,ib)/tl(iz))**beta(ig3,ik3,ib)
          enddo

!               now do rest of stuff in subroutine

          call flxsw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  alb_b,alb_d,solar(ib),amu0,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),asym(1,ib),i,time,vdflx(6))

!               sum the pseudo-band fluxes to get total flux

          fact = wght(ig3,ik3,ib)
          do iz=1,nz
            fu(iz,4) = fu(iz,4)+fact*fu(iz,6)
            fd(iz,4) = fd(iz,4)+fact*fd(iz,6)
          enddo
          vdflx(4) = vdflx(4) + fact * vdflx(6)
      enddo

!           sum the fluxes together: FESFT method

!           do iz=1,nz
!             ufac=max(1.e-14,fu(iz,1))
!             dfac=max(1.e-14,fd(iz,1))
!             fbu=fu(iz,2)/ufac*fu(iz,3)
!    +         /ufac*fu(iz,4)/ufac*fu(iz,1)
!             fbd=fd(iz,2)/dfac*fd(iz,3)
!    +         /dfac*fd(iz,4)/dfac*fd(iz,1)
!             flxsu(iz)=flxsu(iz)+fbu
!             flxsd(iz)=flxsd(iz)+fbd
!           enddo

! this is the stuff that is fixed (old stuff above)

      do iz=1,nz
        ufac=max(1.e-14,fu(iz,1))
        dfac=max(1.e-14,fd(iz,1))
        tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
        tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
        tu3 = max(0.0,min(1.1,fu(iz,4)/ufac))
        td1 = max(0.0,min(1.1,fd(iz,2)/dfac))
        td2 = max(0.0,min(1.1,fd(iz,3)/dfac))
        td3 = max(0.0,min(1.1,fd(iz,4)/dfac))

        fbu = tu1*tu2*tu3*fu(iz,1)
        fbd = td1*td2*td3*fd(iz,1)
        flxsu(iz)=flxsu(iz)+fbu
        flxsd(iz)=flxsd(iz)+fbd
      enddo
      diffuse_fac=max(1.e-14,vdflx(1))
      tdiffuse1 = max(0.0,min(1.1,vdflx(2)/diffuse_fac))
      tdiffuse2 = max(0.0,min(1.1,vdflx(3)/diffuse_fac))
      tdiffuse3 = max(0.0,min(1.1,vdflx(4)/diffuse_fac))
      f_diffuse_d = tdiffuse1*tdiffuse2*tdiffuse3*vdflx(1)
      flx_diff = flx_diff + f_diffuse_d

!-------------------------------------------------------------------------

   else
      write(*,5000)
5000       format('Two is the maximum amount of overlapping',/  &
            'gasses allowed: if you really want to have ',/  &
            'more gasses overlap, they better have few',/  &
            'pseudo-bands for all of them, and in any',/  &
            'case it will cost you. To do that come into',/  &
            'the code and modify it here, look at previous',/  &
            'structures to see how it is done. Its easy',/  &
            'but BEWARE: ITS BOUND TO BE HORRIBLY EXPENSIVE')
      stop
   endif
111     continue

enddo

return
end subroutine swrad

!===============================================================================

! New lwrad (March 8, 1996)

subroutine lwrad(nz,ng,nb,ns,npsb,nuum,            &
           u,pl,tl,dz,vp,                          &
           xp,alpha,beta,wght,prf,trf,ral,         &
           a0,a1,a2,a3,                            &
           ngas,                                   &
           tp,omgp,asym,fu,fd,flxlu,flxld,ulim)

implicit none

integer, parameter :: mb=8,mg=3,mk=7
real, parameter :: top=1800.,tm=1800./296.,gma=0.002,cp=1004.
integer :: nz,ng,nb,ns
!
!     two-stream radiative transfer code.
!     This version will attempt the use of FESFT
!     as outlined in R&G (MWR 1991)
!      JH 5-95
!
!     input

real    u(nz,mg),pl(nz),tl(nz),dz(nz),vp(nz)
real    ulim(mg,mb)
integer npsb(mg,mb),na(mg),nuum(mb),ngas(mg)

!     input parameters

real    ral(mb),wlenlo(mb),wlenhi(mb)
real    xp(mg,mk,mb),alpha(mg,mk,mb),beta(mg,mk,mb)
real    wght(mg,mk,mb),prf(mg,mb),trf(mg,mb)
real    a0(mb),a1(mb),a2(mb),a3(mb)

!     arrays used in fluxes

real    tg(nz),tp(nz,nb),tcr(nz),omgp(nz,nb),src(nz)
real    t(nz),r(nz),tc(nz),sigu(nz),sigd(nz)
real    re(nz),vd(nz),td(nz),vu(nz),asym(nz,nb)
real    fu(nz,6),fd(nz,6),flxlu(nz),flxld(nz)

integer :: ib,iflag,iz,mxa,ig,ik,ig1,ig2,ik1,ik2,nir,ii,ig3,ik3
real :: tf,ewght,uu,chck,fact,ufac,tu1,tu2,fbu,dx,adx,fdx &
       ,dxc,tn1c,tn2c,tn1,tn2,difflxb,fbd,dfac,tu3,tn3c,tn3

!     remember, for this code it is assumed that u, the gaseous
!     absorber amounts, are specified in g/m^3
!----------------------------------------------------------------------------
!     zero out flxsu and flxsd

!     do iz = 1,nz
!      flxlu(iz) = 0.0
!      flxld(iz) = 0.0
!     enddo

!     loop through each band

nir=ns+1
do ib=nir,nb

!        calculate the properties that are grey across the band
!        Planck function at the interfaces, continuum absorption,
!        ...cloud is now done outside of this routine

      iflag=0

!           do surface and top of atmosphere first

     src(1) = a0(ib)+tl(1)*(a1(ib)+tl(1)*  &
                     (a2(ib)+a3(ib)*tl(1)))
     src(2) = a0(ib)+tl(2)*(a1(ib)+tl(2)*  &
                     (a2(ib)+a3(ib)*tl(2)))
     src(nz) =  a0(ib)+tl(nz)*(a1(ib)+tl(nz)*  &
                         (a2(ib)+a3(ib)*tl(nz)))

     tcr(1)=0.
     tcr(2)=0.
     tcr(nz)=0.

     do iz=nz-1,3,-1

!                 get sources at the interface, temp first

            tf = 0.5*(tl(iz)+tl(iz-1))
            src(iz) = a0(ib)+tf*(a1(ib)+tf*  &
                         (a2(ib)+a3(ib)*tf))
            tcr(iz) = 0.
      enddo


      if (nuum(ib).eq.1) then

!              this band has continuum

         do iz=2,nz
            tcr(iz) = ral(ib)*u(iz,1)*  &
                 exp(top/tl(iz)-tm)*(vp(iz)+  &
                 gma*(pl(iz)-vp(iz)))
         enddo

      endif


   do iz=1,nz
         fu(iz,2) = 0.
         fu(iz,3) = 0.
         fu(iz,4) = 0.
         fd(iz,2) = 0.
         fd(iz,3) = 0.
         fd(iz,4) = 0.
         tg(iz)=0.
   enddo

!        determine if, and how many overlaps

   mxa = 0
   do ig=1,mg
      if (npsb(ig,ib).gt.0.and.ngas(ig).eq.1)then
         mxa = mxa+1
         na(mxa) = ig
      endif
   enddo

   if (mxa.eq.0) go to 111
   if (mxa.eq.1) then

!-------------------------------------------------------------------------
!           no overlapping gasses, single aborber

      ig = na(1)

      do ik=1,npsb(ig,ib)

          do iz=nz,2,-1
            uu=min(ulim(ig,ib),u(iz,ig))
            tg(iz) = xp(ig,ik,ib)*uu*  &
               (pl(iz)/prf(ig,ib))**alpha(ig,ik,ib)*  &
               (trf(ig,ib)/tl(iz))**beta(ig,ik,ib)
          enddo

!               now do rest of stuff in subroutine

          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),1.,asym(1,ib))

!               add pseudo-band fluxes to the total flux

          do iz=1,nz
            flxlu(iz)=flxlu(iz)+wght(ig,ik,ib)*fu(iz,6)
            flxld(iz)=flxld(iz)+wght(ig,ik,ib)*fd(iz,6)
          enddo
      enddo
   else if (mxa.eq.2) then

!--------------------------------------------------------------------------
!           overlap of two gasses using the FESFT

      ig1 = na(1)
      ig2 = na(2)

!           do the gray fluxes first

      chck=float(iflag)
      
          tg(1:nz) = 0.
      
          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,1),fd(1,1),chck,asym(1,ib))

!           if there is continuum abs. do it now since tg=0
!
     if(nuum(ib).eq.1)then

          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,5),fd(1,5),1.,asym(1,ib))
     endif

!           do the 1st gas


      do ik1=1,npsb(ig1,ib)

          do iz=2,nz
            uu=min(ulim(ig1,ib),u(iz,ig1))
            tg(iz) = xp(ig1,ik1,ib)*uu*  &
               (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
               (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
          enddo

          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

          fact = wght(ig1,ik1,ib)
          do iz=1,nz
            fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
            fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
          enddo
       enddo

!            do the 2nd gas

       do ik2=1,npsb(ig2,ib)

          do iz=1,nz
            uu=min(ulim(ig2,ib),u(iz,ig2))
            tg(iz) = xp(ig2,ik2,ib)*uu*  &
               (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
               (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
          enddo

          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

          fact = wght(ig2,ik2,ib)
          do iz=1,nz
            fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
            fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
          enddo
      enddo

!           add together the fluxes: FESFT method...
!          NOTE: Here we use FESFT on the net flux NOT
!                the individual fluxes, as noted in R and G.
!                 FESFT seems not only to work for the net flux
!               but also for the upwelling IR. Thus, we may still use
!              FESFT to determine each flux individually.

!           do iz=1,nz
!             ufac=max(1.e-14,fu(iz,1))
!             dfac=max(1.e-14,fd(iz,1))
!             dx=dfac-ufac
!             fbu=fu(iz,2)/ufac*fu(iz,3)
!    +         /ufac*fu(iz,1)*(fu(iz,5)/ufac
!    +         *float(nuum(ib))+1.-float(nuum(ib)))
!             dx1=(fd(iz,2)-fu(iz,2))/dx
!             dx2=(fd(iz,3)-fu(iz,3))/dx
!             dxc=(fd(iz,5)-fu(iz,5))/dx*float(nuum(ib))
!    +         +1.-float(nuum(ib))
!             difflxb=dx1*dx2*dxc*dx
!             fbd=difflxb+fbu
!             flxlu(iz)=flxlu(iz)+fbu
!             flxld(iz)=flxld(iz)+fbd
!           enddo

!  above is old section, below is the new bounded part:

      do iz=1,nz
        ufac=max(1.e-14,fu(iz,1))
        tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
        tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
        fbu = tu1*tu2*fu(iz,1)

        dx=fd(iz,1)-fu(iz,1)
        adx = abs(dx)
!        fdx = dx/adx
        fdx = sign(1.0,dx)
        dxc = max(1.e-14,adx)*fdx
        tn1c = (fd(iz,2)-fu(iz,2))/dxc + 1.
        tn2c = (fd(iz,3)-fu(iz,3))/dxc + 1.
        tn1 = min(max(tn1c,0.),2.) - 1.
        tn2 = min(max(tn2c,0.),2.) - 1.
        difflxb = tn1*tn2*dx
        fbd=difflxb+fbu

        flxlu(iz)=flxlu(iz)+fbu
        flxld(iz)=flxld(iz)+fbd
      enddo


      else if (mxa.eq.3) then

!--------------------------------------------------------------------------
!           overlap of three gasses

      ig1 = na(1)
      ig2 = na(2)
      ig3 = na(3)

!           do the gray fluxes first

      chck=float(iflag)

          tg(1:nz) = 0.
      
          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,1),fd(1,1),chck,asym(1,ib))

!           if there is continuum abs. do it now since tg=0
!
     if(nuum(ib).eq.1)then

          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,5),fd(1,5),1.,asym(1,ib))
     endif

!           do the 1st gas


      do ik1=1,npsb(ig1,ib)

          do iz=2,nz
            uu=min(ulim(ig1,ib),u(iz,ig1))
            tg(iz) = xp(ig1,ik1,ib)*uu*  &
               (pl(iz)/prf(ig1,ib))**alpha(ig1,ik1,ib)*  &
               (trf(ig1,ib)/tl(iz))**beta(ig1,ik1,ib)
          enddo

          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

          fact = wght(ig1,ik1,ib)
          do iz=1,nz
            fu(iz,2) = fu(iz,2)+fact*fu(iz,6)
            fd(iz,2) = fd(iz,2)+fact*fd(iz,6)
          enddo
        enddo

!           do the 2nd gas

      do ik2=1,npsb(ig2,ib)

          do iz=2,nz
            uu=min(ulim(ig2,ib),u(iz,ig2))
            tg(iz) = xp(ig2,ik2,ib)*uu*  &
               (pl(iz)/prf(ig2,ib))**alpha(ig2,ik2,ib)*  &
               (trf(ig2,ib)/tl(iz))**beta(ig2,ik2,ib)
          enddo

          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

          fact = wght(ig2,ik2,ib)
          do iz=1,nz
            fu(iz,3) = fu(iz,3)+fact*fu(iz,6)
            fd(iz,3) = fd(iz,3)+fact*fd(iz,6)
          enddo
        enddo

!          do the 3rd gas

      do ik3=1,npsb(ig3,ib)

          do iz=2,nz
            uu=min(ulim(ig3,ib),u(iz,ig3))
            tg(iz) = xp(ig3,ik3,ib)*uu*  &
               (pl(iz)/prf(ig3,ib))**alpha(ig3,ik3,ib)*  &
               (trf(ig3,ib)/tl(iz))**beta(ig3,ik3,ib)
          enddo

!               now do rest of stuff in subroutine

          call flxlw(nz,tg,tp(1,ib),tcr,omgp(1,ib),  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu(1,6),fd(1,6),chck,asym(1,ib))

!               sum the pseudo-band fluxes to get total flux

          fact = wght(ig3,ik3,ib)
          do iz=1,nz
            fu(iz,4) = fu(iz,4)+fact*fu(iz,6)
            fd(iz,4) = fd(iz,4)+fact*fd(iz,6)
          enddo
      enddo

!           sum the fluxes together: FESFT method

!           do iz=1,nz
!             ufac=max(1.e-14,fu(iz,1))
!             dfac=max(1.e-14,fd(iz,1))
!             dx=dfac-ufac
!             fbu=fu(iz,2)/ufac*fu(iz,3)
!    +         /ufac*fu(iz,4)/ufac*fu(iz,1)*(fu(iz,5)
!    +         /ufac*float(nuum(ib))+1.-float(nuum(ib)))
!             dx1=(fd(iz,2)-fu(iz,2))/dx
!             dx2=(fd(iz,3)-fu(iz,3))/dx
!             dx3=(fd(iz,4)-fu(iz,4))/dx
!             dxc=(fd(iz,5)-fu(iz,5))/dx*float(nuum(ib))
!    +         +1.-float(nuum(ib))
!             difflxb=dx1*dx2*dx3*dxc*dx
!             fbd=difflxb+fbu
!             flxlu(iz)=flxlu(iz)+fbu
!             flxld(iz)=flxld(iz)+fbd
!           enddo

! above is the old section, below is the new bounded section:


      do iz=1,nz
        ufac=max(1.e-14,fu(iz,1))
        dfac=max(1.e-14,fd(iz,1))
        tu1 = max(0.0,min(1.1,fu(iz,2)/ufac))
        tu2 = max(0.0,min(1.1,fu(iz,3)/ufac))
        tu3 = max(0.0,min(1.1,fu(iz,4)/ufac))
        fbu = tu1*tu2*tu3*fu(iz,1)

        dx=fd(iz,1)-fu(iz,1)
        adx = abs(dx)
        fdx = sign(1.0,dx)
!        fdx = dx/adx
        dxc = max(1.e-14,adx)*fdx
        tn1c = (fd(iz,2)-fu(iz,2))/dxc + 1.
        tn2c = (fd(iz,3)-fu(iz,3))/dxc + 1.
        tn3c = (fd(iz,4)-fu(iz,4))/dxc + 1.
        tn1 = min(max(tn1c,0.),2.) - 1.
        tn2 = min(max(tn2c,0.),2.) - 1.
        tn3 = min(max(tn3c,0.),2.) - 1.
        difflxb = tn1*tn2*tn3*dx
        fbd=difflxb+fbu

        flxlu(iz)=flxlu(iz)+fbu
        flxld(iz)=flxld(iz)+fbd
      enddo



!-------------------------------------------------------------------------

   else
      write(*,5000)
5000       format('Two is the maximum amount of overlapping',/  &
            'gasses allowed: if you really want to have ',/  &
            'more gasses overlap, they better have few',/  &
            'pseudo-bands for all of them, and in any',/  &
            'case it will cost you. To do that come into',/  &
            'the code and modify it here, look at previous',/  &
            'structures to see how it is done. Its easy',/  &
            'but BEWARE: ITS BOUND TO BE HORRIBLY EXPENSIVE')
      stop
   endif
111     continue

enddo

return
end subroutine lwrad

!===============================================================================

subroutine flxsw(nz,tg,tp,tcr,omgp,  &
              alb_b,alb_d,slr,amu0,t,r,tc,  &
              sigu,sigd,re,vd,td,vu,  &
              fu,fd,asym,i,time,vdflx)

implicit none
integer :: nz,i
real :: tg(nz),tp(nz),tcr(nz),omgp(nz),  &
     alb_b,alb_d,slr,amu0,asym(nz),time

!     local variables

real :: t(nz),r(nz),tc(nz),sigu(nz),sigd(nz)
real :: re(nz),vd(nz),td(nz),vu(nz),expt1
real(kind=8) :: su, sd

!     output variables

real :: fu(nz),fd(nz)
real, intent(out) :: vdflx

!     at this stage we have gotten rid of all the stuff
!     that are specific absorber gas dependent. The effects
!     of all the absorbing gasses in this specific pseudo-
!     band overlap (or no-overlap) is accumulated in tg.

!     parameter(asym=0.85,diffac=1.66,eps=1.e-6,etmp=0.367879)
real, parameter :: diffac=1.66,eps=1.e-6,etmp=0.367879
real, parameter :: e2=0.5,e3=0.3333,e4=0.25,e5=0.2,e6=0.166667
real, parameter :: e7=0.142857

integer :: iz
real :: pi,tau,omg0,af,fact,beta0,g1,g2,gg,rinf,ggtau,expp1,expp2,denom  &
       ,cc,g3,g4,aa,bb,tcm2,exp1,exp2
real :: amu0i

!     get total optical depth, single scattering albedo
!     and assymetry parameter

pi = 3.14159
amu0i = 1.0 / amu0

   re(nz) = 0.
   vd(nz) = 0.
   expt1 = 1.
   fd(nz) = amu0*slr
   tc(nz) = 0.

do iz=nz,2,-1
   
   tau = tg(iz) + tp(iz) + tcr(iz)
   omg0 = min(.999999,  &
        (tcr(iz) + omgp(iz) * tp(iz)) / max(1.e-20, tau)  &
        )
   af = asym(iz) * omgp(iz) * tp(iz) / (omg0 * max(1.e-20, tau))
   
!           do delta-m scaling (wiscombe)

   fact = af * af

   tau = (1.0 - omg0 * fact) * tau
   omg0 = ( (1.0 - fact) * omg0) / (1.0 - omg0 * fact)

!           determine the ODE matrix coefficients (Ritter and Geleyn)

   beta0 = (4.0 + af) / (8.0 * (1.0 + af))
   g1 = diffac * (1.0 - omg0 * (1.0 - beta0))
   g2 = diffac * omg0 * beta0
   gg = sqrt(g1**2 - g2**2)

!           determine the local (true) reflection and transmission coefficients

   rinf = g2 / (gg + g1)
   ggtau = tau * gg
   expp1 = exp(-ggtau)
   expp2 = expp1**2
   denom = (1.0 - rinf**2 * expp2)
   t(iz) = ((1.0 - rinf**2) * expp1) / denom
   r(iz) = rinf * (1.0 - expp2) / denom
  
!        get the source functions, go from top down to accomodate solar terms

   if ((gg - amu0i) < eps .and. gg > amu0i) then
      fact = 1.0 / (gg**2 - 1.0 / (amu0 + eps)**2)
   elseif ( (amu0i - gg) < eps .and. gg <= amu0i) then
      fact = 1.0 / (gg**2 - 1.0 / (amu0 - eps)**2)
   else
      fact = 1.0 / (gg**2 - amu0i**2)
   endif

   cc = omg0 * slr * fact
   g3 = 0.5 - 0.75 * af * amu0 / (1.0 + af)
   g4 = 1.0 - g3
   aa = g3 * (g1 - amu0i) + g4 * g2
   bb = g4 * (g1 + amu0i) + g3 * g2
   tc(iz-1) = tc(iz) + tau
   tcm2 = tc(iz-1) * amu0i
   exp2 = exp(-tcm2)
   exp1 = expt1
   expt1 = exp2
   sigu(iz) = cc * ( (aa - r(iz) * bb) * exp1 - aa * t(iz) * exp2 )
   sigd(iz) = cc * ( -bb * t(iz) * exp1 + (bb - r(iz) * aa) * exp2)
   fd(iz-1) = amu0 * slr * exp2
   td(iz)   = 1.0 - re(iz) * r(iz)
   re(iz-1) = r(iz) + t(iz)**2 * re(iz) / td(iz)

   vd(iz-1) = sigd(iz) + ( t(iz) * vd(iz) + t(iz) * re(iz) * sigu(iz) ) /   &
        td(iz)
   vu(iz-1) = ( r(iz) * vd(iz) + sigu(iz) ) / td(iz)
enddo

!     specify the boundary conditions

fu(1) = (alb_d * vd(1) + alb_b * slr * amu0 * exp(-tc(1) * amu0i)) /   &
     (1.0 - alb_d * re(1))
vdflx = re(1) * fu(1) + vd(1)

!     do adding, going from top down


!     calculate fluxes going up through the layers

do iz = 2, nz
   fd(iz-1) = re(iz-1) * fu(iz-1) + vd(iz-1) + fd(iz-1)
   fu(iz) = t(iz) * fu(iz-1) / td(iz) + vu(iz-1)
enddo

return
end subroutine flxsw

!===============================================================================

subroutine flxlw(nz,tg,tp,tcr,omgp,  &
                  src,t,r,tc,  &
                  sigu,sigd,re,vd,td,vu,  &
                  fu,fd,chck,asym)
use consts_coms, only: pi1
implicit none
integer :: nz
real :: tg(nz),tp(nz),tcr(nz),omgp(nz),src(nz),  &
     chck,asym(nz)

!     local variables

real :: t(nz),r(nz),tc(nz),sigu(nz),sigd(nz)
real :: re(nz),vd(nz),td(nz),vu(nz)

!     output variables

real :: fu(nz),fd(nz)

  real, parameter :: diffac=1.66
  
  integer :: iz
  real :: tau,omg0,af,fact,beta0,g1,g2,gg,rinf,ggtau,expp1,expp2,aa,bb,cc

!     get total optical depth, single scattering albedo
!     and assymetry parameter

do iz=nz,2,-1

!           calculate the single scattering albedo from both gasses
!           and particulates

   tau = tg(iz) + tp(iz) + tcr(iz) * chck

   if ( tau == 0.0 ) then

      omg0 = 0.
      af = 0.

   else
      
      omg0 = min(.999999,  &
           omgp(iz) * tp(iz) / max(1.e-20, tau))
      af = asym(iz)
      
   endif

!           do delta-m scaling (wiscombe)

   fact = af * af
   tau = (1.0 - omg0 * fact) * tau
   omg0 = ((1.0 - fact) * omg0) / (1.0 - omg0 * fact)

!           determine the ODE matrix coefficients (Ritter and Geleyn)

   beta0 = (4.+af)/(8.*(1.+af))
   g1 = diffac*(1.-omg0*(1.-beta0))
   g2 = diffac*omg0*beta0
   gg = sqrt(g1**2-g2**2)
   
!           determine the local (true) reflection and transmission coefficients

   rinf  = g2/(gg+g1)
   ggtau = gg*tau
   expp1=exp(-ggtau)
   expp2=expp1**2

   t(iz) = (1.0-rinf**2)*expp1 / (1.0-rinf**2*expp2)
   r(iz) = rinf*(1.-expp2) /  &
        (1.-rinf**2*expp2)
   
!        get the source functions, go from top down to accomodate solar terms

   if (tau < 4.e-2) then      !changed June 12 after Jerry's recom.
      sigu(iz) = 0.5 * pi1 * (src(iz) + src(iz-1)) * tau * diffac
      sigd(iz) = sigu(iz)
   else
      aa = (g1 + g2) * (1.0 - r(iz)) - (1.0 + r(iz) - t(iz)) / tau
      bb = -(g1 + g2) * t(iz) + (1.0 + r(iz) - t(iz)) / tau
      cc = diffac * pi1 * (1.0 - omg0) / gg**2
      sigu(iz) = cc*(aa*src(iz)+bb*src(iz-1))
      sigd(iz) = cc*(bb*src(iz)+aa*src(iz-1))
      if (sigu(iz) < 0.0 .or. sigd(iz) < 0.0)then
!         print*,'negative source in flxlw: iz',iz
!         print*,'aa,bb,cc:',aa,bb,cc
!         print*,'src(iz), src(iz-1):',src(iz),src(iz-1)
!         print*,'g1, g2, tau:',g1,g2,tau
!         print*,'r(iz), t(iz):',r(iz),t(iz)
!         print*,'sigu(iz), sigd(iz):',sigu(iz),sigd(iz)
!         stop
         sigu(iz) = 0.5 * pi1 * (src(iz) + src(iz-1)) * tau * diffac
         sigd(iz) = sigu(iz)
      endif
   endif
   
enddo

!     do adding
!     initialize

re(nz) = 0.
vd(nz) = 0.
fd(nz) = 0.
fu(1)  = pi1*src(1)


do iz=nz,2,-1
   td(iz)   = 1. - re(iz)*r(iz)
   re(iz-1) = r(iz) + t(iz)**2*re(iz) /  &
        td(iz)
   vd(iz-1) = sigd(iz) + ( t(iz)*vd(iz) +  &
        t(iz)*re(iz)*sigu(iz) ) /  &
        td(iz)
   vu(iz-1) = ( r(iz)*vd(iz) + sigu(iz) ) /  &
        td(iz)
enddo


!     calculate fluxes going up through the layers

do iz=2,nz
   fd(iz-1) = re(iz-1)*fu(iz-1) + vd(iz-1)
   fu(iz)   = t(iz)*fu(iz-1)/td(iz)+vu(iz-1)
enddo

return
end subroutine flxlw

!===============================================================================

subroutine rayleigh(wlnlo,wlnhi,rayavg)
use mem_harr, only: sun
implicit none

real, intent(in)  :: wlnlo,wlnhi
real, intent(out) :: rayavg
real :: an0,ssum,h1,sum,wl1,f1,fac,f2,h2,wl2,al
integer :: num,i

!     this stuff comes from Houghton 1985 (book), see also Slingo
!     and Schrecker 1982 (QJ)
!     calculate constant part (flux weighted) of rayleigh scattering
!     rayleigh scattering = rayavg*delta z*press/temp
!     written JV May 1993

real, parameter :: an01=1.000064328,an02=0.0294981,an03=.2554e-3
real, parameter :: an0d1=146.,an0d2=41.

an0(al) = an01+an02/(an0d1-al**2)+an03/(an0d2-al**2)

num = int(wlnhi-wlnlo)+1
num = min(max(25,num),5000)
sum = 0.
ssum= 0.
wl1 = wlnlo
f1  = (an0(wl1)**2-1.)**2/wl1**4*sun(1.e4/wl1)
h1  = sun(1.e4/wl1)
do i=2,num
   fac = (real(i)-1)/(real(num)-1.)
   wl2 = wlnlo+(wlnhi-wlnlo)*fac
   f2  = (an0(wl2)**2-1.)**2/wl2**4*sun(1.e4/wl2)
   h2  = sun(1.e4/wl2)
   sum = sum + 0.5*(f1+f2)*(wl2-wl1)
   ssum= ssum + 0.5*(h1+h2)*(wl2-wl1)
   wl1 = wl2
   f1  = f2
   h1  = h2
enddo
rayavg = 0.945319e-2*sum/ssum

return
end subroutine rayleigh

!===============================================================================

subroutine csband(wlnlo,wlnhi,csavg)
implicit none
real :: wlnlo,wlnhi,csavg
real :: cs,wnhi,wnlo,sum,wn1,f1,fac,f2,wn,wn2,plkavg
integer :: num,i
!
!     calculate the blackbody flux (t=296.) weighted self-
!     broadening coefficient for water vapor. See Kniezys
!     et al 1980 (LOWTRAN 5). Units are converted to have
!     the water vapor content input a Pascal.

cs(wn)  = 4.18+5578.*exp(-7.87e-3*wn)

wnhi = 1.e4/wlnlo
wnlo = 1.e4/wlnhi
num  = int(wnhi-wnlo)+1
num  = min(max(25,num),5000)
sum  = 0.
wn1  = wnlo
f1   = cs(wn1)
do i = 2,num
   fac = (real(i)-1)/(real(num)-1.)
   wn2 = wnlo+(wnhi-wnlo)*fac
   f2  = cs(wn2)
   sum = sum + 0.5*(f1+f2)*plkavg(wn1,wn2,296.)
   wn1 = wn2
   f1  = f2
enddo
csavg = sum/(plkavg(wnlo,wnhi,296.)*1013250.*9.81)
return
end subroutine csband

!===============================================================================

real FUNCTION  PLKAVG( WNUMLO, WNUMHI, T )

implicit none

!
!      COMPUTES PLANCK FUNCTION INTEGRATED BETWEEN TWO WAVENUMBERS
!
!      I N P U T
!
!         WNUMLO : LOWER WAVENUMBER ( INV CM ) OF SPECTRAL INTERVAL
!         WNUMHI : UPPER WAVENUMBER
!         T       : TEMPERATURE (K)
!
!      O U T P U T
!
!         PLKAVG : INTEGRATED PLANCK FUNCTION ( WATTS/SQ M )
!
!      REFERENCES-- (1) HOUGHTON,PHYSICS OF ATMOSPHERES,APPENDIX 7
!                   (2) SPECIFICATIONS OF THE PHYSICAL WORLD: NEW VALUE
!                       OF THE FUNDAMENTAL CONSTANTS, DIMENSIONS/N.B.S.
!                       JAN. 1974
!
!      METHOD-- HOUGHTON'S EXPONENTIAL SERIES IS USED FOR V.GT.VCUT
!               ( 'V' IS HOUGHTON'S NOTATION ) AND HIS POWER SERIES
!               IN V FOR V.LE.VCUT.  MORE TERMS ARE TAKEN IN THE
!               EXPONENTIAL SERIES, THE LARGER V IS.  ( NOTE THAT
!               HOUGHTON'S ASSESSMENT THAT THE POWER SERIES IS USEFUL
!               FOR  V.LT.2*PI  IS INCORRECT--VCUT MUST BE LESS THAN
!               2 JUST IN ORDER TO GET 4-5 SIGNIFICANT DIGITS. )
!
!      ACCURACY-- 6 SIGNIFICANT DIGITS
!
!          *** ARGUMENTS
REAL     T, WNUMLO, WNUMHI
!          *** LOCAL VARIABLES
!                 C2       :  SECOND RADIATION CONSTANT
!                 CONA     :  15 / PI**4 / 13305600
!                 CONB     :  440 / 9
!                 D        :  INTEGRAL OF NORMALIZED PLANCK FUNCTION
!                             FROM 0 TO CURRENT WAVENUMBER
!                 EX       :  EXP( - V )
!                 F15PI4   :  15 / PI**4
!                 MMAX     :  NO. OF TERMS TO TAKE IN EXPONENTIAL SERIES
!                 MV       :  MULTIPLES OF *V*
!                 SIGDPI   :  STEFAN-BOLTZMANN CONSTANT DIVIDED BY PI
!                 V        :  H*C*NU / (K*T), WHERE H=PLANCKS CONSTANT,
!                             C=SPEED OF LIGHT, K=BOLTZMANN CONSTANT,
!                             NU=WAVENUMBER
!                 VCUT     :  POWER-SERIES CUTOFF POINT
!                 VCP      :  EXPONENTIAL SERIES CUTOFF POINTS
REAL     C2, CONA, CONB, D( 2 ), EX, F15PI4, MV, SIGDPI,  &
         V, VCUT, VCP( 7 ), VSQ
DATA  C2     / 1.438786    /,  CONA    / 0.11573303E-7 /,  &
      CONB   / 48.888889   /,  F15PI4 / 0.15398973    /,  &
      SIGDPI / 1.804919E-8 /,  VCUT    / 1.5 /,  &
      VCP    / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      
integer :: i,mmax,m
!
!
IF( T.LE.0.0 .OR. WNUMHI.LE.WNUMLO .OR. WNUMLO.LE.0. ) stop 'ahah'
!
DO  50  I = 1, 2
   IF( I.EQ.1 )  V = ( C2 / T ) * WNUMLO
   IF( I.EQ.2 )  V = ( C2 / T ) * WNUMHI
   IF( V.LT.VCUT )  THEN
!
!            *** HOUGHTON'S POWER SERIES (FACTORED USING HORNER'S RULE)
      VSQ = V**2
      D( I ) = 1.0 - CONA *VSQ * V * ( 4435200. + V * ( -1663200.  &
               + V * ( 221760. + VSQ * ( -2640. + VSQ * ( CONB  &
               - VSQ ) ) ) ) )
   ELSE
!
!             *** HOUGHTON'S EXPONENTIAL SERIES
      MMAX = 0
!               *** SET UPPER LIMIT OF SERIES DEPENDING ON VALUE OF V
20       MMAX = MMAX + 1
         IF ( V.LT.VCP( MMAX ) )  GO TO 20
!
      EX = EXP( -V )
      D( I ) = EX * ( 6. + V * ( 6. + V * ( 3. + V ) ) )
!
      DO  30  M = 2, MMAX
         MV = M * V
         D( I ) = D( I ) + EX**M * ( 6. + MV * ( 6. + MV *  &
                           ( 3. + MV ) ) ) / M**4
30       CONTINUE
!
      D( I ) = F15PI4 * D( I )
   END IF
!
50 CONTINUE
!
PLKAVG = SIGDPI * T**4 * ( D( 1 ) - D( 2 ) )
!
RETURN
END FUNCTION PLKAVG

!===============================================================================

real function sunavg(wnumlo, wnumhi, solcon)
use mem_harr, only: sun
implicit none

real :: wnumlo, wnumhi, solcon
real :: scln=1372.844

integer :: num,i
real :: v1,s1,v2,s2,fac2,sum
! Spectral solar flux  (after LOWTRAN7)
!
! On input:
! wnumlo --- lower wavelength of a spectral interval (wavenumbers)
! wnumhi --- upper wavelength of a spectral interval (wavenumbers)
! solcon - value of solar "constant" in Watts/ sq m
!
! On output:
! Total extraterrestrial solar flux (W/sq m) between
! wavelengths wnumlo, wnumhi
! Trapezoidal integration is used with the resolution 1 cm^-1
! or such that there is at least 25 points between wnumhi and
! wnumlo but no more than 5000 points.
!
num = int(wnumhi-wnumlo)+1
num = min(max(25, num),5000)
sum = 0.
v1 = wnumlo
s1 = sun(v1)
do i=2, num
  fac2 = (real(i)-1.)/(real(num)-1.)
  v2 =  wnumlo+(wnumhi-wnumlo)*fac2
  s2 = sun(v2)
  sum = sum + (s1+s2)/2.*(10000./v1-10000./v2)
  v1 = v2
  s1 = s2
enddo
sunavg = (solcon/scln)*sum

return
end function sunavg
