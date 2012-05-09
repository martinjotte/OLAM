subroutine diagvel_t3d(mrl)

  use mem_basic, only: vc, uc, wc, vxe, vye, vze
  use misc_coms, only: meshtype
  implicit none

  integer, intent(in) :: mrl

  if (meshtype == 1) then

     call vel_t3d_tri(mrl, uc, wc, vxe, vye, vze)

  elseif (meshtype == 2) then

     call vel_t3d_hex(mrl, vc, wc, vxe, vye, vze)

  endif

  return
end subroutine diagvel_t3d


!===============================================================================


subroutine vel_t3d_hex(mrl,vs,ws,vxe,vye,vze)

use mem_ijtabs, only: jtab_w, itab_v, itab_w
use mem_grid,   only: mza, mva, mwa, lpw, lpv, vnx, vny, vnz, wnx, wny, wnz
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer, intent(in) :: mrl

real, intent(in)  :: vs(mza,mva)
real, intent(in)  :: ws(mza,mwa)

real, intent(out) :: vxe(mza,mwa)
real, intent(out) :: vye(mza,mwa)
real, intent(out) :: vze(mza,mwa)

integer :: j,iw,npoly,kb,k,jv,iv
real    :: farv2, wst

! Horizontal loop over W columns

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,npoly,kb,k,jv,iv,wst,farv2) 
do j = 1,jtab_w(16)%jend(mrl); iw = jtab_w(16)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   npoly = itab_w(iw)%npoly
   kb = lpw(iw)

! Vertical loop over T levels

   do k = kb, mza-1

! Diagnose 3D earth-velocity vector at T points; W contribution first

      wst = 0.5 * (ws(k-1,iw) + ws(k,iw))

      vxe(k,iw) = wst * wnx(iw)
      vye(k,iw) = wst * wny(iw)
      vze(k,iw) = wst * wnz(iw)

   enddo

! Loop over V neighbors of this W cell

   do jv = 1, npoly

      iv = itab_w(iw)%iv(jv)

      farv2 = 2. * itab_w(iw)%farv(jv)

! Vertical loop over T levels

      do k = kb, mza-1

! Diagnose 3D earth-velocity vector at T points; VC contribution

         vxe(k,iw) = vxe(k,iw) + farv2 * vs(k,iv) * vnx(iv)
         vye(k,iw) = vye(k,iw) + farv2 * vs(k,iv) * vny(iv)
         vze(k,iw) = vze(k,iw) + farv2 * vs(k,iv) * vnz(iv)

      enddo
      
   enddo
   
   vxe(2:kb-1,iw) = vxe(kb,iw)
   vye(2:kb-1,iw) = vye(kb,iw)
   vze(2:kb-1,iw) = vze(kb,iw)

   vxe(mza,iw) = vxe(mza-1,iw)
   vye(mza,iw) = vye(mza-1,iw)
   vze(mza,iw) = vze(mza-1,iw)

enddo
!$omp end parallel do
call rsub('Wa',16)

return
end subroutine vel_t3d_hex


!===============================================================================


subroutine vel_t3d_tri(mrl,us,ws,vxe,vye,vze)

use mem_ijtabs, only: jtab_w, itab_v, itab_w
use mem_grid,   only: mza, mua, mwa, lpw, lpv, vnx, vny, vnz, wnx, wny, wnz
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer, intent(in) :: mrl

real, intent(in)  :: us(mza,mua)
real, intent(in)  :: ws(mza,mwa)

real, intent(out) :: vxe(mza,mwa)
real, intent(out) :: vye(mza,mwa)
real, intent(out) :: vze(mza,mwa)

return
end subroutine vel_t3d_tri




