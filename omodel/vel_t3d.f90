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


subroutine vel_t3d_hex(mrl, vs, ws, vxe, vye, vze)

use mem_ijtabs, only: jtab_w, itab_v, itab_w, jtw_prog
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

if (mrl == 0) return

! Horizontal loop over W columns

call psub()
!----------------------------------------------------------------------
!$omp parallel do private(iw,npoly,kb,k,jv,iv,wst) 
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------
call qsub('W',iw)

   npoly = itab_w(iw)%npoly
   kb = lpw(iw)

! Vertical loop over T levels

   do k = kb, mza-1

! Diagnose 3D earth-velocity vector at T points; W contribution first

      wst = itab_w(iw)%ecvec_w * (ws(k-1,iw) + ws(k,iw))

      vxe(k,iw) = wst * wnx(iw)
      vye(k,iw) = wst * wny(iw)
      vze(k,iw) = wst * wnz(iw)

   enddo

! Loop over V neighbors of this W cell

   do jv = 1, npoly

      iv = itab_w(iw)%iv(jv)

! Vertical loop over T levels

      do k = kb, mza-1

! Diagnose 3D earth-velocity vector at T points; VC contribution

         vxe(k,iw) = vxe(k,iw) + itab_w(iw)%ecvec_v(jv) * vs(k,iv) * vnx(iv)
         vye(k,iw) = vye(k,iw) + itab_w(iw)%ecvec_v(jv) * vs(k,iv) * vny(iv)
         vze(k,iw) = vze(k,iw) + itab_w(iw)%ecvec_v(jv) * vs(k,iv) * vnz(iv)

      enddo
      
   enddo
   
   vxe(1:kb-1,iw) = vxe(kb,iw)
   vye(1:kb-1,iw) = vye(kb,iw)
   vze(1:kb-1,iw) = vze(kb,iw)

   vxe(mza,iw) = vxe(mza-1,iw)
   vye(mza,iw) = vye(mza-1,iw)
   vze(mza,iw) = vze(mza-1,iw)

enddo
!$omp end parallel do
call rsub('Wa',16)

return
end subroutine vel_t3d_hex


!===============================================================================


subroutine vel_t3d_tri(mrl, us, ws, vxe, vye, vze)

use mem_ijtabs, only: jtab_w, itab_u, itab_w, jtw_prog
use mem_grid,   only: mza, mua, mwa, lpw
use misc_coms,  only: io6

!$ use omp_lib

implicit none

integer, intent(in) :: mrl

real, intent(in)  :: us(mza,mua)
real, intent(in)  :: ws(mza,mwa)

real, intent(out) :: vxe(mza,mwa)
real, intent(out) :: vye(mza,mwa)
real, intent(out) :: vze(mza,mwa)

integer :: iu1, iu2, iu3
integer :: j, iw, k, kb
real    :: wm
real    :: vxu1, vxu2, vxu3, vxw
real    :: vyu1, vyu2, vyu3, vyw
real    :: vzu1, vzu2, vzu3, vzw

if (mrl == 0) return

! Loop over all primary W points

!$omp parallel do private (iw,iu1,iu2,iu3,k,wm,vxu1,vxu2,vxu3,vxw,  &
!$omp                      vyu1,vyu2,vyu3,vyw,vzu1,vzu2,vzu3,vzw,kb )
do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)

   iu1 = itab_w(iw)%iu(1)
   iu2 = itab_w(iw)%iu(2)
   iu3 = itab_w(iw)%iu(3)

   vxu1 = itab_w(iw)%vxu(1)
   vxu2 = itab_w(iw)%vxu(2)
   vxu3 = itab_w(iw)%vxu(3)
   vxw  = itab_w(iw)%vxw 

   vyu1 = itab_w(iw)%vyu(1)
   vyu2 = itab_w(iw)%vyu(2)
   vyu3 = itab_w(iw)%vyu(3)
   vyw  = itab_w(iw)%vyw 

   vzu1 = itab_w(iw)%vzu(1)
   vzu2 = itab_w(iw)%vzu(2)
   vzu3 = itab_w(iw)%vzu(3)
   vzw  = itab_w(iw)%vzw 

   kb = lpw(iw)

! Vertical loop over T levels

   do k = kb, mza-1

      wm = ws(k-1,iw) + ws(k,iw)

      vxe(k,iw) = vxu1*us(k,iu1) + vxu2*us(k,iu2) + vxu3*us(k,iu3) + vxw*wm
      
      vye(k,iw) = vyu1*us(k,iu1) + vyu2*us(k,iu2) + vyu3*us(k,iu3) + vyw*wm
      
      vze(k,iw) = vzu1*us(k,iu1) + vzu2*us(k,iu2) + vzu3*us(k,iu3) + vzw*wm

   enddo

   vxe(1:kb-1,iw) = vxe(kb,iw)
   vye(1:kb-1,iw) = vye(kb,iw)
   vze(1:kb-1,iw) = vze(kb,iw)

   vxe(mza,iw) = vxe(mza-1,iw)
   vye(mza,iw) = vye(mza-1,iw)
   vze(mza,iw) = vze(mza-1,iw)

enddo
!$omp end parallel do

return
end subroutine vel_t3d_tri




