subroutine check_nans(icall)

  use mem_basic,  only: sh_w,rho,thil,sh_v,wc,wmc,press,umc,uc,vmc,vc
  use mem_micro,  only: sh_c,sh_r
  use mem_grid,   only: mza,mwa,lpw,volti,mua,mva,lpu,lpv
  use mem_tend,   only: thilt, vmt
  use micro_coms, only: level
  use misc_coms,  only: io6, iparallel, meshtype
  use mem_ijtabs, only: itab_v, itab_w
  use mem_para,   only: myrank
  use var_tables, only: num_scalar, scalar_tab

#ifdef IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif

  implicit none

  integer :: i,k,iv
  integer, intent(in) :: icall

#ifdef IEEE_ARITHMETIC

  do i = 2,mwa
     do k = lpw(i),mza
        if (ieee_is_nan(sh_w (k,i)) .or.  &
             ieee_is_nan(rho  (k,i)) .or.  &
             ieee_is_nan(thil (k,i)) .or.  &
             ieee_is_nan(thilt(k,i)) .or.  &
             ieee_is_nan(press(k,i)) .or.  &
             ieee_is_nan(wc   (k,i)) .or.  &
             ieee_is_nan(wmc  (k,i)) .or.  &
             ieee_is_nan(sh_v (k,i)) .or.  &
             ieee_is_nan(scalar_tab(1)%var_t(k,i)) .or. &
             ieee_is_nan(scalar_tab(2)%var_t(k,i)) .or. &
             thil(k,i) < 100.0) then
           
           write(*,*) 'Node ', myrank
           write(*,*) 'NaN at ', k, i, icall
           write(io6,*) ''
           write(io6,*) 'check_nans',k,i,icall
           write(io6,*) 'sh_w,rho,thil',sh_w(k,i),rho(k,i),thil(k,i)
           write(io6,*) 'thilt,press',thilt(k,i),press(k,i)
           write(io6,*) 'wc, wmc, sh_v',wc(k,i),wmc(k,i),sh_v(k,i)
           stop
        endif

     enddo
  enddo

  if (meshtype == 1) then

     do i = 2,mua
        do k = lpu(i), mza
           if (ieee_is_nan(umc(k,i)) .or.  &
                ieee_is_nan(uc (k,i))) then
              write(io6,*) ''
              write(*,*) 'Node ', myrank
              write(*,*) 'NaN at ', k, i, icall
              write(io6,*) 'check_nans',k,i,icall
              write(io6,*) 'umc, uc: ', umc(k,i), uc(k,i)
              stop
           endif
        enddo
     enddo
     
  else
     
     do i = 2,mva
        do k = lpv(i), mza
           if (ieee_is_nan(vmc(k,i)) .or.  &
                ieee_is_nan(vc (k,i))) then
              write(io6,*) ''
              write(*,*) 'Node ', myrank
              write(*,*) 'NaN at ', k, i, icall
              write(io6,*) 'check_nans',k,i,icall
              write(io6,*) 'vmc, vc: ', vmc(k,i), vc(k,i)
              stop
           endif
        enddo
     enddo

  endif

#endif
end subroutine check_nans

