subroutine check_nans(icall)

  use mem_basic,  only: sh_w,rho,thil,sh_v,wc,wmc,press,vmc,vc
  use mem_micro,  only: sh_c,sh_r
  use mem_grid,   only: mza,mwa,lpw,volti,mva,lpv
  use mem_tend,   only: thilt
  use misc_coms,  only: io6, iparallel
  use mem_ijtabs, only: itab_v, itab_w
  use mem_para,   only: myrank

#ifdef IEEE_ARITHMETIC
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
#endif

  implicit none

  integer :: i,k,iv
  integer, intent(in) :: icall

#ifdef IEEE_ARITHMETIC

  do i = 2,mwa
     do k = lpw(i),mza
        if ( ieee_is_nan(sh_w (k,i)) .or.  &
             ieee_is_nan(rho  (k,i)) .or.  &
             ieee_is_nan(thil (k,i)) .or.  &
             ieee_is_nan(press(k,i)) .or.  &
             ieee_is_nan(wc   (k,i)) .or.  &
             ieee_is_nan(wmc  (k,i)) .or.  &
             ieee_is_nan(sh_v (k,i)) ) then
           
           write(*,*) 'Node ', myrank
           write(*,*) 'NaN at ', k, i, icall
           write(io6,*) ''
           write(io6,*) 'check_nans',k,i,icall
           write(io6,*) 'sh_w,rho',sh_w(k,i),rho(k,i)
           write(io6,*) 'thil,press',thil(k,i),press(k,i)
           write(io6,*) 'wc, wmc, sh_v',wc(k,i),wmc(k,i),sh_v(k,i)
           stop
        endif

     enddo
  enddo

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

#endif
end subroutine check_nans



subroutine check_pos(icall)

  use mem_ijtabs, only: jtab_w, jtw_prog, itab_w, istp, mrl_endl
  use mem_grid,   only: lpw, mza
  use var_tables, only: num_scalar, scalar_tab
  use oname_coms, only: nl
  use mem_para

  implicit none

  integer, intent(in) :: icall
  
  integer :: j, iw, k, mrl, n

  if (nl%iscal_monot < 1) return

  ! Return if this is not the end of the long timestep on any MRL
  mrl = mrl_endl(istp)
  if (mrl == 0) return

  do n = 2, num_scalar

     if (scalar_tab(n)%name == 'Q2') cycle
     if (scalar_tab(n)%name == 'Q6') cycle
     if (scalar_tab(n)%name == 'Q7') cycle

     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
        do k = lpw(iw), mza
           if ( scalar_tab(n)%var_p(k,iw) < 0.0 .or. &
                scalar_tab(n)%var_p(k,iw) /= scalar_tab(n)%var_p(k,iw) ) then
              write(*,*) trim(scalar_tab(n)%name), " equals ",  &
                   scalar_tab(n)%var_p(k,iw), " at ", k, iw, icall, myrank
           endif
        enddo
     enddo
  enddo

end subroutine check_pos
