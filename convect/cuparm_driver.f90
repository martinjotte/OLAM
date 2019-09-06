!===============================================================================
! OLAM was originally developed at Duke University by Robert Walko, Martin Otte,
! and David Medvigy in the project group headed by Roni Avissar.  Development
! has continued by the same team working at other institutions (University of
! Miami (rwalko@rsmas.miami.edu), the Environmental Protection Agency, and
! Princeton University), with significant contributions from other people.

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

!===============================================================================
subroutine cuparm_driver()

  use mem_grid,         only: mwa, mza, lpw, volt, arw0
  use module_cu_g3,     only: grell_driver
  use module_cu_gf,     only: gf_driver
  use module_cu_kfeta,  only: cuparm_kfeta, kf_lutab
  use module_cu_tiedtke,only: cuparm_tiedtke
  use module_cu_emanuel,only: cuparm_emanuel
  use misc_coms,        only: io6, time_istp8, time_istp8p, nqparm, confrq, &
                              dtlong, mstp, runtype, dtlm, iparallel
  use mem_ijtabs,       only: itab_w, jtab_w, mrl_begl, istp, mrls, jtw_prog, jtw_wadj
  use mem_cuparm,       only: thsrc, rtsrc, aconpr, conprr, &
                              kcutop, kcubot, qwcon, iactcu, cbmf, cddf, kddtop, &
                              cu_pcpflx, kudbot, kddmax
  use mem_tend,         only: thilt, rr_wt
  use mem_basic,        only: rho, rr_w, theta, tair, thil
  use olam_mpi_atm,     only: mpi_send_w, mpi_recv_w
  use oname_coms,       only: nl
  use var_tables,       only: num_cumix
  use consts_coms,      only: alvlocp, r8
  use olam_mpi_atm,     only: mpi_recv_w, mpi_send_w

  implicit none

  integer, save :: init_kf = 0
  integer       :: j, iw, k, mrl, mrlw
  real          :: dthmax, dtlong4, confrq4, confrq4i
  integer       :: iwqmax, kqmax, km, km1
  integer       :: ka, kb
  real(r8)      :: qsum, tsum, vtot
  real          :: rt2

! For KF_eta parameterization, initialize scheme if needed

  if ( any(nqparm(1:mrls) == 3) ) then
     if (init_kf == 0) then
        init_kf = 1
        call kf_lutab()
     endif
  endif

  dtlong4  = real(dtlong)

! Check whether it is time to update cumulus parameterization tendencies

  if ((istp == 1 .and. mod(time_istp8p, confrq) < dtlong) .or. &
      (istp == 1 .and. mstp == 0 .and. runtype == 'HISTADDGRID')) then

! Print message that cumulus parameterization is being computed

     write(io6, '(A,f0.2,A)') &
           ' Convective tendencies updated at ', time_istp8/3600., &
           ' hrs into simulation'

     confrq4  = real(confrq)
     confrq4i = 1. / confrq4

     !$omp parallel do
     do iw = 1, mwa
        thsrc(:,iw) = 0.0
        rtsrc(:,iw) = 0.0
        qwcon(:,iw) = 0.0

        cu_pcpflx(:,iw) = 0.0

        conprr(iw) =  0.0
        kcutop(iw) = -1
        kddtop(iw) = -1
        kcubot(iw) = -1
        kudbot(iw) = -1
        kddmax(iw) = -1
        iactcu(iw) =  0
        cddf  (iw) = 0.0
     enddo
     !$omp end parallel do

! Loop over all IW grid cells where cumulus parameterization may be done

!----------------------------------------------------------------------
     !$omp parallel do private(iw,mrlw,km,km1,k,ka,kb,qsum,tsum,vtot) &
     !$omp             schedule( guided )
     do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j) ! jend(1) for mrl = 1
!----------------------------------------------------------------------

! MRL for current IW column

        mrlw = itab_w(iw)%mrlw

! Select cumulus convection scheme based on MRL of current IW column

        if (nqparm(mrlw) == 1) then

! Tiedtke convective parameterization

           km = mza + 1 - lpw(iw) ! km  = # of T levels in cuparm_tiedtke
           km1 = km + 1           ! km1 = # of W levels in cuparm_tiedtke

           call cuparm_tiedtke(iw,km,km1,dtlong4,confrq4,confrq4i)

        elseif (nqparm(mrlw) == 2) then

! Grell deep and shallow convection

           call grell_driver( iw, dtlong4 )

        elseif (nqparm(mrlw) == 3) then

! Kain-Fritsch deep convection

           call cuparm_kfeta(iw,dtlong4)

        elseif (nqparm(mrlw) == 4) then

! Emanuel convective parameterization

           call cuparm_emanuel(iw,dtlong4)

        elseif (nqparm(mrlw) == 5) then

! Grell-Freitas deep and shallow convection

           call gf_driver( iw, dtlong4 )

        endif

! Specify a minimum value of convective cloud water

        if (iactcu(iw) == 1) then
           do k = kcubot(iw), kcutop(iw)
              qwcon(k,iw) = max(qwcon(k,iw),1.e-5)
           enddo
        else
           cbmf(iw) = 0.0
           cddf(iw) = 0.0
        endif

! Make sure convective scheme conserves heat and water
! Add heat/moisture within cloud to balance precipitation

        if (iactcu(iw) == 1) then
           ka = lpw(iw)
           kb = min(kcutop(iw) + 1, mza)

           qsum = sum( rtsrc(ka:kb,iw) * volt(ka:kb,iw) ) + conprr(iw) * arw0(iw)
           tsum = sum( thsrc(ka:kb,iw) * volt(ka:kb,iw) ) - conprr(iw) * arw0(iw) * alvlocp
           vtot = sum( volt(kcubot(iw):kcutop(iw),iw) )

           qsum = qsum / vtot
           tsum = tsum / vtot

           do k = kcubot(iw), kcutop(iw)
              rtsrc(k,iw) = rtsrc(k,iw) - qsum
              thsrc(k,iw) = thsrc(k,iw) - tsum
           enddo
        endif

     enddo
     !$omp end parallel do

! Print maximum heating rate in current parallel sub-domain

     dthmax = 0.
     iwqmax = 0
     kqmax  = 0

     ! this will need to be changed to work with OpenMP
     do j = 1,jtab_w(jtw_prog)%jend(1); iw = jtab_w(jtw_prog)%iw(j)
        if (iactcu(iw) == 1) then
           do k = lpw(iw), mza
              if (thsrc(k,iw) > dthmax) then
                 dthmax = thsrc(k,iw)
                 iwqmax = iw
                 kqmax  = k
              endif
           enddo
        endif
     enddo

     write(io6, '(A,I0,A,I0,A,F0.3,A,F0.4)') " MAX CONVECTIVE HEATING RATE AT IW=",  &
           iwqmax, " K=", kqmax, " IS ", dthmax*86400., " K/DAY"

     if (iparallel == 1) then
        call mpi_send_w(1, i1dvara1=kcutop, i1dvara2=kcubot, i1dvara3=iactcu, &
                           r1dvara1=cbmf, r1dvara2=conprr)
     endif

     ! Update cloud fraction if convection was updated
     call calc_3d_cloud_fraction(1)

  endif

! Add current value of convective tendencies to thilt and rr_wt arrays
! every long timestep (whether cumulus parameterization is updated
! this timestep or not)

!----------------------------------------------------------------------
  mrl = mrl_begl(istp)
  if (mrl > 0) then

     !$omp parallel do private(iw,dtlong4,k,rt2) schedule(guided)
     do j = 1,jtab_w(jtw_prog)%jend(mrl); iw = jtab_w(jtw_prog)%iw(j)
!----------------------------------------------------------------------

        ! Skip if no convection in this column

        if (iactcu(iw) /= 1) cycle

        ! Skip if not doing convection on this mrl

        if (nqparm(itab_w(iw)%mrlw) == 0) then
           iactcu(iw) = 0
           conprr(iw) = 0.0
           cycle
        endif

        dtlong4 = dtlm(itab_w(iw)%mrlw)

        aconpr(iw) = aconpr(iw) + conprr(iw) * dtlong4

        do k = lpw(iw), min(kcutop(iw) + 1, mza)
           if (tair(k,iw) > 253.) then
              thilt(k,iw) = thilt(k,iw) + thsrc(k,iw) * theta(k,iw) / tair(k,iw)
           else
              thilt(k,iw) = thilt(k,iw) + thsrc(k,iw) * thil(k,iw) / tair(k,iw)
           endif

           rt2 = max( rtsrc(k,iw), &
                -0.95 * max(rr_w(k,iw), 0.0) * real(rho(k,iw)) / dtlong4 )
           rr_wt(k,iw) = rr_wt(k,iw) + rt2
        enddo

        ! Compute tracer mixing due to convection

        if ( (nl%conv_tracer_mix > 0 .and. num_cumix > 0) &
           .or. (nl%conv_uv_mix > 0) ) then
           call acmcld_mix( iw, dtlong4 )
        endif

     enddo
     !$omp end parallel do

  endif

  if (iparallel == 1) then
     if ((istp == 1 .and. mod(time_istp8p, confrq) < dtlong) .or. &
         (istp == 1 .and. mstp == 0 .and. runtype == 'HISTADDGRID')) then

        call mpi_recv_w(1, i1dvara1=kcutop, i1dvara2=kcubot, i1dvara3=iactcu, &
                           r1dvara1=cbmf, r1dvara2=conprr)

     endif
  endif

end subroutine cuparm_driver



subroutine reset_cuparm()

  use mem_ijtabs, only: jtab_w, jtw_wadj, itab_w
  use mem_cuparm, only: conprr, kcutop, kcubot, iactcu, cbmf, cddf, kudbot, &
                        kddtop, kddmax
  use misc_coms,  only: nqparm

  implicit none

  integer :: j, iw, mrlw

  !$omp parallel do private(iw,mrlw)
  do j = 1,jtab_w(jtw_wadj)%jend(1); iw = jtab_w(jtw_wadj)%iw(j)

     ! MRL for current IW column
     mrlw = itab_w(iw)%mrlw

     if (nqparm(mrlw) == 0) then

        iactcu(iw) =  0

        if (allocated( conprr )) conprr(iw) =  0.0
        if (allocated( kcutop )) kcutop(iw) = -1
        if (allocated( kcubot )) kcubot(iw) = -1
        if (allocated( kudbot )) kudbot(iw) = -1
        if (allocated( kddtop )) kddtop(iw) = -1
        if (allocated( kddmax )) kddmax(iw) = -1
        if (allocated( cbmf   ))   cbmf(iw) =  0.0
        if (allocated( cddf   ))   cddf(iw) =  0.0

     endif

  enddo
  !$omp end parallel do

end subroutine reset_cuparm
