subroutine lpsource()

  use mem_lp,      only: nlpsrc, lpsrc, nlp, atp, dtlp, relcnt
  use consts_coms, only: pio180, pi2
  use misc_coms,   only: time8, io6
  use map_proj,    only: ps_ec

  implicit none (external, type)

  integer, save :: idum_seed = 1

  real :: sc1, sc3, sc4, denslp, duration, xdisp, ydisp
  real :: radius, angle, rpos, azim

  integer :: ilpsrc, nrelcnt, irelcnt, nlp0

  real, external :: ramran
  external       :: lpe_kiw

  nlp0 = nlp

  ! Loop over all sources

  do ilpsrc = 1,nlpsrc

     ! Skip source if current time (time8) is outside source release period

     if (real(time8) < lpsrc(ilpsrc)%begtime .or. real(time8) > lpsrc(ilpsrc)%endtime) cycle

     sc1 = lpsrc(ilpsrc)%szpwr
     sc3 = lpsrc(ilpsrc)%szmin**(1./sc1)
     sc4 = lpsrc(ilpsrc)%szmax**(1./sc1)
     denslp = 2.0

     ! account for portion of a particle left over from last timestep
     ! and portion of a timestep until completion of release

     if (lpsrc(ilpsrc)%endtime == lpsrc(ilpsrc)%begtime) then ! instantaneous release
        relcnt(ilpsrc) = lpsrc(ilpsrc)%relperdt
     else
        duration = min(real(time8) + dtlp,lpsrc(ilpsrc)%endtime)  &
                 - max(real(time8)       ,lpsrc(ilpsrc)%begtime)
        relcnt(ilpsrc) = relcnt(ilpsrc) + lpsrc(ilpsrc)%relperdt * duration / abs(dtlp)
     endif
     nrelcnt = int(relcnt(ilpsrc))
     relcnt(ilpsrc) = relcnt(ilpsrc) - real(nrelcnt)

     do irelcnt = 1,nrelcnt
        nlp = nlp + 1

        ! random particle placement within source

        if (lpsrc(ilpsrc)%outline == 1) then  ! RECTANGLULAR SOURCE

           xdisp = lpsrc(ilpsrc)%xsize * (ramran(idum_seed) - 0.5)
           ydisp = lpsrc(ilpsrc)%ysize * (ramran(idum_seed) - 0.5)

        elseif (lpsrc(ilpsrc)%outline == 2) then  ! ELLIPTICAL SOURCE

           rpos = sqrt(ramran(idum_seed))
           azim = pi2 * ramran(idum_seed)
           xdisp = rpos * cos(azim) * lpsrc(ilpsrc)%xsize * .5
           ydisp = rpos * sin(azim) * lpsrc(ilpsrc)%ysize * .5

        endif

        radius = sqrt(xdisp**2 + ydisp**2)
        if (lpsrc(ilpsrc)%xsize == 0. .and. lpsrc(ilpsrc)%ysize == 0.) then
           angle = 0.
        else
           angle = atan2(ydisp,xdisp)
        endif
        angle = angle - lpsrc(ilpsrc)%rotate * pio180
        xdisp = radius * cos(angle)
        ydisp = radius * sin(angle)

        atp(nlp)%zp = lpsrc(ilpsrc)%zbot + (lpsrc(ilpsrc)%ztop - lpsrc(ilpsrc)%zbot) * ramran(idum_seed)

        ! Get earth-coordinates for (xdisp,ydisp) point at sea level

        call ps_ec(atp(nlp)%xep0, atp(nlp)%yep0, atp(nlp)%zep0, lpsrc(ilpsrc)%centlat, lpsrc(ilpsrc)%centlon, xdisp, ydisp)

        ! Get particle k,iw from its xep0, yep0, and zep0, and zp

        call lpe_kiw(nlp, lpsrc(ilpsrc)%iwcent)

        ! assign other particle attributes

        atp(nlp)%nsource  = ilpsrc
!e        atp(nlp)%nspecies = species(ilpsrc)
        atp(nlp)%reltime = real(time8)
        if (lpsrc(ilpsrc)%ifall == 1) then
           atp(nlp)%fallvel = -4.e7 * denslp * (sc3 + (sc4 - sc3)  &
                            * ramran(idum_seed))**(2. * sc1)
        else
           atp(nlp)%fallvel = 0.
        endif
!e        atp(nlp)%upp  = 0.
!e        atp(nlp)%vpp  = 0.
!e        atp(nlp)%wpp  = 0.
        atp(nlp)%mass = 0. ! for now
        atp(nlp)%ppm  = 0. ! for now

     enddo

  enddo

  write(io6,'(A,I0)') " LP particles added this timestep: ", nlp - nlp0


end subroutine lpsource

!==============================================================================

subroutine normdist(anorm, nent)

  implicit none (external, type)

  integer, intent(in) :: nent
  real, intent(inout) :: anorm(nent)

  integer :: mident, ip

  real :: rnenti, cnorm, xnorm, znorm

  ! This routine generates a table ANORM of normally-distributed
  ! numbers.  The number of table entries is set by the user in
  ! NENT which must be set to a prime number and must be less than
  ! or equal to MAXENT.

  mident = nent / 2 + 1
  rnenti = 1.0 / real(nent)
  cnorm = 1.0 / sqrt(2.0 * 3.14159265)
  anorm(mident) = 0.0
  anorm(mident-1) = -rnenti / 0.39894
  do ip = mident+1,nent
     xnorm = 1.5 * anorm(ip-1) - 0.5 * anorm(ip-2)
     znorm = cnorm * exp(-0.5 * xnorm * xnorm)
     anorm(ip) = anorm(ip-1) + rnenti / znorm
     anorm(nent+1-ip) = -anorm(ip)
  enddo

end subroutine normdist

!==============================================================================

real function ramran(idum)

  implicit none (external, type)

  ! random number generator with [0,1] uniform distribution
  ! by Knuth subtractive method

  integer, parameter :: mbig = 1000000000, mseed = 161803398, mz = 0
  real, parameter    :: fac = 1. / real(mbig)

  integer       :: idum
  integer       :: i, ii, k, mj, mk
  integer, save :: iff = 0, inext, inextp, ma(55)

  if (idum < 0 .or. iff == 0) then
     iff = 1
     mj = MSEED - abs(idum)
     mj = mod(mj,MBIG)
     ma(55) = mj
     mk = 1
     do i = 1,54
        ii = mod(21 * i, 55)
        ma(ii) = mk
        mk = mj - mk
        if (mk < mz) mk = mk + mbig
        mj = ma(ii)
     enddo

     do k = 1,4
        do i = 1,55
           ma(i) = ma(i) - ma(1 + mod(i + 30, 55))
           if (ma(i) < mz) ma(i) = ma(i) + mbig
        enddo
     enddo
     inext = 0
     inextp = 31
     idum = 1
  endif

  inext = inext + 1
  if (inext == 56) inext = 1
  inextp = inextp + 1
  if (inextp == 56) inextp = 1
  mj = ma(inext) - ma(inextp)
  if (mj < mz) mj = mj + mbig
  ma(inext) = mj
  ramran = mj * fac

end function ramran
