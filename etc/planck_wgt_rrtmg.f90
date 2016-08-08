program planck_wgt
  implicit none

  ! This routine computes the overlap between arbitrary solar bands and
  ! RRTMg's bands. Specify the beginning and end of your band of interest,
  ! and the routine will report how to weight RRTMg's bands to produce your band

  ! UVA
  !real, parameter :: b1_beg = 0.320
  !real, parameter :: b1_end = 0.400

  ! UVB
  !real, parameter :: b1_beg = 0.290
  !real, parameter :: b1_end = 0.320

  ! UVC
  !real, parameter :: b1_beg = 0.100
  !real, parameter :: b1_end = 0.290

  ! PAR
  real, parameter :: b1_beg = 0.400
  real, parameter :: b1_end = 0.700

  ! CMAQ band 1
  !real, parameter :: b1_beg = .2910
  !real, parameter :: b1_end = .2983

  ! CMAQ band 2
  !real, parameter :: b1_beg = .2983
  !real, parameter :: b1_end = .3075

  ! CMAQ band 3
  !real, parameter :: b1_beg = .3075
  !real, parameter :: b1_end = .3125

  ! CMAQ band 4
  !real, parameter :: b1_beg = .3125
  !real, parameter :: b1_end = .3203

  ! CMAQ band 5
  !real, parameter :: b1_beg = .3203
  !real, parameter :: b1_end = .3450

  ! CMAQ band 6
  !real, parameter :: b1_beg = .3450
  !real, parameter :: b1_end = .4125

  ! CMAQ band 7
  !real, parameter :: b1_beg = .4125
  !real, parameter :: b1_end = .8500

  integer, parameter :: nrrtmg = 14

  real,    parameter :: b2_beg(nrrtmg) = (/ 3.077, 2.500, 2.150, 1.942,      &
       1.626, 1.299, 1.242, 0.778, 0.625, 0.442, 0.345, 0.263, 0.200, 3.846 /)

  real,    parameter :: b2_end(nrrtmg) = (/ 3.846, 3.077, 2.500, 2.150,      &
       1.942, 1.626, 1.299, 1.242, 0.778, 0.625, 0.442, 0.345, 0.263, 12.195/)

  integer, parameter :: nwavs = 1000
  integer            :: nwav1, nwav2

  real,    parameter :: bbtemp = 5800.0
  
  real               :: b1s, b1e
  real               :: b1, b2, frac1
  real               :: dwave, wave, sum1, sum2
  integer            :: i, n

  write(*,*) "RRTMg_Band ", " Planck_av ", " Lin_av"

  do n = 1, nrrtmg

     if ( (b1_beg < b2_beg(n) .and. b1_end <= b2_beg(n)) .or. &
          (b1_end > b2_end(n) .and. b1_beg >= b2_end(n)) ) then

        write(*,'(5x,I2,4x,2f10.6)') n, 0.0, 0.0
        
     elseif ((b1_beg <= b2_beg(n)) .and. (b1_end >= b2_end(n))) then

        write(*,'(5x,I2,4x,2f10.6)') n, 1.0, 1.0
        
     else

        b1s = max(b1_beg,b2_beg(n))
        b1e = min(b1_end,b2_end(n))

        b1 = b1e - b1s
        b2 = b2_end(n) - b2_beg(n)

        frac1 = b1 / max(b2,tiny(1.0))
        frac1 = min(frac1,1.0)

        nwav1 = nint(frac1 * real(nwavs))
        nwav1 = max(nwav1,1)
        nwav2 = nwavs

        ! Planck spectral radiance in band 1 overlap with band 2
        dwave = b1 / nwav1
        sum1  = 0.0
        do i = 1, nwav1
           wave = b1s + dwave * real(i-1) + 0.5 * dwave
           sum1 = sum1 + planck( wave, bbtemp ) * dwave
        enddo

        ! Planck spectral radiance in band 2
        dwave = b2 / nwav2
        sum2  = 0.0
        do i = 1, nwav2
           wave = b2_beg(n) + dwave * real(i-1) + 0.5 * dwave
           sum2 = sum2 + planck( wave, bbtemp ) * dwave
        enddo

        write(*,'(5x,I2,4x,2f10.6)') n, sum1 / sum2, b1 / b2

     endif

  enddo

contains

  real function planck(wave,temp)  ! W / m^2 / sr / micron
    implicit none
    real, parameter  :: c1 = 1.1910429e8
    real, parameter  :: c2 = 1.4387770e4
    real, intent(in) :: wave ! (microns)
    real, intent(in) :: temp ! (K)
    planck = c1 / ( wave**5 * (exp(c2/(wave*temp))-1.) )
  end function planck

end program planck_wgt
