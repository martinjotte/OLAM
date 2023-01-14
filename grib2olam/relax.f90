subroutine relax(nxin,nyin,field,mask)

  implicit none

  integer, intent(in)    :: nxin, nyin
  real,    intent(inout) :: field(nxin,nyin)
  integer, intent(in)    :: mask(nxin,nyin)

  real            :: change
  real, parameter :: pi = 3.141592654
  real, parameter :: tolerance = 2.e-4
  real            :: r, omega, fac1, fac2, fac3, utemp, umax

  real, allocatable :: u(:,:)
  integer            :: nx, ny, m, i, j, imin

  ny = nyin
  if (mod(nxin,2) == 0) then
     nx = nxin
  else
     nx = nxin + 1
  endif

  allocate(u(0:nx+1,ny))

  u(1:nxin,:) = field
  if (mod(nx,2) /= 0) u(nx,:) = u(nx-1,:)

  u(0,   :) = u(nx,:)
  u(nx+1,:) = u(1, :)

  r = 0.5 * (cos(pi/real(nx)) + cos(pi/real(ny)))
  omega = min(2.0 / (1.0 + sqrt(1.0 - r*r)), 1.99)
  fac1 = 1.0 - omega
  fac2 = omega * 0.25
  fac3 = omega / 3.

  umax = max(maxval(abs(field),mask=(mask==0)), tolerance)

  do m = 1, nx * ny

     change = 0.0

     ! red squares
     do j = 1, ny

        imin = 1 + mod(j,2)
        do i = imin, nx, 2

           if (mask(i,j) == 0) cycle

           if (j == 1) then
              utemp = fac1 * u(i,j) + fac3 * (u(i-1,j) + u(i+1,j) + u(i,j+1))
           elseif (j == ny) then
              utemp = fac1 * u(i,j) + fac3 * (u(i-1,j) + u(i+1,j) + u(i,j-1))
           else
              utemp = fac1 * u(i,j) + fac2 * (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))
           endif

           change = max(abs(utemp - u(i,j)), change)
           u(i,j) = utemp

        enddo
     enddo

     u(0,   :) = u(nx,:)
     u(nx+1,:) = u(1, :)

     ! black squares
     do j = 1, ny

        imin = 1 + mod(j+1,2)
        do i = imin, nx, 2

           if (mask(i,j) == 0) cycle

           if (j == 1) then
              utemp = fac1 * u(i,j) + fac3 * (u(i-1,j) + u(i+1,j) + u(i,j+1))
           elseif (j == ny) then
              utemp = fac1 * u(i,j) + fac3 * (u(i-1,j) + u(i+1,j) + u(i,j-1))
           else
              utemp = fac1 * u(i,j) + fac2 * (u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1))
           endif

           change = max(abs(utemp - u(i,j)), change)
           u(i,j) = utemp

        enddo
     enddo

     u(0,   :) = u(nx,:)
     u(nx+1,:) = u(1, :)

     if (mod(m,50) == 0) write(*,*) m, change / umax

     if (change / umax <= tolerance) exit

  enddo

  write(*,'(A,I0,A,/)') " Converged after ", m, " iterations"
  field = u(1:nxin,:)

end subroutine relax
