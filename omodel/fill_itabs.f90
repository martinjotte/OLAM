subroutine mdloopf(init,im,j1,j2,j3,j4,j5,j6)

  use mem_delaunay, only: itab_md

  implicit none

  character(1), intent(in) :: init
  integer,      intent(in) :: im
  integer,      intent(in) :: j1,j2,j3,j4,j5,j6

  if (init == 'f') itab_md(im)%loop(:) = .false.

  if (j1 < 0) itab_md(im)%loop(abs(j1)) = .false.
  if (j2 < 0) itab_md(im)%loop(abs(j2)) = .false.
  if (j3 < 0) itab_md(im)%loop(abs(j3)) = .false.
  if (j4 < 0) itab_md(im)%loop(abs(j4)) = .false.
  if (j5 < 0) itab_md(im)%loop(abs(j5)) = .false.
  if (j6 < 0) itab_md(im)%loop(abs(j6)) = .false.

  if (j1 > 0) itab_md(im)%loop(j1) = .true.
  if (j2 > 0) itab_md(im)%loop(j2) = .true.
  if (j3 > 0) itab_md(im)%loop(j3) = .true.
  if (j4 > 0) itab_md(im)%loop(j4) = .true.
  if (j5 > 0) itab_md(im)%loop(j5) = .true.
  if (j6 > 0) itab_md(im)%loop(j6) = .true.

end subroutine mdloopf

!===============================================================================

subroutine udloopf(init,iu,j1,j2,j3,j4,j5,j6)

  use mem_delaunay, only: itab_ud

  implicit none

  character(1), intent(in) :: init
  integer,      intent(in) :: iu
  integer,      intent(in) :: j1,j2,j3,j4,j5,j6

  if (init == 'f') itab_ud(iu)%loop(:) = .false.

  if (j1 < 0) itab_ud(iu)%loop(abs(j1)) = .false.
  if (j2 < 0) itab_ud(iu)%loop(abs(j2)) = .false.
  if (j3 < 0) itab_ud(iu)%loop(abs(j3)) = .false.
  if (j4 < 0) itab_ud(iu)%loop(abs(j4)) = .false.
  if (j5 < 0) itab_ud(iu)%loop(abs(j5)) = .false.
  if (j6 < 0) itab_ud(iu)%loop(abs(j6)) = .false.

  if (j1 > 0) itab_ud(iu)%loop(j1) = .true.
  if (j2 > 0) itab_ud(iu)%loop(j2) = .true.
  if (j3 > 0) itab_ud(iu)%loop(j3) = .true.
  if (j4 > 0) itab_ud(iu)%loop(j4) = .true.
  if (j5 > 0) itab_ud(iu)%loop(j5) = .true.
  if (j6 > 0) itab_ud(iu)%loop(j6) = .true.

end subroutine udloopf

!===============================================================================

subroutine wdloopf(init,iw,j1,j2,j3,j4,j5,j6)

  use mem_delaunay, only: itab_wd

  implicit none

  character(1), intent(in) :: init
  integer,      intent(in) :: iw
  integer,      intent(in) :: j1,j2,j3,j4,j5,j6

  if (init == 'f') itab_wd(iw)%loop(:) = .false.

  if (j1 < 0) itab_wd(iw)%loop(abs(j1)) = .false.
  if (j2 < 0) itab_wd(iw)%loop(abs(j2)) = .false.
  if (j3 < 0) itab_wd(iw)%loop(abs(j3)) = .false.
  if (j4 < 0) itab_wd(iw)%loop(abs(j4)) = .false.
  if (j5 < 0) itab_wd(iw)%loop(abs(j5)) = .false.
  if (j6 < 0) itab_wd(iw)%loop(abs(j6)) = .false.

  if (j1 > 0) itab_wd(iw)%loop(j1) = .true.
  if (j2 > 0) itab_wd(iw)%loop(j2) = .true.
  if (j3 > 0) itab_wd(iw)%loop(j3) = .true.
  if (j4 > 0) itab_wd(iw)%loop(j4) = .true.
  if (j5 > 0) itab_wd(iw)%loop(j5) = .true.
  if (j6 > 0) itab_wd(iw)%loop(j6) = .true.

end subroutine wdloopf

!===============================================================================

subroutine mloopf(init,im,j1,j2,j3,j4,j5,j6)

  use mem_ijtabs, only: itab_m

  implicit none

  character(1), intent(in) :: init
  integer,      intent(in) :: im
  integer,      intent(in) :: j1,j2,j3,j4,j5,j6

  if (init == 'f') itab_m(im)%loop(:) = .false.

  if (j1 < 0) itab_m(im)%loop(abs(j1)) = .false.
  if (j2 < 0) itab_m(im)%loop(abs(j2)) = .false.
  if (j3 < 0) itab_m(im)%loop(abs(j3)) = .false.
  if (j4 < 0) itab_m(im)%loop(abs(j4)) = .false.
  if (j5 < 0) itab_m(im)%loop(abs(j5)) = .false.
  if (j6 < 0) itab_m(im)%loop(abs(j6)) = .false.

  if (j1 > 0) itab_m(im)%loop(j1) = .true.
  if (j2 > 0) itab_m(im)%loop(j2) = .true.
  if (j3 > 0) itab_m(im)%loop(j3) = .true.
  if (j4 > 0) itab_m(im)%loop(j4) = .true.
  if (j5 > 0) itab_m(im)%loop(j5) = .true.
  if (j6 > 0) itab_m(im)%loop(j6) = .true.

end subroutine mloopf

!===============================================================================

subroutine vloopf(init,iv,j1,j2,j3,j4,j5,j6)

  use mem_ijtabs, only: itab_v

  implicit none

  character(1), intent(in) :: init
  integer,      intent(in) :: iv
  integer,      intent(in) :: j1,j2,j3,j4,j5,j6

  if (init == 'f') itab_v(iv)%loop(:) = .false.

  if (j1 < 0) itab_v(iv)%loop(abs(j1)) = .false.
  if (j2 < 0) itab_v(iv)%loop(abs(j2)) = .false.
  if (j3 < 0) itab_v(iv)%loop(abs(j3)) = .false.
  if (j4 < 0) itab_v(iv)%loop(abs(j4)) = .false.
  if (j5 < 0) itab_v(iv)%loop(abs(j5)) = .false.
  if (j6 < 0) itab_v(iv)%loop(abs(j6)) = .false.

  if (j1 > 0) itab_v(iv)%loop(j1) = .true.
  if (j2 > 0) itab_v(iv)%loop(j2) = .true.
  if (j3 > 0) itab_v(iv)%loop(j3) = .true.
  if (j4 > 0) itab_v(iv)%loop(j4) = .true.
  if (j5 > 0) itab_v(iv)%loop(j5) = .true.
  if (j6 > 0) itab_v(iv)%loop(j6) = .true.

end subroutine vloopf

!===============================================================================

subroutine wloopf(init,iw,j1,j2,j3,j4,j5,j6)

  use mem_ijtabs, only: itab_w

  implicit none

  character(1), intent(in) :: init
  integer,      intent(in) :: iw
  integer,      intent(in) :: j1,j2,j3,j4,j5,j6

  if (init == 'f') itab_w(iw)%loop(:) = .false.

  if (j1 < 0) itab_w(iw)%loop(abs(j1)) = .false.
  if (j2 < 0) itab_w(iw)%loop(abs(j2)) = .false.
  if (j3 < 0) itab_w(iw)%loop(abs(j3)) = .false.
  if (j4 < 0) itab_w(iw)%loop(abs(j4)) = .false.
  if (j5 < 0) itab_w(iw)%loop(abs(j5)) = .false.
  if (j6 < 0) itab_w(iw)%loop(abs(j6)) = .false.

  if (j1 > 0) itab_w(iw)%loop(j1) = .true.
  if (j2 > 0) itab_w(iw)%loop(j2) = .true.
  if (j3 > 0) itab_w(iw)%loop(j3) = .true.
  if (j4 > 0) itab_w(iw)%loop(j4) = .true.
  if (j5 > 0) itab_w(iw)%loop(j5) = .true.
  if (j6 > 0) itab_w(iw)%loop(j6) = .true.

end subroutine wloopf
