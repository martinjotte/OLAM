module matrix

contains

!==========================================================================================!
! THESE SUBROUTINES COPIED FROM ED2/UTILS/NUMUTILS.F90 AND RENAMED HERE                    !
!                                                                                          !
! These subroutines will solve the linear system AA . X = Y for given AA and Y, using      !
! the Gaussian elimination method with partial pivoting and back-substitution.             !
!------------------------------------------------------------------------------------------!

subroutine matrix8_NxN(nsiz,AA,Y,X,sing)

  use consts_coms, only: r8
  implicit none

  !----- Arguments. ----------------------------------------------------------------------!
  integer                       , intent(in)  :: nsiz  ! matrix and vector size
  real(r8), dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
  real(r8), dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
  real(r8), dimension(nsiz)     , intent(out) :: X     ! unknown vector
  logical                       , intent(out) :: sing  ! The matrix was singular     [T|F]

  !----- Local variables. ----------------------------------------------------------------!
  real(r8), dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
  real(r8), dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
  real(r8), dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
  real(r8), dimension(nsiz)                   :: diagm1 ! The inverse of the pivot values
  real(r8)                                    :: multip ! Multiplier
  integer                                     :: r      ! Row index
  integer                                     :: b      ! Row below index
  integer                                     :: c      ! Column index
  integer                                     :: p      ! Pivot index
  real(r8)                                    :: dumsca ! Dummy scalar, for row swapping

  !----- Local parameters. ---------------------------------------------------------------!
  real(r8)                      , parameter   :: tinyoff=1.e-20
  !---------------------------------------------------------------------------------------!

  !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
  EE(:,:) = AA(:,:)
  Z (:)   = Y (:)

  !---------------------------------------------------------------------------------------!
  !     We initialise X with a huge, non-sense value, which will become the answer when   !
  ! the matrix is singular.                                                               !
  !---------------------------------------------------------------------------------------!
  X (:)   = -huge(1.)

  !----- We first assume that everything will be fine. -----------------------------------!
  sing    = .false.

  !---------------------------------------------------------------------------------------!
  ! 1. Main elimination loop, done row by row.                                            !
  !---------------------------------------------------------------------------------------!
  elimloop: do r = 1, nsiz-1

     !------ 1a. Finding the largest element, which will become our pivot ----------------!
     p = (r-1) + maxloc(abs(EE(r:nsiz,r)),dim=1)

     !------------------------------------------------------------------------------------!
     ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
     !     singular or almost singular, and we cannot solve it, so we switch the flag and !
     !     return.                                                                        !
     !------------------------------------------------------------------------------------!
     if (abs(EE(p,r)) < tinyoff) then
        sing = .true.
        return
     end if

     !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
     if (p /= r) then
        dumvec(r:nsiz) = EE(r,r:nsiz)
        dumsca         = Z(r)
        EE(r,r:nsiz)   = EE(p,r:nsiz)
        Z(r)           = Z(p)
        EE(p,r:nsiz)   = dumvec(r:nsiz)
        Z(p)           = dumsca
     end if

     ! Store the inverse of the pivot to avoid an extra divide
     diagm1(r) = 1.0_r8 / ee(r,r)

     !------------------------------------------------------------------------------------!
     ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
     !      zero (we won't compute that, but they will be.).                              !
     !------------------------------------------------------------------------------------!
     belowloop: do b=r+1,nsiz
        multip = EE(b,r) * diagm1(r)
        EE(b,r:nsiz) = EE(b,r:nsiz) - multip * EE(r,r:nsiz)
        Z(b)         = Z(b)         - multip * Z(r)
     end do belowloop
  end do elimloop

  !---------------------------------------------------------------------------------------!
  ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
  !    check the last pivot too.                                                          ! 
  !---------------------------------------------------------------------------------------!
  if (abs(EE(nsiz,nsiz)) < tinyoff) then
     sing = .true.
     return
  end if

  !---------------------------------------------------------------------------------------!
  ! 3. We now perform the back-substitution, to find the solution.                        !
  !---------------------------------------------------------------------------------------!
  X(nsiz) = Z(nsiz) / EE(nsiz,nsiz)
  if (nsiz > 1) X(nsiz-1) = (z(nsiz-1) - ee(nsiz-1,nsiz) * x(nsiz)) * diagm1(nsiz-1)

  backsubloop: do r=nsiz-2,1,-1
     X(r) = (Z(r) - dot_product(EE(r,r+1:nsiz),x(r+1:nsiz))) * diagm1(r)
  end do backsubloop

end subroutine matrix8_NxN

!==========================================================================================!

subroutine matrix8_2x2(AA,Y,X,sing)

  use consts_coms, only: r8
  implicit none

  integer, parameter :: nsiz = 2

  !----- Arguments. ----------------------------------------------------------------------!
  real(r8), dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
  real(r8), dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
  real(r8), dimension(nsiz)     , intent(out) :: X     ! unknown vector
  logical                       , intent(out) :: sing  ! The matrix was singular     [T|F]

  !----- Local variables. ----------------------------------------------------------------!
  real(r8), dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
  real(r8), dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
  real(r8), dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
  real(r8), dimension(nsiz)                   :: diagm1 ! The inverse of the pivot values
  real(r8)                                    :: multip ! Multiplier
  integer                                     :: r      ! Row index
  integer                                     :: b      ! Row below index
  integer                                     :: c      ! Column index
  integer                                     :: p      ! Pivot index
  real(r8)                                    :: dumsca ! Dummy scalar, for row swapping

  !----- Local parameters. ---------------------------------------------------------------!
  real(r8)                      , parameter   :: tinyoff=1.e-20
  !---------------------------------------------------------------------------------------!

  !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
  EE(:,:) = AA(:,:)
  Z (:)   = Y (:)

  !---------------------------------------------------------------------------------------!
  !     We initialise X with a huge, non-sense value, which will become the answer when   !
  ! the matrix is singular.                                                               !
  !---------------------------------------------------------------------------------------!
  !X (:)   = -huge(1.)

  !----- We first assume that everything will be fine. -----------------------------------!
  sing    = .false.

  !---------------------------------------------------------------------------------------!
  ! 1. Main elimination loop, done row by row.                                            !
  !---------------------------------------------------------------------------------------!

  !------ 1a. Finding the largest element, which will become our pivot ----------------!
  p = maxloc(abs(EE(1:nsiz,1)),dim=1)

  !------------------------------------------------------------------------------------!
  ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
  !     singular or almost singular, and we cannot solve it, so we switch the flag and !
  !     return.                                                                        !
  !------------------------------------------------------------------------------------!
  if (abs(EE(p,1)) < tinyoff) then
     sing = .true.
     return
  end if

  !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
  if (p /= 1) then
     dumvec(1:nsiz) = EE(1,1:nsiz)
     dumsca         = Z(1)
     EE(1,1:nsiz)   = EE(p,1:nsiz)
     Z(1)           = Z(p)
     EE(p,1:nsiz)   = dumvec(1:nsiz)
     Z(p)           = dumsca
  end if

  ! Store the inverse of the pivot to avoid an extra divide
  diagm1(1) = 1.0_r8 / ee(1,1)

  !------------------------------------------------------------------------------------!
  ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
  !      zero (we won't compute that, but they will be.).                              !
  !------------------------------------------------------------------------------------!
  multip = EE(nsiz,1) * diagm1(1)
  EE(nsiz,1:nsiz) = EE(nsiz,1:nsiz) - multip * EE(1,1:nsiz)
  Z(nsiz)         = Z(nsiz)         - multip * Z(1)

  !---------------------------------------------------------------------------------------!
  ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
  !    check the last pivot too.                                                          ! 
  !---------------------------------------------------------------------------------------!
  if (abs(EE(nsiz,nsiz)) < tinyoff) then
     sing = .true.
     return
  end if

  !---------------------------------------------------------------------------------------!
  ! 3. We now perform the back-substitution, to find the solution.                        !
  !---------------------------------------------------------------------------------------!
  X(nsiz)   = Z(nsiz) / EE(nsiz,nsiz)
  X(nsiz-1) = (z(nsiz-1) - ee(nsiz-1,nsiz) * x(nsiz)) * diagm1(nsiz-1)

end subroutine matrix8_2x2

!==========================================================================================!

subroutine matrix8_3x3(AA,Y,X,sing)

  use consts_coms, only: r8
  implicit none

  integer, parameter :: nsiz = 3

  !----- Arguments. ----------------------------------------------------------------------!
  real(r8), dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
  real(r8), dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
  real(r8), dimension(nsiz)     , intent(out) :: X     ! unknown vector
  logical                       , intent(out) :: sing  ! The matrix was singular     [T|F]

  !----- Local variables. ----------------------------------------------------------------!
  real(r8), dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
  real(r8), dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
  real(r8), dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
  real(r8), dimension(nsiz)                   :: diagm1 ! The inverse of the pivot values
  real(r8)                                    :: multip ! Multiplier
  integer                                     :: r      ! Row index
  integer                                     :: b      ! Row below index
  integer                                     :: c      ! Column index
  integer                                     :: p      ! Pivot index
  real(r8)                                    :: dumsca ! Dummy scalar, for row swapping

  !----- Local parameters. ---------------------------------------------------------------!
  real(r8)                      , parameter   :: tinyoff=1.e-20
  !---------------------------------------------------------------------------------------!

  !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
  EE(:,:) = AA(:,:)
  Z (:)   = Y (:)

  !---------------------------------------------------------------------------------------!
  !     We initialise X with a huge, non-sense value, which will become the answer when   !
  ! the matrix is singular.                                                               !
  !---------------------------------------------------------------------------------------!
  X (:)   = -huge(1.)

  !----- We first assume that everything will be fine. -----------------------------------!
  sing    = .false.

  !---------------------------------------------------------------------------------------!
  ! 1. Main elimination loop, done row by row.                                            !
  !---------------------------------------------------------------------------------------!
  elimloop: do r = 1, nsiz-1

     !------ 1a. Finding the largest element, which will become our pivot ----------------!
     p = (r-1) + maxloc(abs(EE(r:nsiz,r)),dim=1)

     !------------------------------------------------------------------------------------!
     ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
     !     singular or almost singular, and we cannot solve it, so we switch the flag and !
     !     return.                                                                        !
     !------------------------------------------------------------------------------------!
     if (abs(EE(p,r)) < tinyoff) then
        sing = .true.
        return
     end if

     !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
     if (p /= r) then
        dumvec(r:nsiz) = EE(r,r:nsiz)
        dumsca         = Z(r)
        EE(r,r:nsiz)   = EE(p,r:nsiz)
        Z(r)           = Z(p)
        EE(p,r:nsiz)   = dumvec(r:nsiz)
        Z(p)           = dumsca
     end if

     ! Store the inverse of the pivot to avoid an extra divide
     diagm1(r) = 1.0_r8 / ee(r,r)

     !------------------------------------------------------------------------------------!
     ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
     !      zero (we won't compute that, but they will be.).                              !
     !------------------------------------------------------------------------------------!
     belowloop: do b=r+1,nsiz
        multip = EE(b,r) * diagm1(r)
        EE(b,r:nsiz) = EE(b,r:nsiz) - multip * EE(r,r:nsiz)
        Z(b)         = Z(b)         - multip * Z(r)
     end do belowloop
  end do elimloop

  !---------------------------------------------------------------------------------------!
  ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
  !    check the last pivot too.                                                          ! 
  !---------------------------------------------------------------------------------------!
  if (abs(EE(nsiz,nsiz)) < tinyoff) then
     sing = .true.
     return
  end if

  !---------------------------------------------------------------------------------------!
  ! 3. We now perform the back-substitution, to find the solution.                        !
  !---------------------------------------------------------------------------------------!
  X(nsiz)   = Z(nsiz) / EE(nsiz,nsiz)
  X(nsiz-1) = (z(nsiz-1) - ee(nsiz-1,nsiz) * x(nsiz)) * diagm1(nsiz-1)

  backsubloop: do r=nsiz-2,1,-1
     X(r) = (Z(r) - dot_product(EE(r,r+1:nsiz),x(r+1:nsiz))) * diagm1(r)
  end do backsubloop

end subroutine matrix8_3x3

!==========================================================================================!

subroutine matrix8_4x4(AA,Y,X,sing)

  use consts_coms, only: r8
  implicit none

  integer, parameter :: nsiz = 4

  !----- Arguments. ----------------------------------------------------------------------!
  real(r8), dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
  real(r8), dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
  real(r8), dimension(nsiz)     , intent(out) :: X     ! unknown vector
  logical                       , intent(out) :: sing  ! The matrix was singular     [T|F]

  !----- Local variables. ----------------------------------------------------------------!
  real(r8), dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
  real(r8), dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
  real(r8), dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
  real(r8), dimension(nsiz)                   :: diagm1 ! The inverse of the pivot values
  real(r8)                                    :: multip ! Multiplier
  integer                                     :: r      ! Row index
  integer                                     :: b      ! Row below index
  integer                                     :: c      ! Column index
  integer                                     :: p      ! Pivot index
  real(r8)                                    :: dumsca ! Dummy scalar, for row swapping

  !----- Local parameters. ---------------------------------------------------------------!
  real(r8)                      , parameter   :: tinyoff=1.e-20
  !---------------------------------------------------------------------------------------!

  !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
  EE(:,:) = AA(:,:)
  Z (:)   = Y (:)

  !---------------------------------------------------------------------------------------!
  !     We initialise X with a huge, non-sense value, which will become the answer when   !
  ! the matrix is singular.                                                               !
  !---------------------------------------------------------------------------------------!
  X (:)   = -huge(1.)

  !----- We first assume that everything will be fine. -----------------------------------!
  sing    = .false.

  !---------------------------------------------------------------------------------------!
  ! 1. Main elimination loop, done row by row.                                            !
  !---------------------------------------------------------------------------------------!
  elimloop: do r = 1, nsiz-1

     !------ 1a. Finding the largest element, which will become our pivot ----------------!
     p = (r-1) + maxloc(abs(EE(r:nsiz,r)),dim=1)

     !------------------------------------------------------------------------------------!
     ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
     !     singular or almost singular, and we cannot solve it, so we switch the flag and !
     !     return.                                                                        !
     !------------------------------------------------------------------------------------!
     if (abs(EE(p,r)) < tinyoff) then
        sing = .true.
        return
     end if

     !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
     if (p /= r) then
        dumvec(r:nsiz) = EE(r,r:nsiz)
        dumsca         = Z(r)
        EE(r,r:nsiz)   = EE(p,r:nsiz)
        Z(r)           = Z(p)
        EE(p,r:nsiz)   = dumvec(r:nsiz)
        Z(p)           = dumsca
     end if

     ! Store the inverse of the pivot to avoid an extra divide
     diagm1(r) = 1.0_r8 / ee(r,r)

     !------------------------------------------------------------------------------------!
     ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
     !      zero (we won't compute that, but they will be.).                              !
     !------------------------------------------------------------------------------------!
     belowloop: do b=r+1,nsiz
        multip = EE(b,r) * diagm1(r)
        EE(b,r:nsiz) = EE(b,r:nsiz) - multip * EE(r,r:nsiz)
        Z(b)         = Z(b)         - multip * Z(r)
     end do belowloop
  end do elimloop

  !---------------------------------------------------------------------------------------!
  ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
  !    check the last pivot too.                                                          ! 
  !---------------------------------------------------------------------------------------!
  if (abs(EE(nsiz,nsiz)) < tinyoff) then
     sing = .true.
     return
  end if

  !---------------------------------------------------------------------------------------!
  ! 3. We now perform the back-substitution, to find the solution.                        !
  !---------------------------------------------------------------------------------------!
  X(nsiz)   = Z(nsiz) / EE(nsiz,nsiz)
  X(nsiz-1) = (z(nsiz-1) - ee(nsiz-1,nsiz) * x(nsiz)) * diagm1(nsiz-1)

  backsubloop: do r=nsiz-2,1,-1
     X(r) = (Z(r) - dot_product(EE(r,r+1:nsiz),x(r+1:nsiz))) * diagm1(r)
  end do backsubloop

end subroutine matrix8_4x4

!==========================================================================================!

subroutine matrix_NxN(nsiz,AA,Y,X,sing)
  implicit none

  !----- Arguments. ----------------------------------------------------------------------!
  integer                       , intent(in)  :: nsiz  ! matrix and vector size
  real,     dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
  real,     dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
  real,     dimension(nsiz)     , intent(out) :: X     ! unknown vector
  logical                       , intent(out) :: sing  ! The matrix was singular     [T|F]

  !----- Local variables. ----------------------------------------------------------------!
  real,     dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
  real,     dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
  real,     dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
  real,     dimension(nsiz)                   :: diagm1 ! The inverse of the pivot values
  real                                        :: multip ! Multiplier
  integer                                     :: r      ! Row index
  integer                                     :: b      ! Row below index
  integer                                     :: c      ! Column index
  integer                                     :: p      ! Pivot index
  real                                        :: dumsca ! Dummy scalar, for row swapping

  !----- Local parameters. ---------------------------------------------------------------!
  real                          , parameter   :: tinyoff=1.e-16
  !---------------------------------------------------------------------------------------!

  !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
  EE(:,:) = AA(:,:)
  Z (:)   = Y (:)

  !---------------------------------------------------------------------------------------!
  !     We initialise X with a huge, non-sense value, which will become the answer when   !
  ! the matrix is singular.                                                               !
  !---------------------------------------------------------------------------------------!
  X (:)   = -huge(1.)

  !----- We first assume that everything will be fine. -----------------------------------!
  sing    = .false.

  !---------------------------------------------------------------------------------------!
  ! 1. Main elimination loop, done row by row.                                            !
  !---------------------------------------------------------------------------------------!
  elimloop: do r = 1, nsiz-1

     !------ 1a. Finding the largest element, which will become our pivot ----------------!
     p = (r-1) + maxloc(abs(EE(r:nsiz,r)),dim=1)

     !------------------------------------------------------------------------------------!
     ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
     !     singular or almost singular, and we cannot solve it, so we switch the flag and !
     !     return.                                                                        !
     !------------------------------------------------------------------------------------!
     if (abs(EE(p,r)) < tinyoff) then
        sing = .true.
        return
     end if

     !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
     if (p /= r) then
        dumvec(r:nsiz) = EE(r,r:nsiz)
        dumsca         = Z(r)
        EE(r,r:nsiz)   = EE(p,r:nsiz)
        Z(r)           = Z(p)
        EE(p,r:nsiz)   = dumvec(r:nsiz)
        Z(p)           = dumsca
     end if

     ! Store the inverse of the pivot to avoid an extra divide
     diagm1(r) = 1.0 / ee(r,r)

     !------------------------------------------------------------------------------------!
     ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
     !      zero (we won't compute that, but they will be.).                              !
     !------------------------------------------------------------------------------------!
     belowloop: do b=r+1,nsiz
        multip = EE(b,r) * diagm1(r)
        EE(b,r:nsiz) = EE(b,r:nsiz) - multip * EE(r,r:nsiz)
        Z(b)         = Z(b)         - multip * Z(r)
     end do belowloop
  end do elimloop

  !---------------------------------------------------------------------------------------!
  ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
  !    check the last pivot too.                                                          ! 
  !---------------------------------------------------------------------------------------!
  if (abs(EE(nsiz,nsiz)) < tinyoff) then
     sing = .true.
     return
  end if

  !---------------------------------------------------------------------------------------!
  ! 3. We now perform the back-substitution, to find the solution.                        !
  !---------------------------------------------------------------------------------------!
  X(nsiz) = Z(nsiz) / EE(nsiz,nsiz)
  if (nsiz > 1) X(nsiz-1) = (z(nsiz-1) - ee(nsiz-1,nsiz) * x(nsiz)) * diagm1(nsiz-1)

  backsubloop: do r=nsiz-2,1,-1
     X(r) = (Z(r) - dot_product(EE(r,r+1:nsiz),x(r+1:nsiz))) * diagm1(r)
  end do backsubloop

end subroutine matrix_NxN

!==========================================================================================!

subroutine matrix_2x2(AA,Y,X,sing)
  implicit none

  integer, parameter :: nsiz = 2

  !----- Arguments. ----------------------------------------------------------------------!
  real,     dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
  real,     dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
  real,     dimension(nsiz)     , intent(out) :: X     ! unknown vector
  logical                       , intent(out) :: sing  ! The matrix was singular     [T|F]

  !----- Local variables. ----------------------------------------------------------------!
  real,     dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
  real,     dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
  real,     dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
  real,     dimension(nsiz)                   :: diagm1 ! The inverse of the pivot values
  real                                        :: multip ! Multiplier
  integer                                     :: r      ! Row index
  integer                                     :: b      ! Row below index
  integer                                     :: c      ! Column index
  integer                                     :: p      ! Pivot index
  real                                        :: dumsca ! Dummy scalar, for row swapping

  !----- Local parameters. ---------------------------------------------------------------!
  real                          , parameter   :: tinyoff=1.e-16
  !---------------------------------------------------------------------------------------!

  !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
  EE(:,:) = AA(:,:)
  Z (:)   = Y (:)

  !---------------------------------------------------------------------------------------!
  !     We initialise X with a huge, non-sense value, which will become the answer when   !
  ! the matrix is singular.                                                               !
  !---------------------------------------------------------------------------------------!
  !X (:)   = -huge(1.)

  !----- We first assume that everything will be fine. -----------------------------------!
  sing    = .false.

  !---------------------------------------------------------------------------------------!
  ! 1. Main elimination loop, done row by row.                                            !
  !---------------------------------------------------------------------------------------!

  !------ 1a. Finding the largest element, which will become our pivot ----------------!
  p = maxloc(abs(EE(1:nsiz,1)),dim=1)

  !------------------------------------------------------------------------------------!
  ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
  !     singular or almost singular, and we cannot solve it, so we switch the flag and !
  !     return.                                                                        !
  !------------------------------------------------------------------------------------!
  if (abs(EE(p,1)) < tinyoff) then
     sing = .true.
     return
  end if

  !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
  if (p /= 1) then
     dumvec(1:nsiz) = EE(1,1:nsiz)
     dumsca         = Z(1)
     EE(1,1:nsiz)   = EE(p,1:nsiz)
     Z(1)           = Z(p)
     EE(p,1:nsiz)   = dumvec(1:nsiz)
     Z(p)           = dumsca
  end if

  ! Store the inverse of the pivot to avoid an extra divide
  diagm1(1) = 1.0 / ee(1,1)

  !------------------------------------------------------------------------------------!
  ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
  !      zero (we won't compute that, but they will be.).                              !
  !------------------------------------------------------------------------------------!
  multip = EE(nsiz,1) * diagm1(1)
  EE(nsiz,1:nsiz) = EE(nsiz,1:nsiz) - multip * EE(1,1:nsiz)
  Z(nsiz)         = Z(nsiz)         - multip * Z(1)

  !---------------------------------------------------------------------------------------!
  ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
  !    check the last pivot too.                                                          ! 
  !---------------------------------------------------------------------------------------!
  if (abs(EE(nsiz,nsiz)) < tinyoff) then
     sing = .true.
     return
  end if

  !---------------------------------------------------------------------------------------!
  ! 3. We now perform the back-substitution, to find the solution.                        !
  !---------------------------------------------------------------------------------------!
  X(nsiz)   = Z(nsiz) / EE(nsiz,nsiz)
  X(nsiz-1) = (z(nsiz-1) - ee(nsiz-1,nsiz) * x(nsiz)) * diagm1(nsiz-1)

end subroutine matrix_2x2

!==========================================================================================!

subroutine matrix_3x3(AA,Y,X,sing)
  implicit none

  integer, parameter :: nsiz = 3

  !----- Arguments. ----------------------------------------------------------------------!
  real,     dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
  real,     dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
  real,     dimension(nsiz)     , intent(out) :: X     ! unknown vector
  logical                       , intent(out) :: sing  ! The matrix was singular     [T|F]

  !----- Local variables. ----------------------------------------------------------------!
  real,     dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
  real,     dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
  real,     dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
  real,     dimension(nsiz)                   :: diagm1 ! The inverse of the pivot values
  real                                        :: multip ! Multiplier
  integer                                     :: r      ! Row index
  integer                                     :: b      ! Row below index
  integer                                     :: c      ! Column index
  integer                                     :: p      ! Pivot index
  real                                        :: dumsca ! Dummy scalar, for row swapping

  !----- Local parameters. ---------------------------------------------------------------!
  real                          , parameter   :: tinyoff=1.e-16
  !---------------------------------------------------------------------------------------!

  !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
  EE(:,:) = AA(:,:)
  Z (:)   = Y (:)

  !---------------------------------------------------------------------------------------!
  !     We initialise X with a huge, non-sense value, which will become the answer when   !
  ! the matrix is singular.                                                               !
  !---------------------------------------------------------------------------------------!
  X (:)   = -huge(1.)

  !----- We first assume that everything will be fine. -----------------------------------!
  sing    = .false.

  !---------------------------------------------------------------------------------------!
  ! 1. Main elimination loop, done row by row.                                            !
  !---------------------------------------------------------------------------------------!
  elimloop: do r = 1, nsiz-1

     !------ 1a. Finding the largest element, which will become our pivot ----------------!
     p = (r-1) + maxloc(abs(EE(r:nsiz,r)),dim=1)

     !------------------------------------------------------------------------------------!
     ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
     !     singular or almost singular, and we cannot solve it, so we switch the flag and !
     !     return.                                                                        !
     !------------------------------------------------------------------------------------!
     if (abs(EE(p,r)) < tinyoff) then
        sing = .true.
        return
     end if

     !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
     if (p /= r) then
        dumvec(r:nsiz) = EE(r,r:nsiz)
        dumsca         = Z(r)
        EE(r,r:nsiz)   = EE(p,r:nsiz)
        Z(r)           = Z(p)
        EE(p,r:nsiz)   = dumvec(r:nsiz)
        Z(p)           = dumsca
     end if

     ! Store the inverse of the pivot to avoid an extra divide
     diagm1(r) = 1.0 / ee(r,r)

     !------------------------------------------------------------------------------------!
     ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
     !      zero (we won't compute that, but they will be.).                              !
     !------------------------------------------------------------------------------------!
     belowloop: do b=r+1,nsiz
        multip = EE(b,r) * diagm1(r)
        EE(b,r:nsiz) = EE(b,r:nsiz) - multip * EE(r,r:nsiz)
        Z(b)         = Z(b)         - multip * Z(r)
     end do belowloop
  end do elimloop

  !---------------------------------------------------------------------------------------!
  ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
  !    check the last pivot too.                                                          ! 
  !---------------------------------------------------------------------------------------!
  if (abs(EE(nsiz,nsiz)) < tinyoff) then
     sing = .true.
     return
  end if

  !---------------------------------------------------------------------------------------!
  ! 3. We now perform the back-substitution, to find the solution.                        !
  !---------------------------------------------------------------------------------------!
  X(nsiz)   = Z(nsiz) / EE(nsiz,nsiz)
  X(nsiz-1) = (z(nsiz-1) - ee(nsiz-1,nsiz) * x(nsiz)) * diagm1(nsiz-1)

  backsubloop: do r=nsiz-2,1,-1
     X(r) = (Z(r) - dot_product(EE(r,r+1:nsiz),x(r+1:nsiz))) * diagm1(r)
  end do backsubloop

end subroutine matrix_3x3

!==========================================================================================!

subroutine matrix_4x4(AA,Y,X,sing)
  implicit none

  integer, parameter :: nsiz = 4

  !----- Arguments. ----------------------------------------------------------------------!
  real,     dimension(nsiz,nsiz), intent(in)  :: AA    ! matrix
  real,     dimension(nsiz)     , intent(in)  :: Y     ! right-hand side vector
  real,     dimension(nsiz)     , intent(out) :: X     ! unknown vector
  logical                       , intent(out) :: sing  ! The matrix was singular     [T|F]

  !----- Local variables. ----------------------------------------------------------------!
  real,     dimension(nsiz,nsiz)              :: EE     ! Copy of AA, for elimination.
  real,     dimension(nsiz)                   :: Z      ! Copy of Y, for scaling
  real,     dimension(nsiz)                   :: dumvec ! Dummy vector, for row swapping
  real,     dimension(nsiz)                   :: diagm1 ! The inverse of the pivot values
  real                                        :: multip ! Multiplier
  integer                                     :: r      ! Row index
  integer                                     :: b      ! Row below index
  integer                                     :: c      ! Column index
  integer                                     :: p      ! Pivot index
  real                                        :: dumsca ! Dummy scalar, for row swapping

  !----- Local parameters. ---------------------------------------------------------------!
  real                          , parameter   :: tinyoff=1.e-16
  !---------------------------------------------------------------------------------------!

  !----- First thing, we copy AA to EE and Y to Z. ---------------------------------------!
  EE(:,:) = AA(:,:)
  Z (:)   = Y (:)

  !---------------------------------------------------------------------------------------!
  !     We initialise X with a huge, non-sense value, which will become the answer when   !
  ! the matrix is singular.                                                               !
  !---------------------------------------------------------------------------------------!
  X (:)   = -huge(1.)

  !----- We first assume that everything will be fine. -----------------------------------!
  sing    = .false.

  !---------------------------------------------------------------------------------------!
  ! 1. Main elimination loop, done row by row.                                            !
  !---------------------------------------------------------------------------------------!
  elimloop: do r = 1, nsiz-1

     !------ 1a. Finding the largest element, which will become our pivot ----------------!
     p = (r-1) + maxloc(abs(EE(r:nsiz,r)),dim=1)

     !------------------------------------------------------------------------------------!
     ! 1b. Check the pivot and make sure it is a good one.  If not, then this matrix is   !
     !     singular or almost singular, and we cannot solve it, so we switch the flag and !
     !     return.                                                                        !
     !------------------------------------------------------------------------------------!
     if (abs(EE(p,r)) < tinyoff) then
        sing = .true.
        return
     end if

     !----- 1c. If the best pivot is not the current row, we must swap them. -------------!
     if (p /= r) then
        dumvec(r:nsiz) = EE(r,r:nsiz)
        dumsca         = Z(r)
        EE(r,r:nsiz)   = EE(p,r:nsiz)
        Z(r)           = Z(p)
        EE(p,r:nsiz)   = dumvec(r:nsiz)
        Z(p)           = dumsca
     end if

     ! Store the inverse of the pivot to avoid an extra divide
     diagm1(r) = 1.0 / ee(r,r)

     !------------------------------------------------------------------------------------!
     ! 1d.  Eliminate rows below, everything to the left of the (,r) column will become   !
     !      zero (we won't compute that, but they will be.).                              !
     !------------------------------------------------------------------------------------!
     belowloop: do b=r+1,nsiz
        multip = EE(b,r) * diagm1(r)
        EE(b,r:nsiz) = EE(b,r:nsiz) - multip * EE(r,r:nsiz)
        Z(b)         = Z(b)         - multip * Z(r)
     end do belowloop
  end do elimloop

  !---------------------------------------------------------------------------------------!
  ! 2. We may be unlucky and discover that the matrix is singular at the last line, so we !
  !    check the last pivot too.                                                          ! 
  !---------------------------------------------------------------------------------------!
  if (abs(EE(nsiz,nsiz)) < tinyoff) then
     sing = .true.
     return
  end if

  !---------------------------------------------------------------------------------------!
  ! 3. We now perform the back-substitution, to find the solution.                        !
  !---------------------------------------------------------------------------------------!
  X(nsiz)   = Z(nsiz) / EE(nsiz,nsiz)
  X(nsiz-1) = (z(nsiz-1) - ee(nsiz-1,nsiz) * x(nsiz)) * diagm1(nsiz-1)

  backsubloop: do r=nsiz-2,1,-1
     X(r) = (Z(r) - dot_product(EE(r,r+1:nsiz),x(r+1:nsiz))) * diagm1(r)
  end do backsubloop

end subroutine matrix_4x4

end module matrix
