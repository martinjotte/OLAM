Module mem_rayf

  implicit none

! Memory for Rayleigh friction layer

  real :: rayf_zmin
  real :: rayf_expon
  real :: rayf_fact

  real :: rayfw_zmin
  real :: rayfw_expon
  real :: rayfw_fact

  real :: rayfdiv_zmin
  real :: rayfdiv_expon
  real :: rayfdiv_fact

  real :: rayfmix_zmin
  real :: rayfmix_expon
  real :: rayfmix_fact

  real, allocatable :: rayf_cof   (:)
  real, allocatable :: rayf_cofw  (:)
  real, allocatable :: rayf_cofmix(:)

  real, allocatable :: vc03d(:,:)

  integer :: krayf_bot
  integer :: krayfw_bot
  integer :: krayfdiv_bot
  integer :: krayfmix_bot

  logical :: dorayf    = .false.
  logical :: dorayfw   = .false.
  logical :: dorayfdiv = .false.
  logical :: dorayfmix = .false.

Contains

  subroutine rayf_init()

! Initialize Rayleigh friction vertical profile coefficients

    use misc_coms,  only: dtsm, initial
    use mem_grid,   only: mza, mva, lpv, zm, zt
    use oname_coms, only: nl
    use mem_basic,  only: vc

    implicit none

    integer :: k, iv
    real    :: distimi, distim0, dti

    dti = 1.0 / real(dtsm)

    allocate( rayf_cof   (mza) )
    allocate( rayf_cofw  (mza) )
    allocate( rayf_cofmix(mza) )

    rayf_cof    = 0.0
    rayf_cofw   = 0.0
    rayf_cofmix = 0.0

    krayf_bot    = mza + 1
    krayfw_bot   = mza + 1
    krayfdiv_bot = mza + 1
    krayfmix_bot = mza + 1

! RAYF COEFFICIENT FOR THIL AND VMC (ONLY FOR HORIZ. HOMOG. INITIALIZATION)

    if (rayf_fact > 1.e-7 .and. initial == 1) then

       do k = 2, mza
          if (zt(k) > rayf_zmin) then
             krayf_bot = k
             exit
          endif
       enddo

       if (krayf_bot <= mza) then

          dorayf  = .true.
          distimi = rayf_fact * dti

          do k = krayf_bot, mza
             rayf_cof(k) = distimi   &
                  * ((zt(k) - rayf_zmin) / (zm(mza) - rayf_zmin)) ** rayf_expon
          enddo

       endif
    endif

! RAYF coefficient for WMC damping

    if (rayfw_fact > 1.e-7) then

       do k = 2, mza-1
          if (zm(k) > rayfw_zmin) then
             krayfw_bot = k
             exit
          endif
       enddo

       if (krayfw_bot < mza) then

          dorayfw = .true.
          distimi = rayfw_fact * dti

          do k = krayfw_bot, mza-1
             rayf_cofw(k) = distimi   &
                  * ((zm(k) - rayfw_zmin) / (zm(mza) - rayfw_zmin)) ** rayfw_expon
          enddo

       endif
    endif

! RAYF coefficients for Horiz Divergence are now computed in vort_div_damp.f90

    if (rayfdiv_fact > 1.e-7) then

       do k = 2, mza
          if (zt(k) > rayfdiv_zmin) then
             krayfdiv_bot = k
             exit
          endif
       enddo

       if (krayfdiv_bot <= mza) then
          dorayfdiv = .true.
       endif

    endif

! RAYF coefficient for horizontal velocity (VC) mixing

    if (rayfmix_fact > 1.e-7) then

       do k = 2, mza-1
          if (zm(k) > rayfmix_zmin) then
             krayfmix_bot = k
             exit
          endif
       enddo

       if (krayfmix_bot < mza) then

          dorayfmix = .true.

          do k = krayfmix_bot, mza-1
             rayf_cofmix(k) = 0.5 * rayfmix_fact &
                  * ((zm(k) - rayfmix_zmin) / (zm(mza) - rayfw_zmin)) ** rayfmix_expon
          enddo

       endif
    endif

! For horizontally homogeneous run using Rayleigh friction, allocate and fill
! vc03d with initial horiz velocity that the model will relax towards

    if (dorayf) then
       allocate( vc03d(mza,mva) )

       !$omp parallel do private(k)
       do iv = 2, mva
          do k = lpv(iv), mza
             vc03d(k,iv) = vc(k,iv)
          enddo
       enddo
       !$omp end parallel do
    endif

  end subroutine rayf_init



  subroutine rayf_mix_top_vxe( iw, vmxet, vmyet, vmzet )

    use mem_basic,   only: vxe, vye, vze, rho
    use mem_grid,    only: mza, arw0, lpw

    implicit none

    integer, intent(in   ) :: iw
    real,    intent(inout) :: vmxet(mza) ! xe-mom tend scaled by volume
    real,    intent(inout) :: vmyet(mza) ! ye-mom tend scaled by volume
    real,    intent(inout) :: vmzet(mza) ! ze-mom tend scaled by volume

    real    :: vxflux(mza)
    real    :: vyflux(mza)
    real    :: vzflux(mza)
    real    :: fact
    integer :: ka, k

    ka = max(krayfmix_bot, lpw(iw))

    vxflux(mza)  = 0.0
    vxflux(ka-1) = 0.0

    vyflux(mza)  = 0.0
    vyflux(ka-1) = 0.0

    vzflux(mza)  = 0.0
    vzflux(ka-1) = 0.0

    do k = ka, mza-1
       fact = arw0(iw) * rayf_cofmix(k) * real(rho(k+1,iw) + rho(k,iw))

       vxflux(k) = fact * (vxe(k,iw) - vxe(k+1,iw))
       vyflux(k) = fact * (vye(k,iw) - vye(k+1,iw))
       vzflux(k) = fact * (vze(k,iw) - vze(k+1,iw))
    enddo

    do k = ka, mza
       vmxet(k) = vmxet(k) + vxflux(k-1) - vxflux(k)
       vmyet(k) = vmyet(k) + vyflux(k-1) - vyflux(k)
       vmzet(k) = vmzet(k) + vzflux(k-1) - vzflux(k)
    enddo

  end subroutine rayf_mix_top_vxe


end module mem_rayf
