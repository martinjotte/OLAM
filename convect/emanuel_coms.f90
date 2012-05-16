module emanuel_coms

  integer, allocatable ::  nent(:)
  real,    allocatable ::  u1ent(:,:), u2ent(:,:), u3ent(:,:), traent(:,:,:), tratm(:)
  real,    allocatable ::  u1p(:), u2p(:), u3p(:), trap(:,:)
  real,    allocatable ::  m(:), mp(:), ment(:,:), qent(:,:), elij(:,:)
  real,    allocatable ::  sij(:,:), tvp(:), tv(:), water(:)
  real,    allocatable ::  qp(:), ep(:), wt(:), evap(:), clw(:)
  real,    allocatable ::  sigp(:), tp(:), told(:), cpn(:)
  real,    allocatable ::  lv(:), lvcp(:), h(:), hp(:), gz(:), hm(:)
  real,    allocatable ::  qcond(:), nqcond(:), wa(:), ma(:), siga(:), ax(:)

contains

  subroutine alloc_eman(na, ntra)

    implicit none
    integer, intent(in) :: ntra, na

    allocate( nent(na) )
    allocate( u1ent(na,na), u2ent(na,na), u3ent(na,na), traent(na,na,ntra), tratm(na) )
    allocate( u1p(na), u2p(na), u3p(na), trap(na,ntra) )
    allocate( m(na), mp(na), ment(na,na), qent(na,na), elij(na,na) )
    allocate( sij(na,na), tvp(na), tv(na), water(na) )
    allocate( qp(na), ep(na), wt(na), evap(na), clw(na) )
    allocate( sigp(na), tp(na), told(na), cpn(na) )
    allocate( lv(na), lvcp(na), h(na), hp(na), gz(na), hm(na) )
    allocate( qcond(na), nqcond(na), wa(na), ma(na), siga(na), ax(na) )

  end subroutine alloc_eman

end module emanuel_coms
