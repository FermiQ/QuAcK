subroutine sort_ppRPA(nOO,nVV,Om,Z,Om1,X1,Y1,Om2,X2,Y2)

! Compute the metric matrix for pp-RPA

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: Om(nOO+nVV)
  double precision,intent(in)   :: Z(nOO+nVV,nOO+nVV)
  
! Local variables

  integer                       :: pq,ab,ij
! integer                       :: deg1,ab_start,ab_end
! integer                       :: deg2,ij_start,ij_end
  double precision,allocatable  :: M(:,:)
  double precision,allocatable  :: Z1(:,:)
  double precision,allocatable  :: Z2(:,:)
  double precision,allocatable  :: S1(:,:)
  double precision,allocatable  :: S2(:,:)
  double precision,allocatable  :: O1(:,:)
  double precision,allocatable  :: O2(:,:)  
  double precision,allocatable  :: tmp1(:,:)
  double precision,allocatable  :: tmp2(:,:)

  integer,allocatable           :: order1(:)
  integer,allocatable           :: order2(:)

! Output variables

  double precision,intent(out)  :: Om1(nVV)
  double precision,intent(out)  :: X1(nVV,nVV)
  double precision,intent(out)  :: Y1(nOO,nVV)
  double precision,intent(out)  :: Om2(nOO)
  double precision,intent(out)  :: X2(nVV,nOO)
  double precision,intent(out)  :: Y2(nOO,nOO)


! Memory allocation

  allocate(M(nOO+nVV,nOO+nVV),Z1(nOO+nVV,nVV),Z2(nOO+nVV,nOO),order1(nVV),order2(nOO))

! Initializatiom

  Om1(:)  = 0d0
  X1(:,:) = 0d0
  Y1(:,:) = 0d0

  Om2(:)  = 0d0
  X2(:,:) = 0d0
  Y2(:,:) = 0d0

! Compute metric 

  M(:,:) = 0d0
  
  do ab=1,nVV
    M(ab,ab) = 1d0
  end do

  do ij=1,nOO
    M(nVV+ij,nVV+ij) = -1d0
  end do

! Start sorting eigenvectors

  ab = 0
  ij = 0

  do pq=1,nOO+nVV

    if(Om(pq) > 0d0) then 

      ab = ab + 1
      Om1(ab) = Om(pq)
      Z1(1:nOO+nVV,ab) = Z(1:nOO+nVV,pq)

    else

      ij = ij + 1
      Om2(ij) = Om(pq)
      Z2(1:nOO+nVV,ij) = Z(1:nOO+nVV,pq)

    end if

  end do

  if(minval(Om1) < 0d0 .or. ab /= nVV) call print_warning('You may have instabilities in pp-RPA!!')
  if(maxval(Om2) > 0d0 .or. ij /= nOO) call print_warning('You may have instabilities in pp-RPA!!')

  if(nVV > 0) then 

    do ab=1,nVV
      order1(ab) = ab
    end do

    call quick_sort(Om1,order1,nVV)
    call set_order(Z1,order1,nOO+nVV,nVV)

  end if

  if(nOO > 0) then

    do ij=1,nOO
      order2(ij) = ij
    end do

    call quick_sort(Om2,order2,nOO)
    call set_order(Z2,order2,nOO+nVV,nOO)

  end if
 
  allocate(S1(nVV,nVV),S2(nOO,nOO),O1(nVV,nVV),O2(nOO,nOO))
  allocate(tmp1(nOO+nVV,nVV),tmp2(nOO+nVV,nOO))

  if(nVV > 0) call dgemm ('N', 'N', nOO+nVV, nVV, nOO+nVV, 1d0,  M, nOO+nVV, Z1, nOO+nVV, 0d0, tmp1, nOO+nVV)
  if(nVV > 0) call dgemm ('T', 'N', nVV    , nVV, nOO+nVV, 1d0, Z1, nOO+nVV, tmp1, nOO+nVV, 0d0, S1, nVV)
  if(nOO > 0) call dgemm ('N', 'N', nOO+nVV, nOO, nOO+nVV, 1d0,  M, nOO+nVV, -1d0*Z2, nOO+nVV, 0d0, tmp2, nOO+nVV)
  if(nOO > 0) call dgemm ('T', 'N', nOO    , nOO, nOO+nVV, 1d0, Z2, nOO+nVV, tmp2, nOO+nVV, 0d0, S2, nOO)

! S1 = + matmul(transpose(Z1),matmul(M,Z1))
! S2 = - matmul(transpose(Z2),matmul(M,Z2))

  if(nVV > 0) call orthogonalize_matrix(1,nVV,S1,O1)
  if(nOO > 0) call orthogonalize_matrix(1,nOO,S2,O2)

  if(nVV > 0) call dgemm ('N', 'N', nOO+nVV,nVV,nVV, 1d0, Z1, nOO+nVV, O1, nVV,0d0, tmp1, nOO+nVV)
  Z1 = tmp1
  if(nOO > 0) call dgemm ('N', 'N', nOO+nVV,nOO,nOO, 1d0, Z2, nOO+nVV, O2, nOO,0d0, tmp2, nOO+nVV)
  Z2 = tmp2

! Z1 = matmul(Z1,O1)
! Z2 = matmul(Z2,O2)

! Define submatrices X1, Y1, X2, & Y2

  X1(1:nVV,1:nVV) = + Z1(    1:    nVV,1:nVV)
  Y1(1:nOO,1:nVV) = - Z1(nVV+1:nOO+nVV,1:nVV)

  X2(1:nVV,1:nOO) = + Z2(    1:    nVV,1:nOO)
  Y2(1:nOO,1:nOO) = - Z2(nVV+1:nOO+nVV,1:nOO)

! call matout(nVV,nVV,X1)
! call matout(nOO,nVV,Y1)

! call matout(nVV,nOO,X2)
! call matout(nOO,nOO,Y2)

end subroutine 
