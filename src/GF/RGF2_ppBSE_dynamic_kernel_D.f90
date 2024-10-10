subroutine RGF2_ppBSE_dynamic_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOO,lambda,ERI,eGF,OmBSE,KD_dyn,ZD_dyn)

! Compute the resonant part of the dynamic BSE@GF2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: OmBSE
  
! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: i,j,k,l,m
  integer                       :: e
  integer                       :: ij,kl

! Output variables

  double precision,intent(out)  :: KD_dyn(nOO,nOO)
  double precision,intent(out)  :: ZD_dyn(nOO,nOO)

! Initialization

  KD_dyn(:,:) = 0d0
  ZD_dyn(:,:) = 0d0

! Second-order correlation kernel for the block D of the singlet manifold

  if(ispin == 1) then

    ij = 0
    do i=nC+1,nO
      do j=i,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k,nO
            kl = kl + 1
  
            do m=nC+1,nO
              do e=nO+1,nBas-nR
     
                dem = - OmBSE + eGF(k) - eGF(e) + eGF(m) + eGF(j)
                num = 2d0*ERI(i,m,k,e)*ERI(e,j,m,l) - ERI(i,m,k,e)*ERI(e,j,l,m)  &
                     -    ERI(i,m,e,k)*ERI(e,j,m,l) - ERI(i,m,e,k)*ERI(e,j,l,m)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(k) - eGF(e) + eGF(m) + eGF(i)
                num = 2d0*ERI(j,m,k,e)*ERI(e,i,m,l) - ERI(j,m,k,e)*ERI(e,i,l,m)  &
                     -    ERI(j,m,e,k)*ERI(e,i,m,l) - ERI(j,m,e,k)*ERI(e,i,l,m)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(l) - eGF(e) + eGF(m) + eGF(i)
                num = 2d0*ERI(i,e,k,m)*ERI(m,j,e,l) - ERI(i,e,k,m)*ERI(m,j,l,e)  &
                     -    ERI(i,e,m,k)*ERI(m,j,e,l) - ERI(i,e,m,k)*ERI(m,j,l,e)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(l) - eGF(e) + eGF(m) + eGF(j)
                num = 2d0*ERI(j,e,k,m)*ERI(m,i,e,l) - ERI(j,e,k,m)*ERI(m,i,l,e)  &
                     -    ERI(j,e,m,k)*ERI(m,i,e,l) - ERI(j,e,m,k)*ERI(m,i,l,e)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

            KD_dyn(ij,kl) = lambda*KD_dyn(ij,kl)/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))
            
          end do
        end do

      end do
    end do

  end if

! Second-order correlation kernel for the block D of the triplet manifold

  if(ispin == 2) then

    ij = 0
    do i=nC+1,nO
      do j=i+1,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
  
            do m=nC+1,nO
              do e=nO+1,nBas-nR
     
                dem = - OmBSE + eGF(k) - eGF(e) + eGF(m) + eGF(j)
                num = 2d0*ERI(i,m,k,e)*ERI(e,j,m,l) - ERI(i,m,k,e)*ERI(e,j,l,m)  &
                     -    ERI(i,m,e,k)*ERI(e,j,m,l) - ERI(i,m,e,k)*ERI(e,j,l,m)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(k) - eGF(e) + eGF(m) + eGF(i)
                num = 2d0*ERI(j,m,k,e)*ERI(e,i,m,l) - ERI(j,m,k,e)*ERI(e,i,l,m)  &
                     -    ERI(j,m,e,k)*ERI(e,i,m,l) - ERI(j,m,e,k)*ERI(e,i,l,m)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) - num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(l) - eGF(e) + eGF(m) + eGF(i)
                num = 2d0*ERI(i,e,k,m)*ERI(m,j,e,l) - ERI(i,e,k,m)*ERI(m,j,l,e)  &
                     -    ERI(i,e,m,k)*ERI(m,j,e,l) - ERI(i,e,m,k)*ERI(m,j,l,e)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(l) - eGF(e) + eGF(m) + eGF(j)
                num = 2d0*ERI(j,e,k,m)*ERI(m,i,e,l) - ERI(j,e,k,m)*ERI(m,i,l,e)  &
                     -    ERI(j,e,m,k)*ERI(m,i,e,l) - ERI(j,e,m,k)*ERI(m,i,l,e)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) - num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

          end do
        end do

      end do
    end do

  end if

end subroutine 
