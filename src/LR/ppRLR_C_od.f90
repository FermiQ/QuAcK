subroutine ppRLR_C_od(ispin,nBas,nC,nO,nV,nR,nOO,nVV,lambda,ERI,Cpp)

! Compute the C matrix of the pp channel (without diagonal term)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nOO,nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision,external     :: Kronecker_delta

  integer                       :: a,b,c,d,ab,cd

! Output variables

  double precision,intent(out)  :: Cpp(nVV,nVV)
 
! Build C matrix for the singlet manifold

  if(ispin == 1) then

    ab = 0
    do a=nO+1,nBas-nR
     do b=a,nBas-nR
        ab = ab + 1
        cd = 0
        do c=nO+1,nBas-nR
         do d=c,nBas-nR
            cd = cd + 1
 
            Cpp(ab,cd) = lambda*(ERI(a,b,c,d) + ERI(a,b,d,c))/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
 
          end do
        end do
      end do
    end do

  end if

! Build C matrix for the triplet or alpha-alpha manifold

  if(ispin == 2) then

    ab = 0
    do a=nO+1,nBas-nR
     do b=a+1,nBas-nR
        ab = ab + 1
        cd = 0
        do c=nO+1,nBas-nR
         do d=c+1,nBas-nR
            cd = cd + 1
 
            Cpp(ab,cd) = lambda*(ERI(a,b,c,d) - ERI(a,b,d,c))
 
          end do
        end do
      end do
    end do

  end if

! Build the alpha-beta block of the C matrix

  if(ispin == 3) then

    ab = 0
    do a=nO+1,nBas-nR
     do b=nO+1,nBas-nR
        ab = ab + 1
        cd = 0
        do c=nO+1,nBas-nR
         do d=nO+1,nBas-nR
            cd = cd + 1
 
            Cpp(ab,cd) = lambda*ERI(a,b,c,d)
 
          end do
        end do
      end do
    end do

  end if

end subroutine 
