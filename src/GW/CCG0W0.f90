subroutine CCG0W0(maxSCF,thresh,nBas,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! CC-based GW module

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Local variables

  integer                       :: p,q
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

  integer                       :: nSCF
  double precision              :: Conv

  double precision,allocatable  :: OVVO(:,:,:,:)
  double precision,allocatable  :: VOOV(:,:,:,:)

  double precision,allocatable  :: delta_2h1p(:,:,:,:)
  double precision,allocatable  :: delta_2p1h(:,:,:,:)

  double precision,allocatable  :: V_2h1p(:,:,:,:)
  double precision,allocatable  :: V_2p1h(:,:,:,:)

  double precision,allocatable  :: r_2h1p(:,:,:,:)
  double precision,allocatable  :: r_2p1h(:,:,:,:)

  double precision,allocatable  :: t_2h1p(:,:,:,:)
  double precision,allocatable  :: t_2p1h(:,:,:,:)

  double precision,allocatable  :: x_2h1p(:,:)
  double precision,allocatable  :: x_2p1h(:,:)

  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: SigGW(:,:)
  double precision,allocatable  :: cGW(:,:)
  double precision,allocatable  :: Z(:)

  integer,allocatable           :: order(:)

! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* CC-based G0W0 Calculation *'
  write(*,*)'*****************************'
  write(*,*)

! Create integral batches

  allocate(OVVO(nO,nV,nV,nO),VOOV(nV,nO,nO,nV))

  OVVO(:,:,:,:) = ERI(   1:nO   ,nO+1:nOrb,nO+1:nOrb,   1:nO  )
  VOOV(:,:,:,:) = ERI(nO+1:nOrb ,   1:nO  ,   1:nO  ,nO+1:nOrb)
 
! Form energy denominator and guess amplitudes

  allocate(delta_2h1p(nO,nO,nV,nOrb),delta_2p1h(nO,nV,nV,nOrb))
  allocate(V_2h1p(nOrb,nO,nO,nV),V_2p1h(nOrb,nO,nV,nV))
  allocate(t_2h1p(nO,nO,nV,nOrb),t_2p1h(nO,nV,nV,nOrb))
  allocate(x_2h1p(nOrb,nOrb),x_2p1h(nOrb,nOrb))

  do k=nC+1,nO
    do l=nC+1,nO
      do c=1,nV-nR
        do p=nC+1,nOrb-nR

          V_2h1p(p,k,l,c) = sqrt(2d0)*ERI(p,nO+c,k,l)

        end do
      end do
    end do
  end do

  do k=nC+1,nO
    do c=1,nV-nR
      do d=1,nV-nR
        do p=nC+1,nOrb-nR

          V_2p1h(p,k,c,d) = sqrt(2d0)*ERI(p,k,nO+d,nO+c)

        end do
      end do
    end do
  end do

! Initialization

  allocate(r_2h1p(nO,nO,nV,nOrb),r_2p1h(nO,nV,nV,nOrb))
  allocate(eGW(nOrb),SigGW(nOrb,nOrb),cGW(nOrb,nOrb),Z(nOrb))
  allocate(order(nOrb))

  Conv   = 1d0
  nSCF   =  0
  eGW(:) = eHF(:)

  t_2h1p(:,:,:,:) = 0d0
  t_2p1h(:,:,:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------'
  write(*,*)'| CC-based G0W0 calculation                  |'
  write(*,*)'----------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','HOMO','|','LUMO','|','Conv','|'
  write(*,*)'----------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

   ! Increment 

    nSCF = nSCF + 1

    ! Compute energy differences

    do i=nC+1,nO
      do j=nC+1,nO
        do a=1,nV-nR
          do p=nC+1,nOrb-nR
 
            delta_2h1p(i,j,a,p) = eGW(i) + eGW(j) - eGW(nO+a) - eHF(p)
 
          end do
        end do
      end do
    end do

    do i=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR
          do p=nC+1,nOrb-nR
 
            delta_2p1h(i,a,b,p) = eGW(nO+a) + eGW(nO+b) - eGW(i) - eHF(p)
 
          end do
        end do
      end do
    end do

    !  Compute intermediates

    x_2h1p(:,:) = 0d0 

    do p=nC+1,nOrb-nR
      do q=nC+1,nOrb-nR

        do k=nC+1,nO
          do l=nC+1,nO
            do c=1,nV-nR
      
              x_2h1p(p,q) = x_2h1p(p,q) + V_2h1p(q,k,l,c)*t_2h1p(k,l,c,p)
      
            end do
          end do
        end do

      end do
    end do
   
    x_2p1h(:,:) = 0d0 

    do p=nC+1,nOrb-nR
      do q=nC+1,nOrb-nR

        do k=nC+1,nO
          do c=1,nV-nR
            do d=1,nV-nR
      
              x_2p1h(p,q) = x_2p1h(p,q) + V_2p1h(q,k,c,d)*t_2p1h(k,c,d,p)
     
          end do
        end do
      end do

      end do
    end do
   
    ! Compute residual for 2h1p sector

    do i=nC+1,nO
      do j=nC+1,nO
        do a=1,nV-nR

          do p=nC+1,nOrb-nR

            r_2h1p(i,j,a,p) = V_2h1p(p,i,j,a) + delta_2h1p(i,j,a,p)*t_2h1p(i,j,a,p)

            do k=nC+1,nO
              do c=1,nV-nR

                r_2h1p(i,j,a,p) = r_2h1p(i,j,a,p) - 2d0*OVVO(j,c,a,k)*t_2h1p(i,k,c,p)

              end do
            end do

            do q=nC+1,nOrb-nR

              r_2h1p(i,j,a,p) = r_2h1p(i,j,a,p) - t_2h1p(i,j,a,q)*x_2h1p(p,q) - t_2h1p(i,j,a,q)*x_2p1h(p,q)

            end do
 
          end do

        end do
      end do
    end do

    ! Compute residual for 2p1h sector

    do i=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR

          do p=nC+1,nOrb-nR

            r_2p1h(i,a,b,p) = V_2p1h(p,i,a,b) + delta_2p1h(i,a,b,p)*t_2p1h(i,a,b,p)

            do k=nC+1,nO
              do c=1,nV-nR

                r_2p1h(i,a,b,p) = r_2p1h(i,a,b,p) + 2d0*VOOV(a,k,i,c)*t_2p1h(k,c,b,p)

              end do
            end do

            do q=nC+1,nOrb-nR

              r_2p1h(i,a,b,p) = r_2p1h(i,a,b,p) - t_2p1h(i,a,b,q)*x_2h1p(p,q) - t_2p1h(i,a,b,q)*x_2p1h(p,q)

            end do
 
          end do

        end do
      end do
    end do
 
    !  Check convergence 

    Conv = max(maxval(abs(r_2h1p)),maxval(abs(r_2p1h)))
  
    ! Update amplitudes

    t_2h1p(:,:,:,:) = t_2h1p(:,:,:,:) - r_2h1p(:,:,:,:)/delta_2h1p(:,:,:,:)
    t_2p1h(:,:,:,:) = t_2p1h(:,:,:,:) - r_2p1h(:,:,:,:)/delta_2p1h(:,:,:,:)

    ! Compute self-energy

    SigGW(:,:) = 0d0

    do p=nC+1,nOrb-nR

      SigGW(p,p) = SigGW(p,p) + eHF(p)

      do q=nC+1,nOrb-nR

        do i=nC+1,nO
          do j=nC+1,nO
            do a=1,nV-nR

              SigGW(p,q) = SigGW(p,q) + V_2h1p(p,i,j,a)*t_2h1p(i,j,a,q)
 
            end do
          end do
        end do

        do i=nC+1,nO
          do a=1,nV-nR
            do b=1,nV-nR

              SigGW(p,q) = SigGW(p,q) + V_2p1h(p,i,a,b)*t_2p1h(i,a,b,q)
 
            end do
          end do
        end do

      end do
    end do

    !  Diagonalize non-Hermitian matrix

    call diagonalize_general_matrix(nOrb,SigGW,eGW,cGW)

    do p=1,nOrb
      order(p) = p
    end do

    call quick_sort(eGW,order,nOrb)
    call set_order(cGW,order,nOrb,nOrb)

    ! Renormalization factor

    Z(:) = 1d0

    ! Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',eGW(nO)*HaToeV,'|',eGW(nO+1)*HaToeV,'|',Conv,'|'

  end do
  write(*,*)'----------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  end if

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'  CCGW calculation  '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nOrb
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',(eGW(p)-eHF(p))*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
  end do
  write(*,*)'-------------------------------------------------------------------------------'

end subroutine 
