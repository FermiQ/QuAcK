subroutine GGF2_phBSE(TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,EcBSE)

! Compute the second-order Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: dRPA = .false.
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: A_sta(:,:)
  double precision,allocatable  :: B_sta(:,:)
  double precision,allocatable  :: KA_sta(:,:)
  double precision,allocatable  :: KB_sta(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE

! Memory allocation

  allocate(OmBSE(nS),XpY(nS,nS),XmY(nS,nS),A_sta(nS,nS),KA_sta(nS,nS))
  allocate(B_sta(nS,nS),KB_sta(nS,nS))

  EcBSE = 0d0

               call phGLR_A(dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGF,ERI,A_sta)
  if(.not.TDA) call phGLR_B(dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,B_sta)

  ! Compute static kernel

               call GGF2_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,KA_sta)
  if(.not.TDA) call GGF2_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,KB_sta)

               A_sta(:,:) = A_sta(:,:) + KA_sta(:,:)
  if(.not.TDA) B_sta(:,:) = B_sta(:,:) + KB_sta(:,:)

  ! Compute phBSE@GF2 excitation energies

  call phGLR(TDA,nS,A_sta,B_sta,EcBSE,OmBSE,XpY,XmY)
  call print_excitation_energies('phBSE@GGF2','generalized',nS,OmBSE)
  call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY,XmY)

  ! Compute dynamic correction for BSE via perturbation theory

! if(dBSE) &
!   call GGF2_phBSE_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,KA_sta,KB_sta,OmBSE,XpY,XmY)

end subroutine 
