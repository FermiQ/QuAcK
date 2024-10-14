subroutine UGW_phBSE(exchange_kernel,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS, &
                     S,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,cW,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: spin_conserved 
  logical,intent(in)            :: spin_flip

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: cW(nBas,nBas,nspin)
  double precision,intent(in)   :: eW(nBas,nspin)
  double precision,intent(in)   :: eGW(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: isp_W

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)

  integer                       :: nS_aa,nS_bb,nS_sc
  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY_BSE(:,:)
  double precision,allocatable  :: XmY_BSE(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  nS_ab = (nO(1) - nC(1))*(nV(2) - nR(2))
  nS_ba = (nO(2) - nC(2))*(nV(1) - nR(1))
  nS_sf = nS_ab + nS_ba
  
  allocate(OmRPA(nS_sc),XpY_RPA(nS_sc,nS_sc),XmY_RPA(nS_sc,nS_sc),rho_RPA(nBas,nBas,nS_sc,nspin))

!-----!
! TDA !
!-----!

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated in phBSE!'
    write(*,*)
  end if

!--------------------------!
! Spin-conserved screening !
!--------------------------!

  isp_W = 1
  EcRPA = 0d0

  ! Compute spin-conserved RPA screening 

  call phULR(isp_W,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
             eW,ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

  call UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

!----------------------------!
! Spin-conserved excitations !
!----------------------------!

  if(spin_conserved) then

    ispin = 1
    EcBSE(ispin) = 0d0

    allocate(OmBSE(nS_sc),XpY_BSE(nS_sc,nS_sc),XmY_BSE(nS_sc,nS_sc))

    ! Compute spin-conserved BSE excitation energies

    call phULR(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
               eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)
    call print_excitation_energies('phBSE@GW@UHF','spin-conserved',nS_sc,OmBSE)
    call phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nS_aa,nS_bb,nS_sc,dipole_int_aa,dipole_int_bb, &
                                  cW,S,OmBSE,XpY_BSE,XmY_BSE)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) &
      call UGW_phBSE_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nS_aa,nS_bb,nS_sc,nS_sc,    &
                                          eW,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb, &
                                          OmRPA,rho_RPA,OmBSE,XpY_BSE,XmY_BSE)

    deallocate(OmBSE,XpY_BSE,XmY_BSE)

  end if

!-----------------------!
! Spin-flip excitations !
!-----------------------!

 if(spin_flip) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Memory allocation
  
    allocate(OmBSE(nS_sf),XpY_BSE(nS_sf,nS_sf),XmY_BSE(nS_sf,nS_sf))

    ! Compute spin-flip BSE excitation energies

    call phULR(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,nS_sc,1d0, &
               eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcBSE(ispin),       & 
               OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@UHF','spin-flip',nS_sf,OmBSE)
    call phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nS_ab,nS_ba,nS_sf,dipole_int_aa,dipole_int_bb, &
                                  cW,S,OmBSE,XpY_BSE,XmY_BSE)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) &
      call UGW_phBSE_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nS_ab,nS_ba,nS_sf,nS_sc,    &
                                          eW,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb, &
                                          OmRPA,rho_RPA,OmBSE,XpY_BSE,XmY_BSE)

    deallocate(OmBSE,XpY_BSE,XmY_BSE)
 
  end if

  
! Scale properly correlation energy if exchange is included in interaction kernel

  if(exchange_kernel) then

    EcBSE(1) = 0.5d0*EcBSE(1)
    EcBSE(2) = 0.5d0*EcBSE(2)

  else

    EcBSE(2) = 0.0d0

  end if

end subroutine 
