subroutine evRGTeh(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, & 
                   singlet,triplet,linearize,eta,regularize,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform self-consistent eigenvalue-only ehGT calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  logical                       :: print_T = .true.
  logical                       :: dRPA = .false.
  logical                       :: linear_mixing
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond
  double precision              :: Conv
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision              :: alpha
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rhoL(:,:,:)
  double precision,allocatable  :: rhoR(:,:,:)

! Hello world

  write(*,*)
  write(*,*)'**********************************'
  write(*,*)'* Restricted evRGTeh Calculation *'
  write(*,*)'**********************************'
  write(*,*)

! TDA for T

  if(TDA_T) then 
    write(*,*) 'Tamm-Dancoff approximation for eh T-matrix!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Linear mixing

  linear_mixing = .false.
  alpha = 0.2d0

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),eGT(nOrb),eOld(nOrb),Z(nOrb),Sig(nOrb),Om(nS),XpY(nS,nS),XmY(nS,nS), & 
           rhoL(nOrb,nOrb,nS),rhoR(nOrb,nOrb,nS),error_diis(nOrb,max_diis),e_diis(nOrb,max_diis))

! Initialization

  nSCF            = 0
  ispin           = 2
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGT(:)          = eHF(:)
  eOld(:)         = eGT(:)
  Z(:)            = 1d0
  rcond           = 0d0
  
!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

   ! Compute screening

    call phLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGT,ERI,Aph)
    if(.not.TDA_T) call phLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

    call phLR(TDA_T,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

    if(print_T) call print_excitation_energies('phRPA@RHF','triplet',nS,Om)

   ! Compute spectral weights

    call RGTeh_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,XmY,rhoL,rhoR)

    ! Compute correlation part of the self-energy 

    if(regularize) call GTeh_regularization(nOrb,nC,nO,nV,nR,nS,eGT,Om,rhoL,rhoR)

    call RGTeh_self_energy_diag(eta,nOrb,nC,nO,nV,nR,nS,eGT,Om,rhoL,rhoR,EcGM,Sig,Z)

    ! Solve the quasi-particle equation

    eGT(:) = eHF(:) + Sig(:)

    if(linearize) then

       write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
       write(*,*)

       eGT(:) = eGT(:)

    else

       write(*,*) ' *** Quasiparticle energies obtained by root search *** '
       write(*,*)

       call RGTeh_QP_graph(eta,nOrb,nC,nO,nV,nR,nS,eHF,Om,rhoL,rhoR,eOld,eOld,eGT,Z)

    end if

    ! Convergence criteria

    Conv = maxval(abs(eGT - eOld))

    ! Print results

    call print_evRGTeh(nOrb,nO,nSCF,Conv,eHF,ENuc,ERHF,Sig,Z,eGT,EcRPA,EcGM)

    ! Linear mixing or DIIS extrapolation

    if(linear_mixing) then
 
      eGT(:) = alpha*eGT(:) + (1d0 - alpha)*eOld(:)
 
    else

      n_diis = min(n_diis+1,max_diis)
      if(abs(rcond) > 1d-7) then
        call DIIS_extrapolation(rcond,nOrb,nOrb,n_diis,error_diis,e_diis,eGT-eOld,eGT)
      else
        n_diis = 0
      end if

    end if

    ! Save quasiparticles energy for next cycle

    eOld(:) = eGT(:)

    ! Increment

    nSCF = nSCF + 1

  end do
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF+1) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  end if

! Deallocate memory

  deallocate(eOld,Z,Sig,Om,XpY,XmY,rhoL,rhoR,error_diis,e_diis)

! Perform BSE calculation

! if(BSE) then

!   call Bethe_Salpeter(BSE2,TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eGW,eGW,EcBSE)

!   if(exchange_kernel) then

!     EcBSE(1) = 0.5d0*EcBSE(1)
!     EcBSE(2) = 1.5d0*EcBSE(2)

!   end if

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy (singlet) =',EcBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy (triplet) =',EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy           =',EcBSE(1) + EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

!   if(doACFDT) then

!     write(*,*) '------------------------------------------------------'
!     write(*,*) 'Adiabatic connection version of BSE correlation energy'
!     write(*,*) '------------------------------------------------------'
!     write(*,*)

!     if(doXBS) then

!       write(*,*) '*** scaled screening version (XBS) ***'
!       write(*,*)

!     end if

!     call ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,eGW,eGW,EcBSE)

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (singlet) =',EcBSE(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (triplet) =',EcBSE(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy           =',EcBSE(1) + EcBSE(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

! end if

! if(ppBSE) then

!   call Bethe_Salpeter_pp(TDA_W,TDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcppBSE)

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy (singlet) =',EcppBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy (triplet) =',3d0*EcppBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy =',EcppBSE(1) + 3d0*EcppBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 total energy =',ENuc + ERHF + EcppBSE(1) + 3d0*EcppBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

!   nOrb2 = 2*nOrb
!   nO2   = 2*nO
!   nV2   = 2*nV
!   nC2   = 2*nC
!   nR2   = 2*nR
!   nS2   = nO2*nV2
!
!   allocate(seHF(nOrb2),seGW(nOrb2),sERI(nOrb2,nOrb2,nOrb2,nOrb2))
!
!   call spatial_to_spin_MO_energy(nOrb,eHF,nOrb2,seHF)
!   call spatial_to_spin_MO_energy(nOrb,eGW,nOrb2,seGW)
!   call spatial_to_spin_ERI(nOrb,ERI,nOrb2,sERI)
!
!   call Bethe_Salpeter_pp_so(TDA_W,TDA,singlet,triplet,eta,nOrb2,nC2,nO2,nV2,nR2,nS2,sERI,dipole_int,seHF,seGW,EcppBSE)

! end if

! Testing zone

  if(dotest) then

    call dump_test_value('R','evGTeh correlation energy',EcRPA)
    call dump_test_value('R','evGTeh HOMO energy',eGT(nO))
    call dump_test_value('R','evGTeh LUMO energy',eGT(nO+1))

  end if

end subroutine 
