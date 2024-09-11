subroutine SRG_self_energy(nBas,nC,nO,nV,nR,nS,e,Om,rho,EcGM,SigC,Z)

! Compute correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q,r
  integer                       :: m
  double precision              :: Dpim,Dqim,Dpam,Dqam,Diam
  double precision              :: t1,t2
  double precision              :: s
  
! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nBas,nBas)
  double precision,intent(out)  :: Z(nBas)

! SRG flow parameter 

  s = 500d0

! Initialize 

  SigC(:,:) = 0d0

!--------------------!
! SRG-GW self-energy !
!--------------------!

  ! Occupied part of the correlation self-energy

  call wall_time(t1)

  !$OMP PARALLEL &
  !$OMP SHARED(SigC,rho,s,nS,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(m,i,q,p,Dpim,Dqim) &
  !$OMP DEFAULT(NONE)
  !$OMP DO 
  do q=nC+1,nBas-nR
     do p=nC+1,nBas-nR
        do m=1,nS
           do i=nC+1,nO
              Dpim = e(p) - e(i) + Om(m)
              Dqim = e(q) - e(i) + Om(m)
              SigC(p,q) = SigC(p,q) + 2d0*rho(p,i,m)*rho(q,i,m)*(1d0-dexp(-s*Dpim*Dpim)*dexp(-s*Dqim*Dqim)) &
                   *(Dpim + Dqim)/(Dpim*Dpim + Dqim*Dqim)
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

! call wall_time(t2)
! print *, "first loop", (t2-t1)

! Virtual part of the correlation self-energy

 call wall_time(t1)
 !$OMP PARALLEL &
 !$OMP SHARED(SigC,rho,s,nS,nC,nO,nR,nBas,e,Om) &
 !$OMP PRIVATE(m,a,q,p,Dpam,Dqam) &
 !$OMP DEFAULT(NONE)
 !$OMP DO
 do q=nC+1,nBas-nR
    do p=nC+1,nBas-nR
       do m=1,nS
          do a=nO+1,nBas-nR
             Dpam = e(p) - e(a) - Om(m)
             Dqam = e(q) - e(a) - Om(m)
             SigC(p,q) = SigC(p,q) + 2d0*rho(p,a,m)*rho(q,a,m)*(1d0-exp(-s*Dpam*Dpam)*exp(-s*Dqam*Dqam)) &
                  *(Dpam + Dqam)/(Dpam*Dpam + Dqam*Dqam)
          end do
       end do
    end do
 end do
 !$OMP END DO
 !$OMP END PARALLEL

! call wall_time(t2)
! print *, "second loop", (t2-t1)
 

! Initialize

  Z(:)  = 0d0

  do p=nC+1,nBas-nR
     do i=nC+1,nO
        do m=1,nS
           Dpim = e(p) - e(i) + Om(m)
           Z(p) = Z(p)  - 2d0*rho(p,i,m)**2*(1d0-dexp(-2d0*s*Dpim*Dpim))/Dpim**2
        end do
     end do
  end do

  ! Virtual part of the correlation self-energy

  do p=nC+1,nBas-nR
     do a=nO+1,nBas-nR
        do m=1,nS
           Dpam = e(p) - e(a) - Om(m)
           Z(p) = Z(p)  - 2d0*rho(p,a,m)**2*(1d0-dexp(-2d0*s*Dpam*Dpam))/Dpam**2
        end do
     end do
  end do

! Compute renormalization factor from derivative of SigC

  Z(:) = 1d0/(1d0 - Z(:))

! Galitskii-Migdal correlation energy

  EcGM = 0d0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do m=1,nS
        Diam = e(a) - e(i) + Om(m)
        EcGM = EcGM - 4d0*rho(a,i,m)*rho(a,i,m)*(1d0-exp(-2d0*s*Diam*Diam))/Diam 
      end do
    end do
  end do

end subroutine 
