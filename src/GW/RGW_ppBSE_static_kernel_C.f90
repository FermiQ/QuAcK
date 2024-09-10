subroutine RGW_ppBSE_static_kernel_C(ispin,eta,nOrb,nC,nO,nV,nR,nS,nVV,lambda,ERI,Om,rho,KC)

! Compute the VVVV block of the static screening W for the pp-BSE

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)

! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: chi
  double precision              :: eps
  double precision              :: tmp_ab, lambda4, eta2
  integer                       :: a,b,c,d,ab,cd,m
  integer                       :: a0, aa

  double precision, allocatable :: Om_tmp(:)
  double precision, allocatable :: tmp_m(:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)

! Output variables

  double precision,intent(out)  :: KC(nVV,nVV)

!---------------!
! Singlet block !
!---------------!

  if(ispin == 1) then

    a0 = nOrb - nR - nO
    lambda4 = 4.d0 * lambda
    eta2 = eta * eta

    allocate(tmp_m(nOrb,nOrb,nS))
    allocate(tmp(nOrb,nOrb,nOrb,nOrb))

    !$OMP PARALLEL DEFAULT(NONE)         &
    !$OMP          PRIVATE(m, c, a, eps) &
    !$OMP          SHARED(nS, nOrb, eta2, Om, rho, tmp_m)
    !$OMP DO
    do m = 1, nS
      eps = Om(m) / (Om(m)*Om(m) + eta2)
      do c = 1, nOrb
        do a = 1, nOrb
          tmp_m(a,c,m) = eps * rho(a,c,m)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "T", nOrb*nOrb, nOrb*nOrb, nS, 1.d0,       &
               tmp_m(1,1,1), nOrb*nOrb, rho(1,1,1), nOrb*nOrb, &
               0.d0, tmp(1,1,1,1), nOrb*nOrb)

    deallocate(tmp_m)

    !$OMP PARALLEL DEFAULT(NONE)                           &
    !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, tmp_ab) &
    !$OMP          SHARED(nO, nOrb, nR, nS, a0, lambda4, tmp, KC)
    !$OMP DO
    do a = nO+1, nOrb-nR
      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
      do b = a, nOrb-nR
        ab = aa + b

        tmp_ab = lambda4
        if(a .eq. b) then
          tmp_ab = 0.7071067811865475d0 * lambda4
        endif

        cd = 0
        do c = nO+1, nOrb-nR
          do d = c, nOrb-nR
            cd = cd + 1

            KC(ab,cd) = -tmp_ab * (tmp(a,c,b,d) + tmp(a,d,b,c))
            if(c .eq. d) then
              KC(ab,cd) = 0.7071067811865475d0 * KC(ab,cd)
            endif
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmp)


!    do a=nO+1,nOrb-nR
!      do b=a,nOrb-nR
!        ab = ab + 1
!        cd = 0
!        do c=nO+1,nOrb-nR
!          do d=c,nOrb-nR
!            cd = cd + 1
!
!              chi = 0d0
!              do m=1,nS
!                eps = Om(m)**2 + eta**2
!                chi = chi - rho(a,c,m)*rho(b,d,m)*Om(m)/eps &
!                          - rho(a,d,m)*rho(b,c,m)*Om(m)/eps
!              end do


!    --- --- ---
!    OpenMP implementation
!    --- --- ---
!
!    a0 = nOrb - nR - nO
!    lambda4 = 4.d0 * lambda
!    eta2 = eta * eta
!
!    allocate(Om_tmp(nS))
!
!    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(m) SHARED(nS, eta2, Om, Om_tmp)
!    !$OMP DO
!    do m = 1, nS
!      Om_tmp(m) = Om(m) / (Om(m)*Om(m) + eta2)
!    enddo
!    !$OMP END DO
!    !$OMP END PARALLEL
!
!    !$OMP PARALLEL DEFAULT(NONE)                              &
!    !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, m, tmp_ab) &
!    !$OMP          SHARED(nO, nOrb, nR, nS, a0, lambda4, Om_tmp, rho, KC)
!    !$OMP DO
!    do a = nO+1, nOrb-nR
!      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
!      do b = a, nOrb-nR
!        ab = aa + b
!
!        tmp_ab = lambda4
!        if(a .eq. b) then
!          tmp_ab = 0.7071067811865475d0 * lambda4
!        endif
!
!        cd = 0
!        do c = nO+1, nOrb-nR
!          do d = c, nOrb-nR
!            cd = cd + 1
!
!            KC(ab,cd) = 0d0
!            do m = 1, nS
!              KC(ab,cd) = KC(ab,cd) - rho(a,c,m) * rho(b,d,m) * Om_tmp(m) &
!                                    - rho(a,d,m) * rho(b,c,m) * Om_tmp(m)
!            end do
!
!            KC(ab,cd) = tmp_ab * KC(ab,cd)
!            if(c .eq. d) then
!              KC(ab,cd) = 0.7071067811865475d0 * KC(ab,cd)
!            endif
!          enddo
!        enddo
!      enddo
!    enddo
!    !$OMP END DO
!    !$OMP END PARALLEL
!
!    deallocate(Om_tmp)
!    --- --- ---


!    --- --- ---
!    Naive implementation
!    --- --- ---
!
!    ab = 0
!    do a=nO+1,nOrb-nR
!      do b=a,nOrb-nR
!        ab = ab + 1
!        cd = 0
!        do c=nO+1,nOrb-nR
!          do d=c,nOrb-nR
!            cd = cd + 1
!
!              chi = 0d0
!              do m=1,nS
!                eps = Om(m)**2 + eta**2
!                chi = chi - rho(a,c,m)*rho(b,d,m)*Om(m)/eps &
!                          - rho(a,d,m)*rho(b,c,m)*Om(m)/eps
!              end do
!
!              KC(ab,cd) = 4d0*lambda*chi/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
!
!          end do
!        end do
!      end do
!    end do
!    --- --- ---

  end if

!---------------!
! Triplet block !
!---------------!

  if(ispin == 2) then

    a0 = nOrb - nR - nO - 1
    lambda4 = 4.d0 * lambda
    eta2 = eta * eta

    allocate(tmp_m(nOrb,nOrb,nS))
    allocate(tmp(nOrb,nOrb,nOrb,nOrb))

    !$OMP PARALLEL DEFAULT(NONE)         &
    !$OMP          PRIVATE(m, c, a, eps) &
    !$OMP          SHARED(nS, nOrb, eta2, Om, rho, tmp_m)
    !$OMP DO
    do m = 1, nS
      eps = Om(m) / (Om(m)*Om(m) + eta2)
      do c = 1, nOrb
        do a = 1, nOrb
          tmp_m(a,c,m) = eps * rho(a,c,m)
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    call dgemm("N", "T", nOrb*nOrb, nOrb*nOrb, nS, 1.d0,       &
               tmp_m(1,1,1), nOrb*nOrb, rho(1,1,1), nOrb*nOrb, &
               0.d0, tmp(1,1,1,1), nOrb*nOrb)

    deallocate(tmp_m)

    !$OMP PARALLEL DEFAULT(NONE)                   &
    !$OMP          PRIVATE(a, b, aa, ab, c, d, cd) &
    !$OMP          SHARED(nO, nOrb, nR, nS, a0, lambda4, tmp, KC)
    !$OMP DO
    do a = nO+1, nOrb-nR
      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO - 1
      do b = a+1, nOrb-nR
        ab = aa + b

        cd = 0
        do c = nO+1, nOrb-nR
          do d = c+1, nOrb-nR
            cd = cd + 1

            KC(ab,cd) = lambda4 * (-tmp(a,c,b,d) + tmp(a,d,b,c))
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL

    deallocate(tmp)


!    --- --- ---
!    OpenMP implementation
!    --- --- ---
!
!    a0 = nOrb - nR - nO - 1
!    lambda4 = 4.d0 * lambda
!    eta2 = eta * eta
!
!    allocate(Om_tmp(nS))
!
!    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(m) SHARED(nS, eta2, Om, Om_tmp)
!    !$OMP DO
!    do m = 1, nS
!      Om_tmp(m) = Om(m) / (Om(m)*Om(m) + eta2)
!    enddo
!    !$OMP END DO
!    !$OMP END PARALLEL
!
!    !$OMP PARALLEL DEFAULT(NONE)                      &
!    !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, m) &
!    !$OMP          SHARED(nO, nOrb, nR, nS, a0, lambda4, Om_tmp, rho, KC)
!    !$OMP DO
!    do a = nO+1, nOrb-nR
!      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO - 1
!      do b = a+1, nOrb-nR
!        ab = aa + b
!
!        cd = 0
!        do c = nO+1, nOrb-nR
!          do d = c+1, nOrb-nR
!            cd = cd + 1
!
!            KC(ab,cd) = 0d0
!            do m = 1, nS
!              KC(ab,cd) = KC(ab,cd) - rho(a,c,m) * rho(b,d,m) * Om_tmp(m) &
!                                    + rho(a,d,m) * rho(b,c,m) * Om_tmp(m)
!            end do
!
!            KC(ab,cd) = lambda4 * KC(ab,cd)
!          enddo
!        enddo
!      enddo
!    enddo
!    !$OMP END DO
!    !$OMP END PARALLEL
!
!    deallocate(Om_tmp)
!    --- --- ---


!    --- --- ---
!    Naive implementation
!    --- --- ---
!
!    ab = 0
!    do a=nO+1,nOrb-nR
!      do b=a+1,nOrb-nR
!        ab = ab + 1
!        cd = 0
!        do c=nO+1,nOrb-nR
!          do d=c+1,nOrb-nR
!            cd = cd + 1
!
!            chi = 0d0
!            do m=1,nS
!              eps = Om(m)**2 + eta**2
!              chi = chi - rho(a,c,m)*rho(b,d,m)*Om(m)/eps &
!                        + rho(a,d,m)*rho(b,c,m)*Om(m)/eps
!            end do
!           
!            KC(ab,cd) = 4d0*lambda*chi
!
!          end do
!        end do
!      end do
!    end do
!    --- --- ---

  end if

end subroutine 
