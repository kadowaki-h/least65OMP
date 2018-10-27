!  mod_least65_OMP.f90   2017/9/19
!
!**********************************************************************
! Copyright (c) 2018, Hiroaki Kadowaki (kadowaki@tmu.ac.jp)
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without 
! modification, are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, 
! this list of conditions and the following disclaimer.
! 
! 2. Redistributions in binary form must reproduce the above copyright notice, 
! this list of conditions and the following disclaimer in the documentation and/or 
! other materials provided with the distribution.
! 
! 3. Neither the name of the copyright holder nor the names of its contributors may 
! be used to endorse or promote products derived from this software without specific 
! prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
! NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
! POSSIBILITY OF SUCH DAMAGE.
! 
!    (See https://opensource.org/licenses/BSD-3-Clause )
!
!**********************************************************************
!
!! $ gfortran -fopenmp -Wall -O -o test_mod_ls65_OMP.ex mod_mt19937ar.f90 mod_least65_OMP.f90 test_mod_ls65_OMP.f90 -llapack -lblas
!! $ export OMP_NUM_THREADS=4
!! $ ./test_mod_ls65_OMP.ex < tsls65_OMP.dat > tsls65_OMP.out
!! $ ifort -qopenmp -mkl -O -o test_mod_ls65.ex mod_mt19937ar.f90 mod_least65.f90_OMP test_mod_ls65_OMP.f90
!**********************************************************************
!      1102 LINES
!  method of least squares
!   Ver-6.5 f90 module     openMP: all module variables (static) are !!$OMP threadprivate(...)
!       kadowaki@tmu.ac.jp
!**********************************************************************
!  SUBROUTINES: LEAST MARQRT DFMDL DIAGNO DSYEVR DGEQP3 DORMQR DTRTRS
!  MODULE least_var
!**********************************************************************
!  LEAST CALCULATES PA(IPFIT(J)) (J=1-NPFIT)
!  SUMS=SUMMATION(((OBS(I)-CALC(I))/WSQRT(I))**2) (I=1-NS) = minimum
!**********************************************************************
!                   HOW TO USE
!**********************************************************************
!  NS INPUT :NUMBER OF OBSERVED DATA
!     OUTPUT:NO CHANGE
!  OBS(1-NS) INPUT :OBSERVED DATA
!               OUTPUT:NO CHANGE
!  WSQRT(1-NS) INPUT:ERROR OF OBSERVED DATA (ONE SIGMA OF NORMAL
!                    DISTRIBUTION)    IF YOU DO NOT KNOW ERROR OF
!		     OBSERVATION, GIVE APPROXIMATE VALUES.
!              OUTPUT:NO CHANGE
!  NP INPUT :NUMBER OF PARAMETERS
!     OUTPUT:NO CHANGE
!  NPFIT INPUT :NUMBER OF FITTED PARAMETERS (NPFIT >= NP)
!        OUTPUT:NO CHANGE
!  IPFIT(1-NPFIT) INPUT :INDEX OF FITTED PARAMETERS
!                 OUTPUT:NO CHANGE
!  PA(1-NP) INPUT :INITIAL VALUES OF PARAMETERS
!           OUTPUT:FITTED VALUES FOR PA(IPFIT(1-NPFIT))
!  EPSR INPUT :CONVERGENCE CRITERION VALUE
!              FOR NORMAL USE GIVE EPSR = 1.0d-5
!              IF CONVERGENCE IS TOO SLOW, GIVE 1.0d-3 > EPSR > 1.0d-5
!       OUTPUT:EPSR > 1.0d-3 --> EPSR = 1.0d-3
!              EPSR < 1.0d-5 --> EPSR = 1.0d-5
!              OTHERWISE NO CHANGE
!  MAXX INPUT :MAXIMUM ITERATION OF LEAST SQUARE LOOP (MAXX >= 1)
!       OUTPUT:NO CHANGE
!  IDRSW INPUT :SWITCH OF DERIVATIVE
!               =0 --> ANALYTIC DERIVATIVE USING DFMODL
!               =1 --> NUMERICAL DERIVATIVE USING FMODL
!        OUTPUT:NO CHANGE
!  NPR   INPUT :'UNIT' NUMBER FOR PRINT OUT (WRITE(NPR,*))
!        OUTPUT:NO CHANGE
!  IPRSW INPUT :SWITCH OF PRINT OUTPUT
!               =0 --> PRINT OUT STANDARD STATEMENTS
!               =1 --> PRINT OUT ALL STATEMENTS
!               =2 --> PRINT OUT ONLY ERROR MESSAGES
!               =3 --> PRINT OUT NOTHING
!        OUTPUT:NO CHANGE
!  CALC(1-NS) OUTPUT:CALCULATED DATA USING FITTED PA(1-NS)
!  ERPA(1-NP) OUTPUT:ERROR OF PA(1-NS)
!                    ERPA(J)=SQRT(ERMPA(J,J))
!  ERMPA(1-NP,1-NP) OUTPUT:ERROR MATRIX OF PA(1-NP)
!                   ERMPA(I,J)=<(PA(I)-TRUE PA(I))*(PA(J)-TRUE PA(J))>
!  SUMS OUTPUT:SUMMATION(((OBS(I)-CALC(I))/WSQRT(I))**2) (I=1-NS)
!              USING FITTED PA(1-NS)
!              THIS IS KAISQUARE WITH DEGREE OF FREEDOM (NS-NPFIT)
!  SIGST OUTPUT:SIGST=SQRT(SUMS/(NS-NPFIT))
!               IF YOU DO NOT KNOW ERROR OF OBSERVATION (WSQRT(1-NS)),
!               USE THIS FOR NORMARIZATION OF ERROR OF
!               PA(IPFIT(1-NPFIT)).  IN THIS CASE, ERROR OF PA(J) IS
!               SIGST*ERPA(J), AND ERROR MATRIX OF PA(J) IS
!               SIGST**2*ERMPA(I,J).
!  ICONL OUTPUT:CONDITION CODE
!               =0 NO ERROR
!               =10 ITERATION NUMBER > MAXX
!               =20 MARQUARDT NUMBER > UPPER LIMIT
!                   USUALLY FITTING IS OK
!               =30 RANK REDUCTION OF JACOBIAN MATRIX
!                   PARAMETERS ARE NOT INDEPENDENT OR MISTAKES IN
!                   DERIVATIVE
!               =40 IMPOSSIBLE TO DIAGONALIZE TRANSPOSE(AP)*(AP)
!               =50 ERROR IN NUMERICAL DERIVATIVE
!               =100 ERROR IN INPUT
!
!  WRITE FOLLOWING TWO SUBROUTINES WHICH GIVE MODEL FUNCTION (CALC(I))
!  AND ITS PARTIAL DERIVATIVE.
!    SUBROUTINE FMODL
!       INPUT: PA(1-NP)
!       OUTPUT: CALC(1-NS) MODEL FUNCTION
!    SUBROUTINE DFMODL
!       INPUT: PA(1-NP)
!       OUTPUT: AP(I,J)=PARTIAL DERIVATIVE OF CALC(I) BY PA(J)
!                      =DELTA(CALC(I))/DELTA(PA(J))
!                      I=1-NS, J=IPFIT(1-NPFIT)
!       IDRSW=0 --> THIS SUBROUTINE IS INDISPENSABLE
!       IDRSW=1 --> THIS SUBROUTINE IS NOT CALLED, BUT MUST EXIST FOR
!                   LINKAGE. "RETURN" AND "END" IS ENOUGH.
!
!***********   EXAMPLE
!MODULE least_var
!  IMPLICIT NONE
!  integer, parameter :: DP = kind(1.0D0)
!  integer(4), parameter :: set_NSMAX=2048
!  real(DP) :: Q(set_NSMAX)
!END MODULE least_var
!-----
!program main
!!$ use omp_lib
!  use mod_least, only : NPR,NS,NP,NPFIT,IPFIT,PA,IPRSW,IDRSW,OBS,WSQRT,LEAST,least_alloc,least_dealloc &
!  & ,copy_to_init, copy_from_init, copy_to_init_array, copy_from_init_array
!  use least_var, only : Q,set_NSMAX
!  IMPLICIT NONE
! .....
!!$OMP parallel do
!  call least_alloc
!  CALL LEAST(EPSR,MAXX,SUMS,SIGST,ICON)
!  call least_dealloc
!!$OMP end parallel do
! .....
!end program ts_modls63
!-----
!subroutine FMODL
!  use mod_least, only : PA,CALC,NS
!  use least_var, only : Q
!  IMPLICIT NONE
!  integer, parameter :: DP = kind(1.0D0)
!  real(DP) :: xx
!  integer(4) :: I,j
!      DO I=1,NS
!        xx=PA(1)+PA(2)*Q(I)
!        do j=0,18
!         xx=xx+PA(3+j*3)/( ( ( Q(I)-PA(4+j*3) )/PA(5+j*3) )**2 + 1.0d0 )
!        end do
!        CALC(I)=xx
!      end do
!end subroutine FMODL
!
!-----
MODULE mod_least
  IMPLICIT NONE
  private
  public :: LEAST,least_alloc,least_dealloc
  public :: copy_to_init, copy_from_init, copy_to_init_array, copy_from_init_array
  public :: NS,OBS,WSQRT,NP,NPFIT,IPFIT,PA,IDRSW,NPR,IPRSW,CALC,ERPA,ERMPA
  public :: AP
!
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4), save :: NSMAX=200  ! =UPPER LIMIT OF NS
  integer(4), save :: NPMAX=20   ! =UPPER LIMIT OF NP (NPMAX.GE.4)
  integer(4), save :: NSPMAX     ! =NSMAX+NPMAX
  integer(4), save :: NPMAX2     ! =2*NPMAX
!$OMP threadprivate(NSMAX,NPMAX,NSPMAX,NPMAX2)
!
  real(DP), save, ALLOCATABLE :: PA(:),CALC(:),WSQRT(:),AP(:,:) &
     &,OBS(:),ERMPA(:,:),ERPA(:)
  integer(4), save :: NS,NP,NPFIT,IPRSW,IDRSW,NPR=6
  integer(4), save, ALLOCATABLE :: IPFIT(:)
!$OMP threadprivate(PA,CALC,WSQRT,AP,OBS,ERMPA,ERPA)
!$OMP threadprivate(NS,NP,NPFIT,IPRSW,IDRSW,NPR)
!$OMP threadprivate(IPFIT)
!
  integer(4), save, ALLOCATABLE :: IVW(:),IPIV(:)
  real(DP), save, ALLOCATABLE :: PAULI(:),PASCAL(:),DELX(:),VP(:)
  real(DP), save, ALLOCATABLE :: VPP(:),APP(:,:),FD(:)
  real(DP), save, ALLOCATABLE :: ERMX(:,:),ERMXX(:,:)
!$OMP threadprivate(IVW,IPIV)
!$OMP threadprivate(PAULI,PASCAL,DELX,VP)
!$OMP threadprivate(VPP,APP,FD)
!$OMP threadprivate(ERMX,ERMXX)
!
  integer(4), save :: LWORK                      ! DGEQP3 + DORMQR + DTRTRS
  integer(4), save, ALLOCATABLE :: JPVT(:)       ! DGEQP3 + DORMQR + DTRTRS
  real(DP), save, ALLOCATABLE :: TAU(:),WORK(:)  ! DGEQP3 + DORMQR + DTRTRS
!$OMP threadprivate(LWORK)
!$OMP threadprivate(JPVT)
!$OMP threadprivate(TAU,WORK)
!
contains
!
!-----
subroutine copy_to_init( &
 &  NPR     ,IDRSW     ,NS     ,NP     ,NPFIT     ,EPSR     ,MAXX     ,IPRSW      &
 & ,NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,EPSR_init,MAXX_init,IPRSW_init )
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4),intent(in ) :: NPR     ,IDRSW     ,NS     ,NP     ,NPFIT     ,MAXX     ,IPRSW
  integer(4),intent(out) :: NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,MAXX_init,IPRSW_init
  real(DP),intent(in ) :: EPSR
  real(DP),intent(out) :: EPSR_init
!
  NPR_init   = NPR
  IDRSW_init = IDRSW
  NS_init    = NS
  NP_init    = NP
  NPFIT_init = NPFIT
  EPSR_init  = EPSR
  MAXX_init  = MAXX
  IPRSW_init = IPRSW
!
end subroutine copy_to_init
!
!-----
subroutine copy_from_init( &
 &  NPR     ,IDRSW     ,NS     ,NP     ,NPFIT     ,EPSR     ,MAXX     ,IPRSW      &
 & ,NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,EPSR_init,MAXX_init,IPRSW_init )
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4),intent(out) :: NPR     ,IDRSW     ,NS     ,NP     ,NPFIT     ,MAXX     ,IPRSW
  integer(4),intent(in ) :: NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,MAXX_init,IPRSW_init
  real(DP),intent(out) :: EPSR
  real(DP),intent(in ) :: EPSR_init
!
  NPR   = NPR_init
  IDRSW = IDRSW_init
  NS    = NS_init
  NP    = NP_init
  NPFIT = NPFIT_init
  EPSR  = EPSR_init
  MAXX  = MAXX_init
  IPRSW = IPRSW_init
!
end subroutine copy_from_init
!
!-----
subroutine copy_to_init_array( NPFIT, NS     &
     & ,IPFIT     ,OBS     ,WSQRT            &
     & ,IPFIT_init,OBS_init,WSQRT_init       )
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4),intent(in) :: NPFIT, NS
  integer(4),intent(out) :: IPFIT_init(NPFIT)
  integer(4),intent(in ) :: IPFIT(NPFIT)
  real(DP),intent(out) :: OBS_init(NS), WSQRT_init(NS)
  real(DP),intent(in ) :: OBS(NS), WSQRT(NS)
!
  IPFIT_init(1:NPFIT) = IPFIT(1:NPFIT)
  OBS_init(1:NS)      = OBS(1:NS)
  WSQRT_init(1:NS)    = WSQRT(1:NS)
!
end subroutine copy_to_init_array
!
!-----
subroutine copy_from_init_array( NPFIT, NS   &
     & ,IPFIT     ,OBS     ,WSQRT            &
     & ,IPFIT_init,OBS_init,WSQRT_init       )
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4),intent(in) :: NPFIT, NS
  integer(4),intent(in ) :: IPFIT_init(NPFIT)
  integer(4),intent(out) :: IPFIT(NPFIT)
  real(DP),intent(in ) :: OBS_init(NS), WSQRT_init(NS)
  real(DP),intent(out) :: OBS(NS), WSQRT(NS)
!
  IPFIT(1:NPFIT) = IPFIT_init(1:NPFIT)
  OBS(1:NS)      = OBS_init(1:NS)
  WSQRT(1:NS)    = WSQRT_init(1:NS)
!
end subroutine copy_from_init_array
!
!-----
  SUBROUTINE least_alloc
    IMPLICIT NONE
    integer(4) :: ICON                           ! DGEQP3 + DORMQR + DTRTRS
    NSMAX=NS
    NPMAX=max(4,NP)
    NSPMAX=NSMAX+NPMAX
    NPMAX2=2*NPMAX
    if( ALLOCATED(    PA) ) DEALLOCATE(    PA)
    if( ALLOCATED(  CALC) ) DEALLOCATE(  CALC)
    if( ALLOCATED( WSQRT) ) DEALLOCATE( WSQRT)
    if( ALLOCATED(    AP) ) DEALLOCATE(    AP)
    if( ALLOCATED(   OBS) ) DEALLOCATE(   OBS)
    if( ALLOCATED( ERMPA) ) DEALLOCATE( ERMPA)
    if( ALLOCATED(  ERPA) ) DEALLOCATE(  ERPA)
    if( ALLOCATED( IPFIT) ) DEALLOCATE( IPFIT)
    if( ALLOCATED(   IVW) ) DEALLOCATE(   IVW)
    if( ALLOCATED(  IPIV) ) DEALLOCATE(  IPIV)
    if( ALLOCATED( PAULI) ) DEALLOCATE( PAULI)
    if( ALLOCATED(PASCAL) ) DEALLOCATE(PASCAL)
    if( ALLOCATED(  DELX) ) DEALLOCATE(  DELX)
    if( ALLOCATED(    VP) ) DEALLOCATE(    VP)
    if( ALLOCATED(   VPP) ) DEALLOCATE(   VPP)
    if( ALLOCATED(   APP) ) DEALLOCATE(   APP)
    if( ALLOCATED(    FD) ) DEALLOCATE(    FD)
    if( ALLOCATED(  ERMX) ) DEALLOCATE(  ERMX)
    if( ALLOCATED( ERMXX) ) DEALLOCATE( ERMXX)
    if( ALLOCATED(  JPVT) ) DEALLOCATE(  JPVT)  ! DGEQP3 + DORMQR + DTRTRS
    if( ALLOCATED(   TAU) ) DEALLOCATE(   TAU)  ! DGEQP3 + DORMQR + DTRTRS
    if( ALLOCATED(  WORK) ) DEALLOCATE(  WORK)  ! DGEQP3 + DORMQR + DTRTRS
    ALLOCATE(    PA(NPMAX)        )
    ALLOCATE(  CALC(NSMAX)        )
    ALLOCATE( WSQRT(NSMAX)        )
    ALLOCATE(    AP(NSPMAX,NPMAX) )
    ALLOCATE(   OBS(NSMAX)        )
    ALLOCATE( ERMPA(NPMAX,NPMAX)  )
    ALLOCATE(  ERPA(NPMAX)        )
    ALLOCATE( IPFIT(NPMAX)        )
    ALLOCATE(   IVW(NPMAX)        )
    ALLOCATE(  IPIV(NPMAX)        )
    ALLOCATE( PAULI(NPMAX)        )
    ALLOCATE(PASCAL(NPMAX)        )
    ALLOCATE(  DELX(NPMAX)        )
    ALLOCATE(    VP(NSPMAX)       )
    ALLOCATE(   VPP(NSPMAX)       )
    ALLOCATE(   APP(NSPMAX,NPMAX) )
    ALLOCATE(    FD(NPMAX2)       )
    ALLOCATE(  ERMX(NPMAX,NPMAX)  )
    ALLOCATE( ERMXX(NPMAX,NPMAX)  )
    ALLOCATE(  JPVT(NPMAX)        )  ! DGEQP3 + DORMQR + DTRTRS
    ALLOCATE(   TAU(NPMAX)        )  ! DGEQP3 + DORMQR + DTRTRS
!
    ALLOCATE(WORK(1))                ! DGEQP3 + DORMQR + DTRTRS
      LWORK=-1
      JPVT(1:NPMAX) = 0
      APP(1:NSPMAX,1:NPMAX) = 0.0d0
      ICON=0
      CALL DGEQP3(NSPMAX,NPMAX,APP,NSPMAX,JPVT,TAU,WORK,LWORK,ICON)
!      WRITE(*,'(/,A,I10,G15.7)')'  LWORK=-1  DGEQP3 ;  ICON, WORK(1) = ',ICON,WORK(1)
      LWORK=NINT(WORK(1))
    DEALLOCATE(WORK)
    ALLOCATE(WORK(LWORK))            ! DGEQP3 + DORMQR + DTRTRS
!
  END SUBROUTINE least_alloc
!-----
  SUBROUTINE least_dealloc
    IMPLICIT NONE
    if( ALLOCATED(    PA) ) DEALLOCATE(    PA)
    if( ALLOCATED(  CALC) ) DEALLOCATE(  CALC)
    if( ALLOCATED( WSQRT) ) DEALLOCATE( WSQRT)
    if( ALLOCATED(    AP) ) DEALLOCATE(    AP)
    if( ALLOCATED(   OBS) ) DEALLOCATE(   OBS)
    if( ALLOCATED( ERMPA) ) DEALLOCATE( ERMPA)
    if( ALLOCATED(  ERPA) ) DEALLOCATE(  ERPA)
    if( ALLOCATED( IPFIT) ) DEALLOCATE( IPFIT)
    if( ALLOCATED(   IVW) ) DEALLOCATE(   IVW)
    if( ALLOCATED(  IPIV) ) DEALLOCATE(  IPIV)
    if( ALLOCATED( PAULI) ) DEALLOCATE( PAULI)
    if( ALLOCATED(PASCAL) ) DEALLOCATE(PASCAL)
    if( ALLOCATED(  DELX) ) DEALLOCATE(  DELX)
    if( ALLOCATED(    VP) ) DEALLOCATE(    VP)
    if( ALLOCATED(   VPP) ) DEALLOCATE(   VPP)
    if( ALLOCATED(   APP) ) DEALLOCATE(   APP)
    if( ALLOCATED(    FD) ) DEALLOCATE(    FD)
    if( ALLOCATED(  ERMX) ) DEALLOCATE(  ERMX)
    if( ALLOCATED( ERMXX) ) DEALLOCATE( ERMXX)
    if( ALLOCATED(  JPVT) ) DEALLOCATE(  JPVT)  ! DGEQP3 + DORMQR + DTRTRS
    if( ALLOCATED(   TAU) ) DEALLOCATE(   TAU)  ! DGEQP3 + DORMQR + DTRTRS
    if( ALLOCATED(  WORK) ) DEALLOCATE(  WORK)  ! DGEQP3 + DORMQR + DTRTRS
  END SUBROUTINE least_dealloc
!
!----
  SUBROUTINE LEAST(EPSR,MAXX,SUMS,SIGST,ICONL)
!
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    real(DP), intent(inout) :: EPSR
    real(DP), intent(out) :: SUMS,SIGST
    integer(4), intent(in) :: MAXX
    integer(4), intent(out) :: ICONL
    real(DP) :: SUMSN,WSQRTI(NSMAX),yy
    integer(4) :: J,I,ICON,NDG,NPP
!
!-----------CHECK INPUT
    IF(     NS > NSMAX   .OR. NS <= 0                    &
     & .OR. NP > NPMAX   .OR. NP <= 0                    &
     & .OR. NPFIT > NS   .OR. NPFIT > NP .OR. NPFIT <= 0 &
     & .OR. EPSR < 0.0d0 .OR. MAXX <= 0                  &
     & .OR. IPRSW >= 4   .OR. IPRSW <= -1                &
     & .OR. IDRSW >= 2   .OR. IDRSW <= -1                ) THEN
      IF(IPRSW /= 3) write(NPR,'(A)')'   INVALID NS, NP, NPFIT, EPSR, MAXX, IPRSW OR IDRSW'
      ICONL=100 ;  RETURN
    ENDIF
    DO J=1,NPFIT
      IF(IPFIT(J) <= 0 .OR. IPFIT(J) > NP) THEN
        IF(IPRSW /= 3) write(NPR,'(A,I3,A,I3)')'   INVALID IPFIT(',J,')=',IPFIT(J)
        ICONL=100 ;  RETURN
      ENDIF
      IF(J < NPFIT) THEN
        DO I=J+1,NPFIT
          IF(IPFIT(J) == IPFIT(I)) THEN
            IF(IPRSW /= 3) write(NPR,'(A,I3,A,I3,A)')'   INVALID IPFIT(',J,')=IPFIT(',I,')'
            ICONL=100 ;  RETURN
          ENDIF
        end do
      ENDIF
    end do
    DO I=1,NS
      IF(WSQRT(I) <= 0.0d0) THEN
        IF(IPRSW /= 3) WRITE(NPR,'(A,I5)')'   INVALID WSQRT(I) <= 0   I=',I
        ICONL=100 ;  RETURN
      ENDIF
    end do
!----------- INPUT OK
    IF(IPRSW == 1) write(NPR,'(A)')'   --------SUBROUTINE LEAST START---------'
    ICONL=0
    IF(EPSR > 1.0d-3) THEN
      EPSR=1.0d-3
      IF(IPRSW == 1 .OR. IPRSW == 0) write(NPR,'(A)')'   EPSR IS CHANGED TO 1.0d-3'
    ENDIF
    IF(EPSR < 1.0d-5) THEN
      EPSR=1.0d-5
      IF(IPRSW == 1 .OR. IPRSW == 0) write(NPR,'(A)')'   EPSR IS CHANGED TO 1.0d-5'
    ENDIF
    ERMPA(:,:)=0.0d0
    ERPA(:)=0.0d0
    WSQRTI(1:NS)=1.0d0/WSQRT(1:NS)
    IF(IPRSW == 1) write(NPR,'(A,G15.7)')'   CONVERGENCE CRITERION VALUE EPSR=',EPSR
!-------- LEAST SQUARE PROGRAM ---- MARQUARDT METHOD
    CALL MARQRT(EPSR,MAXX,SUMS,WSQRTI,ICON)
!
    IF(ICON /= 0 .AND. IPRSW /= 3) write(NPR,'(A,I5)')'   AT MARQRT ICON=',ICON
    ICONL=ICON
    IF(ICON == 50) RETURN
    IF(ICON /= 30) THEN
!--      CALCULATE ERROR MATRIX = ERMX
      APP(1:NS,1:NPFIT)=AP(1:NS,1:NPFIT)  ! copy
      call set_ERMX_RAMC('ERMX',ICONL,yy)
    ENDIF
!-------- DIAGNOSYS OF MATRIX AP
    if(IPRSW==0 .OR. IPRSW==1) then
      CALL DIAGNO(SUMS,WSQRTI,ICONL,ICON)
      IF(ICON /= 0) THEN
        IF(IPRSW /= 3) WRITE(NPR,'(A,I5)')'   AT DIAGNO ICON=',ICON
        ICONL=40
      ENDIF
    endif
!---------
    IF(ICONL == 30) RETURN
    DO J=1,NPFIT
    DO I=1,NPFIT
     ERMPA(IPFIT(I),IPFIT(J))=ERMX(I,J)*PASCAL(I) *PASCAL(J)
    end do
    end do
    DO I=1,NPFIT
      ERPA(IPFIT(I))=SQRT(ERMPA(IPFIT(I),IPFIT(I)))
    end do
!----------
!----CHI SQUARE
    NDG=NS-NPFIT
    IF(NS > NPFIT) THEN
      SUMSN=SUMS/NDG
!--   IF ABSOLUTE ERROR IS NOT KNOWN THE FOLLOWING ERROR SHOULD BE USED
      SIGST=SQRT(SUMSN)
    ELSE
      SUMSN=0.0d0
      SIGST=0.0d0
    ENDIF
!---------PRINT OUT RESULTS
    IF(IPRSW == 0 .OR. IPRSW == 1) THEN
      write(NPR,'(/,A,G15.7)') '   CHI SQUARE            SUMS=',SUMS
      write(NPR,'(A,I5)'   ) '   DEGREE OF FREEDOM NS-NPFIT=',NDG
      write(NPR,'(A,G15.7)') '              SUMS/(NS-NPFIT)=',SUMSN
      write(NPR,'(/,A)')'    I     PA FINAL      ERROR   ERROR*SQRT(SUMS/(NP-NPFIT))'
      do I=1,NP
        write(NPR,'(2X,I3,3G15.7)') I,PA(I),ERPA(I),SIGST*ERPA(I)
      end do
    ENDIF
!
    IF(IPRSW == 1) THEN
      NPP= MIN(10,NPFIT)
      write(NPR,'(A)')'   ERROR MATRIX OF PA'
      write(NPR,'(5X,10(I5,6X))')(IPFIT(J),J=1,NPP)
      DO I=1,NPFIT
        WRITE(NPR,'(I5,10G11.3)') IPFIT(I),(ERMPA(IPFIT(I),IPFIT(J)),J=1,NPP)
      end do
      IF(NPFIT > 10) THEN
        NPP= MIN(20,NPFIT)
        write(NPR,'(A)')'   ERROR MATRIX next PART'
        write(NPR,'(5X,10(I5,6X))')(IPFIT(J),J=11,NPP)
        DO I=1,NPFIT
          WRITE(NPR,'(I5,10G11.3)') IPFIT(I),(ERMPA(IPFIT(I),IPFIT(J)),J=11,NPP)
        end do
      ENDIF
    ENDIF
!
    IF(IPRSW == 1) THEN
      write(NPR,'(/,4X,A,9X,A,10X,A)')'I    OBS  ','CALC ','OBS-CALC'
      do I=1,NS
        write(NPR,'(I5,3G15.7)') I,OBS(I),CALC(I),OBS(I)-CALC(I)
      end do
    ENDIF
!
    IF(IPRSW == 1) write(NPR,'(A)')'   ----------LEAST END-----------'
!
  END SUBROUTINE LEAST
!
!-----
  SUBROUTINE MARQRT(EPSR,MAXX,SUMS,WSQRTI,ICONM)
!
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    real(DP), intent(in) :: EPSR,WSQRTI(*)
    real(DP), intent(out) :: SUMS
    integer(4), intent(in) :: MAXX
    integer(4), intent(out) :: ICONM
    real(DP) :: RAM,XX,RMAX,RAMC,RAMUL,EDELS,SUMSN,RSUMS
    integer(4) :: ITN,ITNS,J,I,IER,K,ISW,N,ICON
    integer(4) :: NRHS                              ! DGEQP3 + DORMQR + DTRTRS
!
    ICONM=0
    ITN=0
    ITNS=0
    RAM=0.0d0
    RAMUL=0.0d0
    IF(IDRSW == 1) THEN
!     SET INITIAL DELPA(J)=ERMX(J,1) AND PAP(J)=ERMX(J,2) FOR DFMDL
      DO J=1,NPFIT
        IF( ABS(PA(IPFIT(J))) > 1.0d-30) THEN
          ERMX(J,1)=0.0004d0*ABS(PA(IPFIT(J)))
        ELSE
          ERMX(J,1)=0.0004d0
        ENDIF
        ERMX(J,2)=PA(IPFIT(J))+4.0d0*ERMX(J,1)
      end do
    ENDIF
    IF(IPRSW == 1) write(NPR,'(A)')'   ---------SOBROUTINE MARQRT START---------'
!-----CALCULATE SUMS VP AP -------
    CALL FMODL
!
    VP(1:NS) = WSQRTI(1:NS)*(OBS(1:NS)-CALC(1:NS))
    SUMS = DOT_PRODUCT(VP(1:NS),VP(1:NS))
!-----------------
    CALL DFMDL(WSQRTI,IER)
    IF(IER == 1) THEN
      ICONM=50
      RETURN
    ENDIF
!
!-----   SET PASCAL, PAULI, AND AP
    DO J=1,NPFIT
      XX=DOT_PRODUCT(AP(1:NS,J),AP(1:NS,J))
      PASCAL(J)=SIGN(10.0d0/SQRT(XX) , PA(IPFIT(J)) )
      PAULI(J)=PA(IPFIT(J))-PASCAL(J)
      AP(1:NS,J)=PASCAL(J)*AP(1:NS,J)
    end do
!-----START ITERATION OF MARQUARDT METHOD
    IF(IPRSW == 0 .OR. IPRSW== 1) then
      WRITE(NPR,'(/,A)')'   ITN     SUMSN        MAX(ABS(DELX))     MARQUARDT NUMBER'
    endif
!-----ITN=NUMBER OF ITERATION
 40 ITN=ITN+1
!-----SOLVE VP=AP*DELX LINEAR LEAST SQUARE SOLUTION
!-----SET VP AND AP
      VP(NS+1:NS+NPFIT)=0.0d0
      AP(NS+1:NS+NPFIT,1:NPFIT)=0.0d0
      XX=SQRT(RAM)
      DO J=1,NPFIT
        AP(NS+J,J)=XX
      end do
      APP(1:NS+NPFIT,1:NPFIT)=AP(1:NS+NPFIT,1:NPFIT)  ! copy
      VPP(1:NS+NPFIT)=VP(1:NS+NPFIT)  ! copy
!
!    Solve the least squares problem min( norm2(VPP - APP X) ) for X
      K=NSPMAX
      ISW=1
      N=NS+NPFIT
!
!------LLS  HOUSEHOLDER  QR factorization
      JPVT(1:NPFIT) = 0                                     ! DGEQP3 + DORMQR + DTRTRS
      CALL DGEQP3(N,NPFIT,APP,K,JPVT,TAU,WORK,LWORK,ICON)   !  DGEQP3 + DORMQR + DTRTRS
!!!!!      CALL LAXL(APP,K,N,NPFIT,VPP,ISW,FD,IVW,ICON)     !  LAXL = ULALH + ULALB
!-----
      IF(ICON /= 0) THEN
        IF(IPRSW /= 3) then
          write(NPR,'(A,I7)')'   error AT DGEQP3 (or LAXL) ICON=',ICON
        endif
        ICONM=30 ;  RETURN
      ENDIF
!
      NRHS=1                                                            !  DGEQP3 + DORMQR + DTRTRS
      CALL DORMQR('L','T',N,NRHS,NPFIT,APP,K,TAU,VPP,K,WORK,LWORK,ICON) !  DGEQP3 + DORMQR + DTRTRS
      CALL DTRTRS('U','N','N',NPFIT,NRHS,APP,K,VPP,K,ICON)              !  DGEQP3 + DORMQR + DTRTRS
      FD(JPVT(1:NPFIT)) = VPP(1:NPFIT)                                  !  DGEQP3 + DORMQR + DTRTRS
      VPP(1:NPFIT) = FD(1:NPFIT)                                        !  DGEQP3 + DORMQR + DTRTRS
!
!------CONVERGENCE CRITERION VALUE RMAX=MAX(ABS(DELX(I=1-NPFIT)))--
      DELX(1:NPFIT)=VPP(1:NPFIT)
      RMAX=MAXVAL(ABS(DELX(1:NPFIT)))
!
      IF ( ( MOD(ITNS,5) == 0 .AND. ITNS /= 0) .OR. (RMAX > 0.1d0) .OR. (ITN == 1) ) THEN
!--       CALCULATE STANDARD MARQUARDT NUMBER = RAMC
        APP(1:NS,1:NPFIT)=AP(1:NS,1:NPFIT)   !  copy
        call set_ERMX_RAMC('RAMC',IER,RAMC)  !  calc RAMC using APP
        RAMUL=RAMC*1.0d8                     !  upper limit of RAM
      ENDIF
!------PARAMETER INDICATING NONLINEALITY RSUMS ---
!------CONVERGENCE CRITERION VALUE RMAX=MAX(ABS(DELX(I)))--
!------    RSUMS=(SUMSN-SUMS)/EDELS
!------    EDELS=EXPECTED DECREASE OF SUMS
!------    EDELS=-2*RAM*NORM(DELX)**2-NORM(AP*DELX)**2
      PA(IPFIT(1:NPFIT))=PAULI(1:NPFIT)+(1.0d0+DELX(1:NPFIT))*PASCAL(1:NPFIT)
      EDELS = -2.0d0 * RAM * DOT_PRODUCT( DELX(1:NPFIT) , DELX(1:NPFIT) )
!----------
      CALL FMODL
!----------
      DO I=1,NS
        XX= DOT_PRODUCT( AP(I,1:NPFIT) , DELX(1:NPFIT) )
        EDELS=EDELS-XX*XX
      end do
      VPP(1:NS)= ( OBS(1:NS) - CALC(1:NS) ) * WSQRTI(1:NS)
      SUMSN = DOT_PRODUCT( VPP(1:NS) , VPP(1:NS) )
      IF(EDELS < 0.0d0) THEN
        RSUMS=(SUMSN-SUMS)/EDELS
      ELSE
        RSUMS=0.0d0
      ENDIF
!------
      IF(SUMSN < SUMS) THEN
        ITNS=ITNS+1
        IF(IPRSW == 0 .OR. IPRSW == 1) write(NPR,'(I5,G20.12,2G15.7)')ITN,SUMSN,RMAX,RAM
      ENDIF
!------CHECK NEW SUMS (.LE. OR .GT.) OLD SUMS
      IF(SUMS <= SUMSN) THEN
!-------PA NOT IMPROVED. SET OLD VALUES
        PA(IPFIT(1:NPFIT))=PAULI(1:NPFIT)+PASCAL(1:NPFIT)
      ELSE
!-------SET SUMS VP AP PAULI PASCAL USING IMPROVED PA
        SUMS=SUMSN
        VP(1:NS)=VPP(1:NS)
!---
        CALL DFMDL(WSQRTI,IER)
        IF(IER == 1) THEN
          ICONM=50 ;  RETURN
        ENDIF
!---
        DO J=1,NPFIT
          XX=DOT_PRODUCT(AP(1:NS,J),AP(1:NS,J))
          PASCAL(J)=SIGN(10.0d0/SQRT(XX) , PA(IPFIT(J)) )
          PAULI(J)=PA(IPFIT(J))-PASCAL(J)
          AP(1:NS,J)=PASCAL(J)*AP(1:NS,J)
        end do
      ENDIF
!-----   CHECK CONVERGENCE
      IF(RMAX <= EPSR) THEN
!-----  CONVERGENCE   EXIT
        GOTO 250
      ELSE
!-----  NO CONVERGENCE
!-----  CHECK ITN<MAXX AND RAM<RAMUL, OTHERWISE EXIT
        IF(ITN >= MAXX) THEN
          ICONM=10 ;  GOTO 250
        ELSEIF(RAM > RAMUL) THEN
          ICONM=20 ;  GOTO 250
        ENDIF
!C-----   CONTINUE ITERATION     ADJUST RAM
        IF(RSUMS >= 0.75d0) THEN
          RAM=0.5d0*RAM
          IF(RAM < RAMC) RAM=0.0d0
        ELSEIF(RSUMS <= 0.25d0 .AND. RSUMS >= 0.0d0) THEN
          RAM=2.0d0*RAM
          IF(RAM < RAMC) RAM=RAMC
        ELSEIF(RSUMS < 0.0d0) THEN
          RAM=MIN(2.0d0-RSUMS,10.0d0)*RAM
          IF(RAM < RAMC) RAM=RAMC
        ENDIF
    GOTO 40
    ENDIF
!-----     RETURN
 250  continue
    CALC(1:NS)=-VP(1:NS)/WSQRTI(1:NS)+OBS(1:NS)
!
    IF((ITNS == 0) .AND. (IPRSW /= 3)) write(NPR,'(A)')'   I CANNOT IMPROVE PARAMETERS SOMETHING WRONG?'
    IF((IPRSW == 0) .OR. (IPRSW == 1)) write(NPR,'(A,I5)')'   ITERATION NUMBER ITN=',ITN
    IF(IPRSW == 1) THEN 
      write(NPR,'(A,G15.7)')'   STANDARD MARQUARDT NUMBER RAMC=',RAMC
      write(NPR,'(A,G15.7)')'   LAST MARQUARDT NUMBER RAM=',RAM
      write(NPR,'(A)')'   -----  MARQRT END -----'
    ENDIF
!
  END SUBROUTINE MARQRT
!
!-----
  SUBROUTINE DFMDL(WSQRTI,IER)
!
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    real(DP), intent(in) :: WSQRTI(*)
    integer(4), intent(out) :: IER
    integer(4) :: J,I,IPFJ,ISW,ISWB,K
    integer(4), parameter :: KMAX=20
    real(DP), parameter :: EPSRP=1.0d-7,EPSRC=1.0d-6,FACD=0.25d0,FACI=4.0d0
    real(DP) :: PRC(2),XX,XXX,XXXX,XXXXX,FACT,CONV,YY
!
!!      DATA EPSRP,EPSRC,FACD,FACI,KMAX /1.E-7,1.E-6, .25,  4.,20  /
!     IDRSW=0-ANALYTIC DEREVATIVE =1-NUMERICAL DERIVATIVE
!     APP(I,ISW)=DELTA(CALC(I)*WSQRTI(I))/DELTA(PA(IPFIT(J)))
!     DELTA( )=( )(PA(IPFJ)+DELPA(J))-( )(PA(IPFJ)-DELPA(J))
!     APP(I,ISWB)=PREVIOUS APP(I,ISW) BEFORE DELPA(J)=FACT*DELPA(J)
!     ISW,ISWB=1,2 OR 2,1
!     FACT=FACD OR FACI
!     PRC(ISW)=MEASURE OF ROUNDOFF ERROR OF APP(I,ISW)
!     CONV=MEASURE OF CONVERGENCE
!         =SQRT(SUM((APP(I,1)-APP(I,2))**2)/SUM(APP(I,ISW)**2))
!     PAP(J)=PREVIOUS PA(IPFFIT(J))
!     IER= ERROR CODE =0-NO ERROR =1-ERROR
    IER=0
    IF(IDRSW == 0) THEN
      CALL DFMODL
      APP(1:NS,1:NPFIT)=AP(1:NS,IPFIT(1:NPFIT))
      DO J=1,NPFIT
        AP(1:NS,J)=WSQRTI(1:NS)*APP(1:NS,J)
      end do
      RETURN
    ELSEIF(IDRSW == 1) THEN
      DO J=1,NPFIT
        IPFJ=IPFIT(J)
        XX=PA(IPFJ)
        ISW=1
        PA(IPFJ)=XX+ERMX(J,1)  ! PA(IPFJ)=XX+DELPA(J)
        CALL FMODL
        APP(1:NS,3)=CALC(1:NS)*WSQRTI(1:NS)
        PA(IPFJ)=XX-ERMX(J,1)  ! PA(IPFJ)=XX-DELPA(J)
        CALL FMODL
        XXX=0.5d0/ERMX(J,1)    !  XXX=.5/DELPA(J)
        PRC(ISW)=0.0d0
        DO I=1,NS
          APP(I,4)=CALC(I)*WSQRTI(I)
          XXXX=APP(I,3)-APP(I,4)
          APP(I,ISW)=XXXX*XXX
          XXXXX=( ABS(APP(I,3)) + ABS(APP(I,4)) )*0.5d0
          IF(XXXXX > 1.0d-60) THEN
            PRC(ISW)=MAX(PRC(ISW),ABS(XXXX/XXXXX))
          ENDIF
        end do
        IF(ABS(XX-ERMX(J,2)) <= ERMX(J,1) )THEN  ! IF(ABS(XX-PAP(J)).LE.DELPA(J))THEN
          AP(1:NS,J)=APP(1:NS,ISW)
          GOTO 180
        ENDIF
        ISW=2
        ISWB=1
        FACT=FACD
        DO K=1,KMAX
          ERMX(J,1)=FACT*ERMX(J,1)  !    DELPA(J)=FACT*DELPA(J)
          PA(IPFJ)=XX+ERMX(J,1)     !    PA(IPFJ)=XX+DELPA(J)
          CALL FMODL
          APP(1:NS,3)=CALC(1:NS)*WSQRTI(1:NS)
          PA(IPFJ)=XX-ERMX(J,1)     !    PA(IPFJ)=XX-DELPA(J)
          CALL FMODL
          XXX=0.5d0/ERMX(J,1)       !    XXX=.5/DELPA(J)
          PRC(ISW)=0.0d0
          CONV=0.0d0
          YY=0.0d0
          DO I=1,NS
            APP(I,4)=CALC(I)*WSQRTI(I)
            XXXX=APP(I,3)-APP(I,4)
            APP(I,ISW)=XXXX*XXX
            XXXXX=( ABS(APP(I,3)) + ABS(APP(I,4)) )*0.5d0
            IF(XXXXX > 1.0d-60) THEN
              PRC(ISW)=MAX(PRC(ISW),ABS(XXXX/XXXXX))
            ENDIF
            CONV=CONV+(APP(I,1)-APP(I,2))**2
            YY=YY+APP(I,ISW)**2
          end do
          IF(YY < 1.0d-60) GOTO 400
          CONV=SQRT(CONV/YY)
          IF(K == 1) THEN
            IF(PRC(1) < EPSRP .AND. PRC(2) < EPSRP) GOTO 200
          ENDIF
          IF(CONV <= EPSRC .AND. PRC(ISWB) >= EPSRP) THEN
!           CONVERGENCE NEXT J
            IF(PRC(ISW) >= EPSRP) THEN
              AP(1:NS,J)=APP(1:NS,ISW)
            ELSE
              AP(1:NS,J)=APP(1:NS,ISWB)
            ENDIF
            ERMX(J,1)=ERMX(J,1)/FACT  !    DELPA(J)=DELPA(J)/FACT
            GOTO 180
          ELSEIF(PRC(1) >= EPSRP .AND. PRC(2) >= EPSRP) THEN
!           SMALLER DELPA(J) NEXT K
            IF(ISW == 1) THEN
              ISW=2
              ISWB=1
            ELSE
              ISW=1
              ISWB=2
            ENDIF
          ELSE
!           ERROR RETURN
            IF(IPRSW /= 3) THEN
              WRITE(NPR,'(A)') '   no convergence of derivative IN DFMDL'
              WRITE(NPR,'(A)') '   decreasing DELPA(J) = delta PA(IPFIT(J))'
            ENDIF
            GOTO 300
          ENDIF
        end do
!       ERROR K=KMAX RETURN
        IF(IPRSW /= 3) THEN
          WRITE(NPR,'(A)') '   no convergence of derivative IN DFMDL'
          WRITE(NPR,'(A)') '   decreasing DELPA(J) = delta PA(IPFIT(J))'
          WRITE(NPR,'(A,I5)') '   K > KMAX=',KMAX
        ENDIF
        GOTO 300
!
 200    continue
        ERMX(J,1)=ERMX(J,1)/FACT    ! DELPA(J)=DELPA(J)/FACT
        FACT=FACI
        DO K=1,KMAX
          ERMX(J,1)=FACT*ERMX(J,1)  !   DELPA(J)=FACT*DELPA(J)
          PA(IPFJ)=XX+ERMX(J,1)     !   PA(IPFJ)=XX+DELPA(J)
          CALL FMODL
          APP(1:NS,3)=CALC(1:NS)*WSQRTI(1:NS)
          PA(IPFJ)=XX-ERMX(J,1)     !   PA(IPFJ)=XX-DELPA(J)
          CALL FMODL
          XXX=0.5d0/ERMX(J,1)       !   XXX=.5/DELPA(J)
          PRC(ISW)=0.0d0
          CONV=0.0d0
          YY=0.0d0
          DO I=1,NS
            APP(I,4)=CALC(I)*WSQRTI(I)
            XXXX=APP(I,3)-APP(I,4)
            APP(I,ISW)=XXXX*XXX
            XXXXX=( ABS(APP(I,3)) + ABS(APP(I,4)) )*0.5d0
            IF(XXXXX > 1.0d-60) THEN
              PRC(ISW)=MAX(PRC(ISW),ABS(XXXX/XXXXX))
            ENDIF
            CONV=CONV+(APP(I,1)-APP(I,2))**2
            YY=YY+APP(I,ISW)**2
          end do
          IF(YY < 1.0d-60) GOTO 400
          CONV=SQRT(CONV/YY)
          IF(CONV <= EPSRC .AND. PRC(ISW) >= EPSRP) THEN
!           CONVERGENCE NEXT J
            IF(PRC(ISWB) >= EPSRP) THEN
              AP(1:NS,J)=APP(1:NS,ISWB)
            ELSE
              AP(1:NS,J)=APP(1:NS,ISW)
            ENDIF
            GOTO 180
          ELSEIF(PRC(1) < EPSRP .AND. PRC(2) < EPSRP) THEN
!           LARGER DELPA(J) NEXT K
            IF(ISW == 1) THEN
              ISW=2
              ISWB=1
            ELSE
              ISW=1
              ISWB=2
            ENDIF
          ELSE
!           ERROR RETURN
            IF(IPRSW /= 3) THEN
              WRITE(NPR,'(A)') '   no convergence of derivative IN DFMDL'
              WRITE(NPR,'(A)') '   increasing DELPA(J) = delta PA(IPFIT(J))'
            ENDIF
            GOTO 300
          ENDIF
        end do
!       ERROR K=KMAX RETURN
        IF(IPRSW /= 3) THEN
          WRITE(NPR,'(A)') '   no convergence of derivative IN DFMDL'
          WRITE(NPR,'(A)') '   decreasing DELPA(J) = delta PA(IPFIT(J))'
          WRITE(NPR,'(A,I5)') '   K > KMAX=',KMAX
        ENDIF
        GOTO 300
!
 180    continue
        PA(IPFJ)=XX
        ERMX(J,2)=XX    ! PAP(J)=XX
      end do
      RETURN
    ENDIF
    RETURN
!
 300 continue
    IER=1
    PA(IPFJ)=XX
    IF(IPRSW /= 3) THEN
      WRITE(NPR,'(A)')'   ERROR IN DFMDL'
      WRITE(NPR,'(A)')'   K,J,IPFIT(J),PA(IPFIT(J)),DELPA(J)'
      WRITE(NPR,'(3I5,2G15.7)')K,J,IPFJ,XX,ERMX(J,1)  ! ' /3I5,2G15.7)')K,J,IPFJ,XX,DELPA(J)
      WRITE(NPR,'(A)')'   ISW,FACT,CONV,PRC(ISW),PRC(ISWB)'
      WRITE(NPR,'(I5,F7.3,3G15.7)')ISW,FACT,CONV,PRC(ISW),PRC(ISWB)
    ENDIF
    RETURN
!
 400 continue
    IER=1
    PA(IPFJ)=XX
    IF(IPRSW /= 3) THEN
      WRITE(NPR,'(A)')'   ERROR IN DFMDL'
      WRITE(NPR,'(A,I3,A)')'   ZERO DERIVATIVE - PA(',IPFJ,')'
    ENDIF
    RETURN
  END SUBROUTINE DFMDL
!
!-----
  SUBROUTINE DIAGNO(SUMS,WSQRTI,ICONL,ICOND)
!
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    real(DP), intent(in) :: SUMS
    real(DP), intent(in) :: WSQRTI(*)
    integer(4), intent(in) :: ICONL
    integer(4), intent(out) :: ICOND
    real(DP) :: CONDAP,S,SS,XX,SSUM             !  DSYEVR
!!!!!    real(DP) :: AMACH,CONDAP,S,SS,XX,SSUM  !  HSHOUS
    integer(4) :: I,J,K,ICON
    real(DP) :: AMYU(NPFIT),VORTH(NPFIT,NPFIT)        !  DSYEVR
!!!!!    real(DP) :: AMYU(NPMAX),VORTH(NPMAX,NPMAX)   !  HSHOUS
    INTEGER(4) :: il,iu,liwork,lwork,m,isuppz(2*NPFIT)  !  DSYEVR
    REAL(DP) :: abstol,vl,vu                            !  DSYEVR
    REAL(DP), Allocatable :: work(:)                    !  DSYEVR
    Integer(4), Allocatable :: iwork(:)                 !  DSYEVR
!!!!!    real(DP) :: W1(NPMAX),W2(NPMAX),W3(NPMAX),W4(NPMAX),W5(NPMAX),W6(NPMAX) ! HSHOUS
!!!!!    real(DP) :: W7(1)  ! dummy                                              ! HSHOUS
!!!!!    AMACH=1.0d-6   !      DATA AMACH/1.E-6/        ! HSHOUS
    ICOND=0
    IF(IPRSW == 1) WRITE(NPR,'(A)')'   ------SOBROUTINE DIAGNO START---------'
!
!-----DIAGNOSIS OF AP
    IF(IPRSW == 1) THEN
      write(NPR,'(A)')'   DIAGNOSIS OF JACOBIAN MATRIX (AP) FOR THE FOLLOWING PARAMETERS'
      write(NPR,'(A)')'   I=IPFIT(J)  PA(I)=PAULI(J)+X(J)*PASCAL(J)  X(J)=1'
      write(NPR,'(A)')'     I   J     PA(I)      PAULI(J)        PASCAL(J) '
      DO I=1,NP
        K=0
        DO J=1,NPFIT
          IF(I == IPFIT(J)) THEN
            WRITE(NPR,'(1X,2I5,3G15.7)')I,J,PA(I),PAULI(J),PASCAL(J)
            K=1 ;  exit
          ENDIF
        end do
        if(K==0) WRITE(NPR,'(2X,I5,5X,G15.7)')I,PA(I)
      end do
    ENDIF
    ERMXX(1:NPFIT,1:NPFIT)=MATMUL( TRANSPOSE(AP(1:NS,1:NPFIT)) , AP(1:NS,1:NPFIT) )
!-----DIAGONALIZE ERMXX=TRANSPOSE(AP)*AP
!-----                 =VORTH*DIAG(AMYU(I)**2)*TRANS(VORTH)
    IF(NPFIT == 1) THEN
      AMYU(1)=SQRT(ERMXX(1,1))
      VORTH(1,1)=1.0d0
      CONDAP=1.0d0
    ELSEIF(NPFIT >= 2) THEN
      K=NPMAX
      I=1
!         Solve the symmetric eigenvalue problem
  abstol=0.0d0  ;  lwork = -1  ;  liwork = -1        !  DSYEVR
  allocate(work(1),iwork(1))
    call DSYEVR('V','A','U',NPFIT,ERMXX,NPMAX,vl,vu,il,iu,abstol,m,AMYU,VORTH,NPFIT,isuppz,work,lwork,iwork,liwork,ICON)
    lwork = nint(work(1))  ;  liwork = iwork(1)
  deallocate(work,iwork)
  allocate(work(lwork),iwork(liwork))
    Call DSYEVR('V','A','U',NPFIT,ERMXX,NPMAX,vl,vu,il,iu,abstol,m,AMYU,VORTH,NPFIT,isuppz,work,lwork,iwork,liwork,ICON)
  deallocate(work,iwork)                             !  DSYEVR
!!!!!      CALL NSHOUS(ERMXX,K,NPFIT,NPFIT,NPFIT,AMACH,I,AMYU,VORTH,ICON,W1,W2,W3,W4,W5,W6,W7)    !  HSHOUS
      IF(ICON /= 0) THEN
        IF(IPRSW /= 3) write(NPR,'(A,I6)')'   DSYEVR  ICON=',ICON         !  DSYEVR
!!!!!        IF(IPRSW /= 3) write(NPR,'(A,I6)')'   NSHOUS  ICON=',ICON    !  HSHOUS
        ICOND=50
        GOTO 9000
      ENDIF
      AMYU(1:NPFIT)=SQRT(AMYU(1:NPFIT))
      CONDAP=AMYU(NPFIT)/AMYU(1)         !  DSYEVR
!!!!!      CONDAP=AMYU(1)/AMYU(NPFIT)    !  HSHOUS
    ENDIF
!-----PRINT OUT
    IF(IPRSW==0 .OR. IPRSW==1) write(NPR,'(A,G15.7)')'   CONDITION NUMBER CONDAP=',CONDAP
    IF(IPRSW == 1) THEN
      WRITE(NPR,'(A)')'  K SINGULAR VALUE     EIGEN VECTOR (DELTA_PA(IPFIT(J)))'
      DO K=1,NPFIT
        WRITE(NPR,'(I4,6G11.3,/,20(15X,5G11.3,/))') K,AMYU(NPFIT+1-K),(PASCAL(J)*VORTH(J,NPFIT+1-K),J=1,NPFIT)  !  DSYEVR
!!!!!        WRITE(NPR,'(I4,6G11.3,/,10(15X,5G11.3,/))') K,AMYU(K),(PASCAL(J)*VORTH(J,K),J=1,NPFIT)             !  HSHOUS
      end do
    ENDIF
!-----   STEEPEST DECRESING DIRECTION DELX=TRANS(AP)*VP/NORM
    SSUM=0.0D0
    DO J=1,NPFIT
      XX=DOT_PRODUCT(AP(1:NS,J),VP(1:NS))
      SSUM=SSUM+XX**2
      DELX(J)=XX
    end do
    SSUM=SQRT(SSUM)
    IF(SSUM > 0.0D0)THEN
      DELX(1:NPFIT)=DELX(1:NPFIT)/SSUM
    ENDIF
    IF(IPRSW == 1) then
      write(NPR,'(A)')'   - GRADIENT SUMS (DELTA_PA(IPFIT(J)))'
      write(NPR,'(20(3X,5G11.3,/))')(PASCAL(I)*DELX(I),I=1,NPFIT)
    endif
!
    IF(ICONL /= 0 .AND. ICONL /= 10 .AND. ICONL /= 20)GOTO 9000
!-----  SUMS(X+VORTH/AMYU)-SUMS(X)
!-----  SUMS(X-VORTH/AMYU)-SUMS(X)
    VPP(1:NS)=CALC(1:NS)
    IF(IPRSW==0 .OR. IPRSW==1) THEN
      write(NPR,'(/,A)')'   CHECK VALIDITY OF ERRORS'
      write(NPR,'(A)')'   IF SOME OF THE FOLLOWING VALUES ARE MUCH DIFFERENT FROM 1 (FACTOR 2),'
      write(NPR,'(A)')'   THE CALCULATED ERRORS ARE NOT RELIABLE.'
      write(NPR,'(A)')'    K  SUMS(X+-VORTH(K)/AMYU(K))-SUMS(X) =1, IF LINEAR'
    ENDIF
    DO K=1,NPFIT
      PA(IPFIT(1:NPFIT))=PAULI(1:NPFIT)+(1.0d0+VORTH(1:NPFIT,NPFIT+1-K)/AMYU(NPFIT+1-K))*PASCAL(1:NPFIT)  !  DSYEVR
!!!!!      PA(IPFIT(1:NPFIT))=PAULI(1:NPFIT)+(1.0d0+VORTH(1:NPFIT,K)/AMYU(K))*PASCAL(1:NPFIT)             !  HSHOUS
      CALL FMODL
      S=0.0d0
      DO I=1,NS
        S=S+( WSQRTI(I)*(OBS(I)-CALC(I)) )**2
      end do
      S=S-SUMS
      PA(IPFIT(1:NPFIT))=PAULI(1:NPFIT)+(1.0d0-VORTH(1:NPFIT,NPFIT+1-K)/AMYU(NPFIT+1-K))*PASCAL(1:NPFIT)  !  DSYEVR
!!!!!      PA(IPFIT(1:NPFIT))=PAULI(1:NPFIT)+(1.0d0-VORTH(1:NPFIT,K)/AMYU(K))*PASCAL(1:NPFIT)             !  HSHOUS
      CALL FMODL
      SS=0.0d0
      DO I=1,NS
        SS=SS+( WSQRTI(I)*(OBS(I)-CALC(I)) )**2
      end do
      SS=SS-SUMS
      IF(IPRSW==0 .OR. IPRSW==1) WRITE(NPR,'(I5,2G15.7)')K,S,SS
    end do
    CALC(1:NS)=VPP(1:NS)
    PA(IPFIT(1:NPFIT))=PAULI(1:NPFIT)+PASCAL(1:NPFIT)
!
 9000 IF(IPRSW == 1) WRITE(NPR,'(A)')'   ----------DIAGNO END----------'
!
  END SUBROUTINE DIAGNO
!
!-----
  subroutine set_ERMX_RAMC(EorR,ICONL,RAMC)
!
    IMPLICIT NONE
    integer, parameter :: DP = kind(1.0D0)
    character(4), intent(in) :: EorR  ! 'ERMX' (set error matrix) or 'RAMC' (set STANDARD MARQUARDT NUMBER)
    real(DP), intent(out) :: RAMC
    integer(4), intent(out) :: ICONL
    integer(4) :: J,I,ICON,K
    real(DP) :: XX
!     EorR='ERMX'  -->  set error matrix = ERMX
!         ='RAMC'  -->  set STANDARD MARQUARDT NUMBER = RAMC
!
!--      Compute the QR factorization of APP
    K=NSPMAX
    ICON=0
!
    JPVT(1:NPFIT) = 0                                       !  DGEQP3
    CALL DGEQP3(NS,NPFIT,APP,K,JPVT,TAU,WORK,LWORK,ICON)    !  DGEQP3
!!!!!    CALL ULALH(APP,K,NS,NPFIT,FD,IVW,ICON)             !  ULALH
!
    IF(ICON /= 0) THEN
      if(EorR == 'ERMX') then
        IF(IPRSW /= 3) write(NPR,'(A,I7)')'    at DGEQP3 (ULALH) in set_ERMX_RAMC ICON=',ICON
        ICONL=40
        ERMX(1:NPFIT,1:NPFIT)=0.0d0  ! set error matrix = ERMX = 0
        return
      elseif(EorR == 'RAMC')then
        IF(IPRSW /= 3) write(NPR,'(A,I7)')'    AT DGEQP3 (ULALH) IN set_ERMX_RAMC ICON=',ICON
      endif
    endif
!
!       inverting R = APP = upper triangular matrix
!                -->  APP=INV(R)
    DO I=1,NPFIT
!!!!!      APP(I,I)=1.0d0/FD(I)  !  ULALH   diagonal elements of APP=INV(R)
      APP(I,I)=1.0d0/APP(I,I)    !  DGEQP3
    end do
    IF(NPFIT >= 2) THEN
      DO I=NPFIT-1,1,-1
        DO J=NPFIT,I+1,-1
          APP(I,J)=-APP(I,I)*DOT_PRODUCT( APP(I,I+1:J) , APP(I+1:J,J) )  !  APP=INV(R) upper part
          APP(J,I)=0.0d0  !  APP=INV(R) lower part
        end do
      end do
    ENDIF
!
    if(EorR == 'ERMX') then
!--     CALCULATE ERROR MATRIX = ERMXX  USING  R
!--     ERMXX=INV(R)*TRANSPOSE(INV(R))=APP * TRANSPOSE(APP)=ERROR MATRIX (not permuted)
      ERMXX(1:NPFIT,1:NPFIT)=MATMUL( APP(1:NPFIT,1:NPFIT), TRANSPOSE(APP(1:NPFIT,1:NPFIT)) )
!
!--     CALCULATE ERROR MATRIX  ERMX
!--     PERFORM PERMUTATION ERMX=INV(TRANSOSE(AP)*AP)=ERROR MATRIX OF X
!--     SET IPIV (permuation) using  PIVOTING info
      IPIV(JPVT(1:NPFIT))=(/ (i,i=1,NPFIT) /)    !  DGEQP3
!!!!!      IPIV(1:NPFIT)= (/ (I, I=1,NPFIT) /)   !  ULALH
!!!!!      DO I=NPFIT,1,-1
!!!!!        IF(I /= IVW(I)) THEN
!!!!!          J=IPIV(I)
!!!!!          IPIV(I)=IPIV(IVW(I))
!!!!!          IPIV(IVW(I))=J
!!!!!        ENDIF
!!!!!      end do
      ERMX(1:NPFIT,1:NPFIT)=ERMXX(IPIV(1:NPFIT),IPIV(1:NPFIT))
!
    elseif(EorR == 'RAMC')then
!--   CALCULATE RAMC using R     RAMC = 1.0d0/trace( INV(R) TRANSPOSE(INV(R)) )
!--   XX=trace( INV(R) TRANSPOSE(INV(R)) )
!       =trace( APP TRANSPOSE(APP) )  APP=INV(R)=UPPER TRIANGULAR MATRIX
      XX=0.0D0
      DO I=1,NPFIT
        XX=XX + DOT_PRODUCT( APP(I,I:NPFIT) , APP(I,I:NPFIT) )
      end do
      RAMC=ABS(1.0d0/XX)
    endif
!
  end subroutine set_ERMX_RAMC
!
END MODULE mod_least
!

