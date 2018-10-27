! test_mod_ls65_OMP.f90   2018/10/27  2017/9/19
!
! ifort -qopenmp -mkl -static -O -o test_mod_ls65_OMPstatic.ex mod_mt19937ar.f90 mod_least65_OMP.f90 test_mod_ls65_OMP.f90
! export OMP_NUM_THREADS=6
! time ./test_mod_ls65_OMP.ex < tsls65_OMP.dat > tsls65_OMP.out6
!
! ifort -qopenmp -mkl -O -o test_mod_ls65_OMP.ex mod_mt19937ar.f90 mod_least65_OMP.f90 test_mod_ls65_OMP.f90
! gfortran -fopenmp -Wall -O -o test_mod_ls65_OMP.ex mod_mt19937ar.f90 mod_least65_OMP.f90 test_mod_ls65_OMP.f90 -llapack -lblas
!
!-----
MODULE least_var
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4), parameter :: set_NSMAX=2048
  real(DP), save :: Q(set_NSMAX)
!
END MODULE least_var
!
!-----
program main
!
!$ use omp_lib
  use mod_least, only : NPR,NS,NP,NPFIT,IPFIT,PA,IPRSW,IDRSW,OBS,WSQRT,LEAST,least_alloc,least_dealloc &
  & ,copy_to_init, copy_from_init, copy_to_init_array, copy_from_init_array
  use least_var, only : Q,set_NSMAX
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4) :: MAXX,ICON
  real(DP) :: EPSR,SUMS,SIGST
  integer(4) :: i,k,kk,k_begin,k_end
  integer(4) :: imax_init
  integer(4) :: NPR_init,NS_init,NP_init,NPFIT_init,IPRSW_init,IDRSW_init,MAXX_init
  integer(4), ALLOCATABLE :: IPFIT_init(:)
  real(DP) :: EPSR_init
  real(DP), ALLOCATABLE :: PA_init(:,:)
  real(DP), ALLOCATABLE :: OBS_init(:), WSQRT_init(:)
  integer(4), ALLOCATABLE :: ICON_fin(:)
  real(DP), ALLOCATABLE :: PA_fin(:,:),SUMS_fin(:)
  integer(4), parameter :: parallel_do_unit=60, chunk=1   ! parallel_do_unit=48
!  parallel_do_unit=processing unit of  !!$OMP parallel do  e.g. 10 * OMP_NUM_THREADS
!
!$ call omp_set_dynamic(.false.)
!$ write(*,'(A)')'  omp_set_dynamic(.false.)'
!$ write(*,'(A,I10)')'  omp_get_max_threads() = ',omp_get_max_threads() ! export OMP_NUM_THREADS=6
!$ write(*,'(A,I10)')'  omp_get_num_procs() = ',omp_get_num_procs()
!$ write(*,'(A,L10)')'  omp_get_nested() = ',omp_get_nested()           ! omp_get_nested() = F
!$ write(*,'(A,L10)')'  omp_get_dynamic() = ',omp_get_dynamic()         ! omp_get_dynamic() = F
!
  call init_rand_MT    !  Initialize Mersenne-Twister pseudo-random number generator
  READ(*,*) imax_init  !  number of initial parameter set
  WRITE(*,'(A,I10)')'  imax_init = ',imax_init
  if(imax_init < 1) then
    write(*,*)' invalid imax_init'  ;  stop
  endif
!
  NPR = 6    ! write(NPR,'(A)')'  write output'
  IDRSW = 1  ! NUMERICAL DERIVATIVE
!
  READ(*,*) NS, NP, NPFIT
  WRITE(*,'(a,3I5)')'  NS, NP, NPFIT=', NS, NP, NPFIT
  if(NS > set_NSMAX) then
    write(*,*)' invalid input'  ;  stop
  endif
!
  call least_alloc  ! allocate IPFIT(1:NPFIT),OBS(1:NS),WSQRT(1:NS) etc
!
  READ(*,*) IPFIT(1:NPFIT)
  WRITE(*,'(a)')'  IPFIT'
  WRITE(*,'(5I3,2X,5I3)') IPFIT(1:NPFIT)
  READ(*,*) EPSR, MAXX, IPRSW
  WRITE(*,'(a,g15.7,I10)')'  EPSR,MAXX=',EPSR,MAXX
  WRITE(*,'(a,I10)')'  IPRSW=',IPRSW
  WRITE(*,'(a)')'  IDRSW=1   NUMERICAL DERIVATIVE'
  do i=1,NS
    READ(*,*)Q(i),OBS(i),WSQRT(i)
  end do
!
!      copy NPR     ,IDRSW     ,NS     ,NP     ,NPFIT     ,EPSR     ,MAXX     ,IPRSW
!        to NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,EPSR_init,MAXX_init,IPRSW_init
  call copy_to_init( &
   &  NPR     ,IDRSW     ,NS     ,NP     ,NPFIT     ,EPSR     ,MAXX     ,IPRSW      &
   & ,NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,EPSR_init,MAXX_init,IPRSW_init )
!
  ALLOCATE( IPFIT_init(NPFIT) )
  ALLOCATE(   OBS_init(NS)    )
  ALLOCATE( WSQRT_init(NS)    )
  ALLOCATE( PA_init(NP, parallel_do_unit) )
  ALLOCATE(  PA_fin(NP, parallel_do_unit) )
  ALLOCATE( ICON_fin(parallel_do_unit)    )
  ALLOCATE( SUMS_fin(parallel_do_unit)    )
!
!      copy IPFIT     ,OBS     ,WSQRT
!        to IPFIT_init,OBS_init,WSQRT_init
  call copy_to_init_array( NPFIT, NS       &
       & ,IPFIT     ,OBS     ,WSQRT        &
       & ,IPFIT_init,OBS_init,WSQRT_init   )
!
  call least_dealloc  ! deallocate IPFIT(1:NPFIT),OBS(1:NS),WSQRT(1:NS) etc
!
  write(*,'(A)')'  k  ICON  SUMS  PA(1:NP)'
  do k_begin = 1, imax_init, parallel_do_unit
    k_end = MIN(k_begin + parallel_do_unit -1, imax_init)
    kk = 1 + k_end - k_begin
           ! set initial PA(1:NP) and save in PA_init(1:NP, 1:kk )
    call set_PA_init(NP,PA_init,kk)
!
!$OMP parallel do default(none) &
!$OMP private(k,kk, EPSR,MAXX,SUMS,SIGST,ICON) &
!$OMP shared(k_begin,k_end, NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,IPFIT_init) &
!$OMP shared(EPSR_init,MAXX_init,IPRSW_init,OBS_init,WSQRT_init,PA_init) &
!$OMP shared(ICON_fin,SUMS_fin,PA_fin) &
!$OMP schedule(dynamic,chunk)
    do k = k_begin, k_end
      kk= 1 + k - k_begin
!          copy back NPR     ,IDRSW     ,NS     ,NP     ,NPFIT     ,EPSR     ,MAXX     ,IPRSW
!               from NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,EPSR_init,MAXX_init,IPRSW_init
      call copy_from_init( &
   &  NPR     ,IDRSW     ,NS     ,NP     ,NPFIT     ,EPSR     ,MAXX     ,IPRSW      &
   & ,NPR_init,IDRSW_init,NS_init,NP_init,NPFIT_init,EPSR_init,MAXX_init,IPRSW_init )
!
      call least_alloc    ! allocate threadprivate IPFIT(1:NPFIT),OBS(1:NS),WSQRT(1:NS) etc
!
!          copy back IPFIT     ,OBS     ,WSQRT
!               from IPFIT_init,OBS_init,WSQRT_init
      call copy_from_init_array( NPFIT, NS       &
             & ,IPFIT     ,OBS     ,WSQRT        &
             & ,IPFIT_init,OBS_init,WSQRT_init   )
!          copy back PA from PA_init
      PA(1:NP) = PA_init(1:NP, kk)
!
      CALL LEAST(EPSR,MAXX,SUMS,SIGST,ICON)
      ICON_fin(kk) = ICON         ! keep ICON in ICON_fin
      SUMS_fin(kk) = SUMS         ! keep SUMS in SUMS_fin
      PA_fin(1:NP,kk) = PA(1:NP)  ! keep PA   in PA_fin
!
      call least_dealloc    ! deallocate threadprivate IPFIT(1:NPFIT),OBS(1:NS),WSQRT(1:NS) etc
    end do
!$OMP end parallel do
!
    do k = k_begin, k_end      ! print results ICON_fin, SUMS_fin, PA_fin
      kk= 1 + k - k_begin
      write(*,'(2I5,100G15.7)') k, ICON_fin(kk), SUMS_fin(kk), PA_fin(1:NP_init,kk)
    end do
!
  end do
!
  DEALLOCATE( PA_init    )
  DEALLOCATE( ICON_fin   )
  DEALLOCATE( SUMS_fin   )
  DEALLOCATE( PA_fin     )
  DEALLOCATE( IPFIT_init )
  DEALLOCATE( OBS_init   )
  DEALLOCATE( WSQRT_init )
!
end program main
!-----
subroutine set_PA_init(NP,PA_init,kk)
!
  use mod_mt19937ar , only : genrand_res53
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
!
  integer(4),intent(in ) :: NP,kk
  real(DP),intent(out) :: PA_init(NP,kk)
  integer(4) :: i,j
!
  do i = 1, kk
    PA_init(1, i)= 1.0d0 * (1.0d0 + ( -1.0d0 + 2.0d0 * genrand_res53() ) * 0.1d0 )
    PA_init(2, i)= 1.0d0 * (1.0d0 + ( -1.0d0 + 2.0d0 * genrand_res53() ) * 0.1d0 )
    do j=3,10000,6
      if( j > NP ) exit
      PA_init(j  , i)=  0.25d0 * (1.0d0 + ( -1.0d0 + 2.0d0 * genrand_res53() ) * 0.3d0 )
      PA_init(j+1, i)= -0.90d0 + ((j-3)/6)*0.2d0 + ( -1.0d0 + 2.0d0 * genrand_res53() ) * 0.03d0
      PA_init(j+2, i)=  0.02d0 * (1.0d0 + ( -1.0d0 + 2.0d0 * genrand_res53() ) * 0.3d0 )
      PA_init(j+3, i)= -0.25d0 * (1.0d0 + ( -1.0d0 + 2.0d0 * genrand_res53() ) * 0.3d0 )
      PA_init(j+4, i)= -0.80d0 + ((j-3)/6)*0.2d0 + ( -1.0d0 + 2.0d0 * genrand_res53() ) * 0.03d0
      PA_init(j+5, i)=  0.02d0 * (1.0d0 + ( -1.0d0 + 2.0d0 * genrand_res53() ) * 0.3d0 )
    end do
  end do
!
end subroutine set_PA_init
!
!-----
subroutine FMODL
!
  use mod_least, only : PA,CALC,NS
  use least_var, only : Q
!
  IMPLICIT NONE
  integer, parameter :: DP = kind(1.0D0)
  real(DP) :: xx
  integer(4) :: I,j
!
      DO I=1,NS
        xx=PA(1)+PA(2)*Q(I)
        do j=0,18
         xx=xx+PA(3+j*3)/( ( ( Q(I)-PA(4+j*3) )/PA(5+j*3) )**2 + 1.0d0 )
        end do
        CALC(I)=xx
      end do
!
end subroutine FMODL
!
!-----
subroutine DFMODL
  IMPLICIT NONE
!
  continue
  return
!
end subroutine DFMODL
!
!-----
!  Initialize mod_mt19937ar
!             Mersenne Twister pseudo-random number generator
!    see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/mt19937ar.html
!        https://en.wikipedia.org/wiki/Mersenne_Twister
subroutine init_rand_MT
!
  use mod_mt19937ar , only : init_by_array, genrand_res53
!
  implicit none
  integer, parameter :: DP = kind(1.0D0)
!
  real(DP) :: xx
  integer(4) :: j
  integer :: init(4)=(/ 291,564,837,1110 /)    ! default of initial state of mod_mt19937ar
  integer :: length=4
!
  read(*,*) init(1:4)                          ! initial state of mod_mt19937ar
  write(*,'(/,A)')'  call init_by_array(init,length)'
  write(*,'(A,4I15)')'     init(1:4) = ',init(1:4)
  write(*,'(A,I10)') '     length = ',length
  call init_by_array(init,length)              ! set initial state of mod_mt19937ar
  write(*,'(/,A)')'  four  genrand_res53() = (0,1) with 53-bit resolution'
  do j=1,4
    xx=genrand_res53()
    write(*,'(10X,G15.7)')xx                   ! first 4 random numbers of genrand_res53()
  end do
!
end subroutine init_rand_MT
!

