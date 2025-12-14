module mymod
  implicit none

contains

function linspace(a, b, n) result(x)
    implicit none
    real(kind=8), intent(in) :: a, b
    integer, intent(in)      :: n
    real(kind=8) :: x(n)

    integer :: i
    real(kind=8) :: step

    if (n < 1) then
        ERROR STOP 'linspace: n must be at least 1'
    endif

    if (n == 1) then
        x(1) = a
    else
        step = (b - a) / real(n - 1, kind=8)
        do i = 1, n
            x(i) = a + step * real(i - 1, kind=8)
        enddo
    endif
end function linspace

function return_fn(aprime,a,z,r,w) result(res)
    implicit none
    ! Declare input arguments
    real(8), intent(in) :: aprime, a, z
    real(8), intent(in) :: r, w
    ! Declare function result
    real(8) :: res
    ! Declare local variables
    real(8) :: consumption

    consumption = (1.0d0+r)*a + w*z - aprime
    if (consumption > 0.0d0) then
        res = log(consumption)
    else
        res = -1.0d10
    endif
end function return_fn

end module mymod

program main
  use mymod,   only: return_fn, linspace
  use omp_lib, only: omp_get_wtime
  implicit none

  integer, parameter :: n_a = 1000, n_z = 100
  real(8), parameter :: r = 0.04d0, w = 1.0d0

  integer :: iap, ia, iz, alloc_status
  real(8) :: a_grid(n_a), z_grid(n_z)
  real(8), allocatable :: payoff(:,:,:), payoff2(:,:,:)
  real(8) :: t_start, t_end, err
  integer :: par_fortran

  ! Toggle OpenMP region on/off (set to 0 to force serial)
  par_fortran = 1

  ! Build vectors for a and z
  a_grid = linspace(0.001d0, 50.0d0, n_a)
  z_grid = linspace(0.5d0, 1.5d0, n_z)

  ! Allocate 3D arrays to hold results
  allocate(payoff(n_a, n_a, n_z), payoff2(n_a, n_a, n_z), stat=alloc_status)
  if (alloc_status /= 0) then
    print *, "Error allocating memory for payoff arrays"
    stop
  endif

  ! -------------------------
  ! Serial computation timing
  ! -------------------------
  t_start = omp_get_wtime()

  do iz = 1, n_z
    do ia = 1, n_a
      do iap = 1, n_a
        payoff(iap, ia, iz) = return_fn(a_grid(iap), a_grid(ia), z_grid(iz), r, w)
      enddo
    enddo
  enddo

  t_end = omp_get_wtime()
  print *, "Serial wall time (s): ", t_end - t_start

  ! -------------------------
  ! OpenMP computation timing
  ! -------------------------
  t_start = omp_get_wtime()

  !$omp parallel default(shared) private(iz, ia, iap) if (par_fortran == 1)
  !$omp do collapse(3) schedule(static)
  do iz = 1, n_z
    do ia = 1, n_a
      do iap = 1, n_a
        payoff2(iap, ia, iz) = return_fn(a_grid(iap), a_grid(ia), z_grid(iz), r, w)
      enddo
    enddo
  enddo
  !$omp end do
  !$omp end parallel

  t_end = omp_get_wtime()
  print *, "OpenMP wall time (s): ", t_end - t_start

  ! Check results are correct
  err = maxval(abs(payoff - payoff2))

  if (err > 1.0d-10) then
    print *, "Error: results do not match! max error = ", err
  else
    print *, "Results match! max error = ", err
  endif

  deallocate(payoff, payoff2)

end program main


