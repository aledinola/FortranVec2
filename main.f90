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
use mymod, only: return_fn, linspace
implicit none
 
integer, parameter :: n_a = 1000, n_z = 50
real(8), parameter :: r = 0.04d0, w = 1.0d0 
integer :: iap,ia,iz, alloc_status
real(8) :: a_grid(n_a), z_grid(n_z)
real(8), allocatable :: payoff(:,:,:)
real(8) :: t_start, t_end

! Build vectors for a and z
a_grid = linspace(0.001d0, 50.0d0, n_a)
z_grid = linspace(0.5d0, 1.5d0, n_z)

! Allocate 3D array to hold results
allocate(payoff(n_a, n_a, n_z), stat=alloc_status)
if (alloc_status /= 0) then
    print *, "Error allocating memory for payoff array"
    stop
endif 

call cpu_time(t_start)
! Construct 3D array to hold results using nested loops
do iz = 1, n_z
  do ia = 1, n_a
    do iap = 1, n_a
        payoff(iap,ia,iz) = return_fn(a_grid(iap),a_grid(ia),z_grid(iz),r,w)
    enddo
  enddo
enddo

call cpu_time(t_end)
print *, "CPU time (s): ", t_end - t_start

!do iap=1,5
!  do ia=1,5
!    do iz=1,5
!      print *, "payoff(",iap,",",ia,",",iz,") = ", payoff(iap,ia,iz)
!    enddo
!  enddo
!enddo

deallocate(payoff)

end program main

