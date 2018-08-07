program enmtest
use tools
implicit none

integer, parameter :: dp=selected_real_kind(15,300)

! ======= Test parameter ===============================
integer :: M  ! Number of purity parameters x to try
integer :: N  ! Number of random density matrices per x value
integer :: S  ! Number of samples of each measurement to

! counters
integer :: i, j, k

real(kind=dp) :: x_start ! Specify purity parameter x range
real(kind=dp) :: x_end

real(kind=dp), allocatable, dimension(:) :: non_physical

!real(kind=dp), allocatable, dimension(:,:) ::  

M = 20
x_start = 0
x_end = 1
N = 500
S = 500
! simulate for each density matrix 
! ======================================================

allocate(non_physical(0:M-1))
non_physical=0.0_dp

do i=0,M-1
  non_physical(i) = i
end do
!
do i=0,M-1
  print*, "loop", i
end do




!
n=thing(1)
call doubleit(n)


print*,"This is the function", n
write(*,*) non_physical


deallocate(non_physical)
end program enmtest 
 


