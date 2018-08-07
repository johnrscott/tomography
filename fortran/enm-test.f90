program enmtest
use tools
implicit none

integer, parameter :: dp=selected_real_kind(15,300)

! ======= Test parameter ===============================
integer :: M  ! Number of purity parameters x to try
real(kind=dp) :: x_start ! Specify purity parameter x range
real(kind=dp) :: x_end
integer :: N  ! Number of random density matrices per x value
integer :: S  ! Number of samples of each measurement to
real(kind=dp), dimension(1, 2000) :: non_physical

M = 2000
x_start = 0
x_end = 1
N = 500
S = 500
! simulate for each density matrix 
! ======================================================


n=thing(1)
call doubleit(n)


print*,"This is the function", n



end program enmtest 
 


