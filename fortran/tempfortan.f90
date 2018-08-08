program enmtest
implicit none

integer, parameter :: dp=selected_real_kind(15,300)

! ======= Test parameter ===============================
integer :: M  ! Number of purity parameters x to try
integer :: N  ! Number of random density matrices per x value
integer :: S  ! Number of samples of each measurement to

! counters
integer :: j, k

real(kind=dp) :: x_start ! Specify purity parameter x range
real(kind=dp) :: x_end

real(kind=dp), allocatable, dimension(:) :: non_physical

complex(kind=dp), dimension(2,2) :: I, X, Y, Z

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

I(1,1)=(1, 0)
I(2,2)=(1,0)

X(1,2)=(1,0)
X(2,1)=(1,0)

Y(1,2)=-(0,1)
Y(2,1)=(0,1)

Z(1,1)=(1,0)
Z(2,2)=-(1,0)

x=transpose(x)
y=transpose(y)
z=transpose(z)


100 format(4F10.2)

write(*,*) y
end program enmtest

