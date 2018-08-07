program enmtest
use functions
implicit none

integer, parameter :: dp=selected_real_kind(15,300)

! ======= Test parameter ===============================
integer :: M  ! Number of purity parameters x to try
integer :: N  ! Number of random density matrices per x value
integer :: S  ! Number of samples of each measurement to

! counters
integer :: j, k, ierr

real(kind=dp) :: x_start ! Specify purity parameter x range
real(kind=dp) :: x_end

real(kind=dp), allocatable, dimension(:) :: non_physical
real(kind=dp), dimension(2) :: outcomes_x, outcomes_y, outcomes_z

complex(kind=dp), dimension(2,2,2) :: proj_x, proj_y, proj_z
!real(kind=dp), allocatable, dimension(:,:) ::  
complex(kind=dp), dimension(4,2,2) :: operators
complex(kind=dp), dimension(2,2) :: I, X, Y, Z

M = 20
x_start = 0
x_end = 1
N = 500
S = 500
! simulate for each density matrix 
! ======================================================

allocate(non_physical(0:M-1))
non_physical=0.0_dp

!! Get an output file ready
open(unit=15, file='linear_test_1_fortran.dat', status='replace', iostat=ierr)
        if (ierr/=0) stop 'Error opening outfile'

write(15,*) "Distances between estimated and original density matrices"
write(15,*) "using various distances. Go to the end of the file for the running time."

  !! Write the simulation parameters
write(15,*) "Number of purity values tried = ",  M
write(15,*) "Number of density matrices per purity parameter = ", N
write(15,*) "Total number of measurements for each of X, Y and Z = ", S 
write(15,*)

write(15,*) "PURITY, \tOPERATOR, \tTRACE, \t\tFIDELITY, \tNON PHYSICAL"

!! Set precision
100 format(4F10.2)


!! Preliminaries: define measurement operators

! I
I(1,1)=(1, 0)
I(2,2)=(1,0)

X(1,2)=(1,0)
X(2,1)=(1,0)

Y(1,2)=-(0,1)
Y(2,1)=(0,1)

Z(1,1)=(1,0)
Z(2,2)=-(1,0)


!! I
!operators(1,1,1)=(1, 0)
!operators(1,2,2)=(1,0)
!! X
!operators(2,1,2)=(1,0)
!operators(2,2,1)=(1,0)
!! Y
!operators(3,1,2)=-(0,1)
!operators(3,2,1)=(0,1)
!! Z
!operators(4,1,1)=(1,0)
!operators(4,2,2)=-(1,0)

! transpose as column majored
x=transpose(x)
y=transpose(y)
z=transpose(z)

!allocate( (proj1 proj2) (rows) (cols))

proj_x=0.0_dp
proj_y=0.0_dp
proj_z=0.0_dp

outcomes_x=0.0_dp
outcomes_y=0.0_dp
outcomes_z=0.0_dp


print*,"This is the function", n
write(*,*) non_physical
call make_projector()

! X
call printvectors(x)
call checklapack(x, proj_x, outcomes_x)

write(*,*) "projector 1 for x value", outcomes_X(1)
call printvectors(proj_x(1,:,:))
print*, " projector 2 for x value", outcomes_x(2)
call printvectors(proj_x(2,:,:))


!! y
call printvectors(y)
call checklapack(y, proj_y, outcomes_y)

write(*,*) "projector 1 for y", outcomes_y(1)
call printvectors(proj_y(1,:,:))
print*, " projector 2 for y", outcomes_y(2)
call printvectors(proj_y(2,:,:))

!! z
call printvectors(z)
call checklapack(z, proj_z, outcomes_z)

write(*,*) "projector 1 for z", outcomes_z(1)
call printvectors(proj_z(1,:,:))
print*, " projector 2 for z", outcomes_z(2)
call printvectors(proj_z(2,:,:))

!!!!!!!!!!!!!!! end of projectors 



deallocate(non_physical)
end program enmtest 
 


