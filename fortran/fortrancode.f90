program enmtest
use functions
use olis_fstdlib
implicit none

integer, parameter :: dp=selected_real_kind(15,300)

! ======= Test parameter ===============================
integer :: M  ! Number of purity parameters x to try
integer :: N  ! Number of random density matrices per x value
integer :: S  ! Number of samples of each measurement to

! counters
integer :: j, k, ierr
! for random seed
integer, dimension(:), allocatable :: seed

real :: rand
real(kind=dp) :: x_start ! Specify purity parameter x range
real(kind=dp) :: x_end, pur, non_physical_count

real(kind=dp), allocatable, dimension(:) :: non_physical
real(kind=dp), allocatable, dimension(:,:) :: dist
real(kind=dp), dimension(2) :: outcomes_x, outcomes_y, outcomes_z

complex(kind=dp), dimension(2,2,2) :: proj_x, proj_y, proj_z
complex(kind=dp), dimension(4,2,2) :: operators
complex(kind=dp), dimension(2,2) :: I, X, Y, Z, paulix, pauliz, pauliy
complex(kind=dp), dimension(2,2) :: r_unitary, dens, dens_est 
! operator has a %matrix and
! %proj which has -> %matrix & %outcome
!
! op(1)%matrix is x matrix
! op(1)%proj(1)%matrix is the first projector matrix
! op(1)%proj(2)%outcome is the 2nd projectors outcome val
type(operator), dimension(3) :: op

! consists of x, y, z elements where each is a dim(3) vector
!type(basisvects) :: pauli, stradle_v
complex(kind=dp), dimension(3,3) :: pauli, strad_vects

call randseed(seed)

M = 2000
x_start = 0
x_end = 1
N = 500
S = 500
! simulate for each density matrix 
! ======================================================

allocate(non_physical(M))
allocate(dist(3,n))

non_physical=0.0_dp
dist=0.0_dp

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
I=0.0_dp

I(1,1)=(1,0)
I(2,2)=(1,0)

! pauli matrices 
paulix=0.0_dp; pauliy=0.0_dp; pauliz=0.0_dp
paulix(1,2)=(1.0_dp, 0.0_dp);	paulix(2,1)=(1.0_dp,0.0_dp)
pauliy(1,2)=-(0.0_dp,1.0_dp);	pauliy(2,1)=(0.0_dp,1.0_dp)
pauliz(1,1)=(1.0_dp,0.0_dp);	pauliz(2,2)=-(1.0_dp,0.0_dp)


! the pauli basis vectors, keep static here
pauli(1,:)=(/(1,0),(0,0),(0,0)/)
pauli(2,:)=(/(0,0),(1,0),(0,0)/)
pauli(3,:)=(/(0,0),(0,0),(1,0)/)

! initially make paulis, x,y,z,
! their projectors and matching outcomes
!call makeops(op(2),pauli%y)
!call makeops(op(3),pauli%z)

call constructoperators(op,pauli,s)

!end do						
!!!!!!!!!!!!!!! end of projectors 

! to reiterate, ops(1,2,3) have % matrix
! have %proj(1,2) which has %matrix or %outcome
! example
!call printvectors(op(1)%proj(1)%matrix, 'xop projector 1 matrix')

! previously called x...
pur=0.0_dp

do j=1,1! m
	non_physical_count = 0.0_dp
	do k=1,1 !n
		! purity 
		pur=x_start + j*(x_end-x_start)/real(m,kind=dp)	
       	
        ! random density matrix 
		dens=rand_density(pur) 
        
        ! for s samples 
        call simulate(dens,op,s)
        
        dens_est=linear_estimate(op)
        call printvectors(dens, 'Density matrix')
        call printvectors(dens_est, 'Estimated Density matrix')
        
        ! dist is a (3,n) array
        ! dist(1) is operator dist
        ! dist(2) is trace dist
        ! dist(3) is fidelity dist
        dist(1,k)= distance_op(dens,dens_est)
        dist(2,k)= distance_trace(dens,dens_est) 
   !     dist(3,k)= distance_fidelity(dens,dens_est)
    
    end do
    !write(*,*) op(1)%data(:)
end do

write(*,*) size(op(1)%data)


rand=0.0
r_unitary=rand_unitary()
!call printvectors(r_unitary, 'random unitary is')

!dens=rand_density(pur)
!call printvectors(dens, 'rand density')

!print*, linear_estimate(op)

deallocate(non_physical)
end program enmtest 
 


