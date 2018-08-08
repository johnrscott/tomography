!! Generate projector
!!
!! Generator the projectors and outcomes for the Hermitian
!! matrix A. The results are stored in two arrays. The
!! order of the projectors in the proj_A array corresponds to
!! the order of the outcomes in the outcomes_A array
!!
!! At the moment the function only works for 2x2 A.

module functions
use olis_fstdlib
implicit none
! private set precision for reals and complx
integer, parameter, private :: dp=selected_real_kind(15,300)
real(kind=dp), parameter :: pi=4*atan(1.0_dp)
!
type projector
	! each projector has a matrix and eigen val
	complex(kind=dp), dimension(2,2) :: matrix
	real(kind=dp) :: outcome
end type projector

type operator
	! each operator has a matrix
	! and projector, outcome values
	! and data 
	complex(kind=dp), dimension(2,2) :: matrix
	type (projector), dimension(2) :: proj
	real(kind=dp), dimension(:), allocatable :: data
end type operator

type basisvects
	! basisvects%x, %y, %z
	! for the three basis vectors
	complex(kind=dp), dimension(3) :: x,y,z
end type basisvects

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine constructoperators(ops, basisvects, samples)
type(operator), dimension(3) :: ops
complex(kind=dp), dimension(3,3), intent(in) :: basisvects
integer, intent(in) :: samples
integer :: j, k

do j=1,size(ops)
	write(15,*) "Op",j
	call makeops(ops(j),basisvects(j,:))
	call printvectors(ops(j)%matrix, 'Op matrix', 15)
	do k=1,2
		write(15,*) "projectors", k
		write(15,*) "outcome", ops(j)%proj(k)%outcome
		call printvectors(ops(j)%proj(k)%matrix, 'Projector matrix', 15)
	end do
	! allocate the 3 basis vectors data arrays for sampling results
	allocate(ops(j)%data(samples))
    print*, 'allocated ops data'
end do			
end subroutine constructoperators


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! make operators from the 3 basis vectors given and 
!construct projectors and outcome values
subroutine makeops(op, vect)
	complex(kind=dp), dimension(3) :: vect
	complex(kind=dp), dimension(2,2) :: paulix, pauliy, pauliz, ident
	type (operator) :: op

	paulix=0.0_dp; pauliy=0.0_dp; pauliz=0.0_dp
	paulix(1,2)=(1.0_dp, 0.0_dp);	paulix(2,1)=(1.0_dp,0.0_dp)
	pauliy(1,2)=-(0.0_dp,1.0_dp);	pauliy(2,1)=(0.0_dp,1.0_dp)
	pauliz(1,1)=(1.0_dp,0.0_dp);	pauliz(2,2)=-(1.0_dp,0.0_dp)

	op%matrix=vect(1)*paulix+vect(2)*pauliy + vect(3)*pauliz
	call makeprojector(op)
end subroutine makeops

! make projectors and outcomes for each basis vector operator
subroutine makeprojector(op)
implicit none 

type(operator) :: op
integer :: n, j
integer :: lda, ldvl, ldvr

! matrix in is a
complex(kind=dp), dimension(2,2) :: a
complex(kind=dp), allocatable, dimension(:,:) :: projectors

! left & right vectors
complex(kind=dp), allocatable, dimension(:,:) :: vl
! eigen values are w
complex(kind=dp), allocatable, dimension(:) :: w 
complex(kind=dp), dimension(2,2,2) :: proj

real(kind=dp), dimension(2) :: outcome

	a=op%matrix

	! use size of input matrix
	n=size(a,1)
	ldvl=size(a,1)

	! eigen vectors, eigen vals & temp arrays
	allocate(vl( ldvl, n ))
    allocate(w(n))
	allocate(projectors(n,n))
    w=0.0_dp
    vl=0.0_dp
    !!!! end 

    call complexeigenvects(a,w, vl) 
	
    !for the 2 eigen vectors construct projectors
	! assign matching eigen vals
	do j=1,2
	proj(j,:,:) = outerproduct(vl(j,:),conjg(vl(j,:)))
	outcome(j)=w(j)
	end do

	do j=1,2
		op%proj(j)%matrix(:,:)=proj(j,:,:)
		op%proj(j)%outcome=outcome(j)
	end do
end subroutine makeprojector


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

function rand_unitary()
implicit none
integer :: i
real(kind=dp) :: xi, phi
real(kind=dp), dimension(3) :: alpsichi
complex(kind=dp) :: glob, a, b
complex(kind=dp), dimension(3) :: imagphases
complex(kind=dp), dimension(2,2) :: rand_unitary

	!places a random number arg 
	xi=0.0_dp
	alpsichi=0.0_dp
	! 3 numbers between 0 & 2Pi
	! 1 number between 0 & 1
	call random_number(xi)
	call random_number(alpsichi)
	alpsichi=alpsichi*2.0_dp*pi

	!print*, 'alpha', char(9), 'psi', char(9),'chi',char(9),' xi'
	!print*, alpsichi, xi 

	! convert to imaginary nums
	do i=1,3
	imagphases(i)=complex(0.0_dp,alpsichi(i))
	end do

	phi=asin(sqrt(xi))

	! matrix elements
	glob=0.0_dp
	glob=exp(imagphases(1))
	!print*, glob
	a=exp(imagphases(2)) * cos(phi) 
	b=exp(imagphases(3)) * sin(phi)

	!print*, 'a', a, 'b', b
	!print*, '-conjg(b)', -conjg(b), 'conjg(a)', conjg(a)
	rand_unitary=0.0_dp
	rand_unitary=reshape((/a,b,-conjg(b), conjg(a)/), shape(rand_unitary))
	!call printvectors(rand_unitary, 'no')
	rand_unitary=glob*rand_unitary
	!write(*,*) imagphases
	!call printvectors(rand_unitary,desc='unitary')

	!print*, 'check unitary'

	!call printvectors(matmul(rand_unitary,conjg(transpose(rand_unitary))), 'what is this?')

end function rand_unitary

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function rand_density(purity)
complex(kind=dp), dimension(2,2) :: rand_density, unitary
real(kind=dp), dimension(2,2) :: temp
real(kind=dp), intent(in) :: purity

	rand_density=0.0_dp
	rand_density(1,1)=cmplx(purity, 0.0_dp); rand_density(2,2)=cmplx(1.0_dp-purity, 0.0_dp)

	unitary=rand_unitary()
	rand_density=matmul(unitary,(matmul(rand_density,conjg(transpose(unitary)))))

end function rand_density

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine simulate(dens,ops, samples)

type(operator), dimension(3) :: ops

integer :: i,j, k, samples
real(kind=dp) :: r
real(kind=dp), dimension(2) :: p
complex(kind=dp), dimension(2,2) :: dens
!complex(kind=dp), dimension(2,2,2) :: proj

do i=1,size(ops)
    do j=1, 2
        ! from oli's std lib
        p(j) = abs(complextrace(matmul(dens,ops(i)%proj(j)%matrix))) 
    end do
    do k=1,samples
        call random_number(r)
        if (r .lt. p(1)) then
            ops(i)%data(k)=ops(i)%proj(1)%outcome
            !print*, r,'outcome 1 ', 'p1',p(1), 'p2', p(2)
        elseif (r.gt.p(1)) then
            ops(i)%data(k)=ops(i)%proj(2)%outcome
            !print*, r, 'outcome 2 ', 'p1', p(1), 'p2', p(2)
        end if
        ! write the outcome in
        
        if (abs(1.0_dp-(p(1)+p(2))) > 1e-4_dp ) stop 'err'
    end do
end do
end subroutine simulate 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function linear_estimate(ops)
    type(operator), dimension(3) :: ops
    complex(kind=dp), dimension(2,2) :: paulix,pauliy,pauliz, ident, linear_estimate
    real(kind=dp), dimension(3) :: means
    integer :: i
    
    ! I
    Ident=0.0_dp
    Ident(1,1)=(1,0); Ident(2,2)=(1,0)

    ! pauli matrices 
    paulix=0.0_dp; pauliy=0.0_dp; pauliz=0.0_dp
    paulix(1,2)=(1.0_dp, 0.0_dp);	paulix(2,1)=(1.0_dp,0.0_dp)
    pauliy(1,2)=-(0.0_dp,1.0_dp);	pauliy(2,1)=(0.0_dp,1.0_dp)
    pauliz(1,1)=(1.0_dp,0.0_dp);	pauliz(2,2)=-(1.0_dp,0.0_dp)

    do i=1,size(ops)
        means(i)=sum(ops(i)%data(:))/real(size(ops(i)%data(:)), kind=dp)
    end do
    linear_estimate=(means(1)*paulix+means(2)*pauliy+means(3)*pauliz + ident)/2.0_dp
end function linear_estimate

function distance_op(a,b)
    real(kind=dp) :: distance_op
    complex(kind=dp), dimension(2,2) :: a,b, c
    real(kind=dp), dimension(2) :: sigma
    sigma=0.0_dp
    c=a-b
call complexsvd(c,sigma)
   print*,'sigma is', sigma
   distance_op=sigma(1)
end function distance_op

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module functions

