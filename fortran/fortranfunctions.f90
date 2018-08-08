!! Generate projector
!!
!! Generator the projectors and outcomes for the Hermitian
!! matrix A. The results are stored in two arrays. The
!! order of the projectors in the proj_A array corresponds to
!! the order of the outcomes in the outcomes_A array
!!
!! At the moment the function only works for 2x2 A.

module functions
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
	complex(kind=dp), dimension(2,2) :: matrix
	type (projector), dimension(2) :: proj
end type operator

type basisvects
	! basisvects%x, %y, %z
	! for the three basis vectors
	complex(kind=dp), dimension(3) :: x,y,z
end type basisvects
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
contains

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
integer, parameter   :: lwmax = 1000 

! matrix in is a
complex(kind=dp), dimension(2,2) :: a
complex(kind=dp), allocatable, dimension(:,:) :: projectors

! left & right vectors
complex(kind=dp), allocatable, dimension(:,:) :: vl, vr
! eigen values are w
complex(kind=dp), allocatable, dimension(:) :: w, work 
complex(kind=dp), dimension(2,2,2) :: proj

! temp scalars
integer :: info, lwork

!temp arrays
real(kind=dp), allocatable, dimension(:) ::  rwork
real(kind=dp), dimension(2) :: outcome

a=op%matrix

! use size of input matrix
n=size(a,1)
ldvl=size(a,1)
ldvr=size(a,1)
lda=size(a,1)

! eigen vectors, eigen vals & temp arrays
allocate(vl( ldvl, n ))
allocate(vr( ldvr, n ))
allocate(w( n ))
allocate( work( lwmax ))
allocate(rwork(2*n))

allocate(projectors(n,n))

!     .. eXECUTABLE sTATEMENTS ..
!     qUERY THE OPTIMAL WORKSPACE.
lwork = -1
call zgeev( 'V', 'N', n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
lwork = min( lwmax, int( work( 1 ) ) )

!     sOLVE EIGENPROBLEM.
call zgeev( 'v', 'n', n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )

!     cHECK FOR CONVERGENCE.
if( info.gt.0 ) then
write(*,*)'tHE ALGORITHM FAILED TO COMPUTE EIGENVALUES.'
stop
end if

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
deallocate(vl)
deallocate(vr)
deallocate(w)
deallocate(work)
deallocate(rwork)
deallocate(projectors)
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
print*, glob
a=exp(imagphases(2)) * cos(phi) 
b=exp(imagphases(3)) * sin(phi)

!print*, 'a', a, 'b', b
!print*, '-conjg(b)', -conjg(b), 'conjg(a)', conjg(a)
rand_unitary=0.0_dp
rand_unitary=reshape((/a,b,-conjg(b), conjg(a)/), shape(rand_unitary))
call printvectors(rand_unitary, 'no')
rand_unitary=glob*rand_unitary
!write(*,*) imagphases
call printvectors(rand_unitary,desc='unitary')

print*, 'check unitary'

call printvectors(matmul(rand_unitary,conjg(transpose(rand_unitary))), 'what is this?')

end function rand_unitary

function rand_density(purity)
complex(kind=dp), dimension(2,2) :: rand_density, unitary
real(kind=dp), dimension(2,2) :: temp
real(kind=dp), intent(in) :: purity

rand_density=0.0_dp
rand_density(1,1)=cmplx(purity, 0.0_dp); rand_density(2,2)=cmplx(1.0_dp-purity, 0.0_dp)

unitary=rand_unitary()
rand_density=matmul(unitary,(matmul(rand_density,conjg(transpose(unitary)))))

end function rand_density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 
subroutine printvectors(vect, desc, f)
implicit none
integer :: i, j, m, n
integer, intent(in), optional :: f 
character(len=*), intent(in), optional :: desc
complex(kind=dp), dimension(:,:), intent(in) :: vect
n=size(vect,1)
m=size(vect,1)

if (present(f)) then
        if (present(desc)) then
        write(f,*) desc
        else 
        write(f,*)
        endif 
do i=1, m
        write(f,9998) (vect (i,j), j=1,n)
end do
write(f,*)
else     
        if (present(desc)) then
        write(*,*) desc
        else 
        write(*,*)
        endif 
do i=1, m
        write(*,9998) (vect (i,j), j=1,n)
end do
write(*,*)
end if

9998 format( 11(:,1x,'(',F6.2,',',F6.2,')'))
end subroutine printvectors


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function outerproduct(a,b)
implicit none
complex(kind=dp), dimension(2,2) :: outerproduct
complex(kind=dp), dimension(:), intent(in) :: a, b
integer :: n, j ,k

! the return value is the function name
n=size(a)
do j=1,n
	do k=1, n
		outerproduct(j,k)=a(j)*b(k)
	end do
end do
end function outerproduct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module functions

