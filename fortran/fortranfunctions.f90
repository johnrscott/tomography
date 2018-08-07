!! Generate projector
!!
!! Generator the projectors and outcomes for the Hermitian
!! matrix A. The results are stored in two arrays. The
!! order of the projectors in the proj_A array corresponds to
!! the order of the outcomes in the outcomes_A array
!!
!! At the moment the function only works for 2x2 A.
!!

module functions
implicit none

integer, parameter :: dp1=selected_real_kind(15,300)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine makeprojectors(a, proj, outcome)
implicit none 

integer :: n, j, k
integer :: lda, ldvl, ldvr
integer, parameter   :: lwmax = 1000 

complex(kind=dp1), dimension(:,:), intent(in) :: a
complex(kind=dp1), allocatable, dimension(:,:) :: projectors

complex(kind=dp1), allocatable, dimension(:,:) :: vl, vr
complex(kind=dp1), allocatable, dimension(:) :: w, work 
complex(kind=dp1), dimension(:,:,:) :: proj

! temp scalars
integer :: info, lwork

!temp arrays
real(kind=dp1), allocatable, dimension(:) ::  rwork
real(kind=dp1), dimension(2) :: outcome

n=size(a,1)
ldvl=size(a,1)
ldvr=size(a,1)
lda=size(a,1)

allocate(vl( ldvl, n ))
allocate(vr( ldvr, n ))
allocate(w( n ))
allocate( work( lwmax ))
allocate(rwork(2*n))

allocate(projectors(n,n))
!     .. eXECUTABLE sTATEMENTS ..
write(*,*)'zgeev eXAMPLE pROGRAM rESULTS'

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

!     pRINT EIGENVALUES.
!write(*,*)'eIGENVALUES', 1, n, w, 1 
!call printvectors((w))


!     pRINT LEFT EIGENVECTORS.
!write(*,*) 'lEFT EIGENVECTORS'!, n, n, real(vl), ldvl 
!call printvectors(vl)

!     pRINT RIGHT EIGENVECTORS.
!write(*,*) 'rIGHT EIGENVECTORS'!, n, n, real(vr), ldvr 
!call printvectors(vr)

!write(*,*) 'print vector all!'
!write(*,*) vl, char(10)

!write(*,*) 'print vectors 0,1'
!write(*,*) vl(1,:), char(10)

!print*, 'print vectors 2,3'
!write(*,*) vl(2,:), char(10)

!print*, ' vectors * vectors.T'
!write(*,*) (vl(1,:)*vl(1,:))

! projectors=matmul(vl,transpose(vl))

do j=1,2
write(*,*)w(j)
write(*,*) vl(j,:)
!proj(j,:,:)=spread(vl(j,:),dim=2,ncopies=2)*spread(vl(:,j),dim=1,ncopies=2)
proj(j,:,:) = outerproduct(vl(j,:),conjg(vl(j,:)))
!write(*,*) proj(j,:,:)
end do
!write(*,*) vl(2,:)

! make 
!proj(1,:,:)=(:,:)
outcome(1)=w(1)

!proj(2,:,:)=projectors(:,:)
outcome(2)=w(2)

end subroutine makeprojectors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11 
subroutine printvectors(vect)
integer :: i, j, m, n
complex(kind=dp1), dimension(:,:) :: vect
n=size(vect,1)
m=size(vect,1)

write(*,*)
write(*,*) 'description'
do i=1, m
        write(*,9998) (vect (i,j), j=1,n)
end do
write(*,*)

9998 format( 11(:,1x,'(',F6.2,',',F6.2,')'))
end subroutine printvectors

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine make_projector() 
print*, "This is make_projector"
end subroutine make_projector

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function outerproduct(a,b)
complex(kind=dp1), dimension(2,2) :: outerprod, outerproduct
complex(kind=dp1), dimension(:), intent(in) :: a, b
integer :: n, j ,k

n=size(a)
!allocate(outerproduct(n,n))

do j=1,n
	do k=1, n
		outerproduct(j,k)=a(j)*b(k)
	end do
end do
end function outerproduct
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module functions

