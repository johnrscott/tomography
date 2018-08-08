module olis_fstdlib
implicit none

integer, parameter, private :: dp=selected_real_kind(15,300)
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine randseed(seed)
!!!!! random seed
integer :: values(1:8), seedsize
integer, dimension(:), allocatable :: seed

call date_and_time(values=values)
call random_seed(size=seedsize)

allocate(seed(1:seedsize))

seed(:) = values(8)

call random_seed(put=seed)
end subroutine randseed

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

function complextrace(a)
        complex(kind=dp), dimension(:,:) :: a
        complex(kind=dp) :: complextrace 
        integer :: i

        complextrace=0.0_dp
        do i=1, size(a,1)
                complextrace=complextrace+a(i,i)
        end do
end function complextrace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine complexeigenvects(a, w, vl)
! matrix a in, eigenvals, eigenvects out
implicit none 

integer :: n, j
integer :: lda, ldvl, ldvr
integer, parameter   :: lwmax = 1000 

! matrix in is a
complex(kind=dp), dimension(:,:), intent(in) :: a

! left & right vectors
complex(kind=dp), allocatable, dimension(:,:) :: vl, vr
! eigen values are w
complex(kind=dp), allocatable, dimension(:) :: w, work 

! temp scalars
integer :: info, lwork

!temp arrays
real(kind=dp), allocatable, dimension(:) ::  rwork

! use size of input matrix
n=size(a,1)
ldvl=size(a,1)
ldvr=size(a,1)
lda=size(a,1)

! eigen vectors, eigen vals & temp arrays
!allocate(vl( ldvl, n ))
allocate(vr( ldvr, n ))
!allocate(w( n ))
allocate( work( lwmax ))
allocate(rwork(2*n))

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
      

deallocate(vr)
deallocate(work)
deallocate(rwork)

end subroutine complexeigenvects

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine complexsvd(a, s)
! matrix a in, eigenvals, eigenvects out
implicit none 

integer :: n, j
integer :: lda, ldu, ldvt
integer, parameter   :: lwmax = 1000 

! matrix in is a
complex(kind=dp), dimension(:,:), intent(inout) :: a
real(kind=dp), dimension(:) :: s
! left & right vectors
complex(kind=dp), allocatable, dimension(:,:) :: u, vt
! eigen values are w
complex(kind=dp), allocatable, dimension(:) ::  work 

! temp scalars
integer :: info, lwork

!temp arrays
real(kind=dp), allocatable, dimension(:) ::  rwork

! use size of input matrix
n=size(a,1)
ldu=size(a,1)
ldvt=size(a,1)
lda=size(a,1)

! eigen vectors, eigen vals & temp arrays
allocate(u( ldu, n ))
allocate(vt( ldvt, n ))
!allocate(w( n ))
allocate( work( lwmax ))
allocate(rwork(2*n))

!no left and right col vectors
! rows m =size(a,1)
! cols n= size(a,2)
! a is matrix
! lda =size(a,1)
! s vector svd
! u matrix
! ldu = m 
! vt matrix hermitian conjg
! ldvt = n 
! work 
! lwork
! rwork
call printvectors(a, 'dens -dens_est')

! quiery the workspace size
lwork = -1
call zgesvd('S','N', n, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info )

! do svd
lwork = min( lwmax, int( work( 1 ) ) )
call zgesvd('S','N', n, n, a, lda, s, u, ldu, vt, ldvt,  work, lwork, rwork, info )

!     cHECK FOR CONVERGENCE.
if( info.gt.0 ) then
write(*,*)'tHE ALGORITHM FAILED TO COMPUTE EIGENVALUES.'
stop
end if

print*, s
end subroutine complexsvd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!matrix norms
function matrixnorm(c)
    complex(kind=dp), dimension(:,:) :: c
    real(kind=dp) :: matrixnorm, zlange
    ! temp
    real(kind=dp), dimension(:), allocatable :: work
    integer :: m,n,lda, lwmax=1000 
 
    m=size(c,1)
    n=size(c,2)
    lda=m
    allocate(work(lwmax))
    matrixnorm= zlange('F', m,n,c,lda,work )
end function matrixnorm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
end module olis_fstdlib
