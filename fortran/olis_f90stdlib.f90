module olis_f90stdlib
implicit none

integer, parameter, private :: dp=selected_real_kind(15,300)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! make temp arrays for complex_eigenvects
integer, private :: lda, ldvl, ldvr
integer, parameter, private   :: lwmax = 1000 

! temp scalars
integer, private :: info, lwork

!temp arrays
real(kind=dp), allocatable, dimension(:), private ::  rwork_eigen, rwork_svd
complex(kind=dp), allocatable, dimension(:), private :: work_eigen, work_svd





! make temp arrays for complex_svd


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine alloc_complex_eigenvects(matrix, eigenvals, u, v)
! complex diag
complex(kind=dp), dimension(:,:), intent(in) :: matrix
complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: u,v
complex(kind=dp), dimension(:), allocatable, intent(inout) :: eigenvals
integer :: n

n=size(matrix,1)

! u, v and eigenvals
allocate(u(n,n))
allocate(v(n,n))
allocate(eigenvals(n))

! temp work arrays
allocate( work_eigen( lwmax ))
allocate(rwork_eigen(2*size(matrix,1)))

print*, 'Allocated temp work arrays for DIAG'
end subroutine alloc_complex_eigenvects

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alloc_complex_svd(matrix, sigma, u, vt)
! complex SVD
complex(kind=dp), dimension(:,:), intent(in) :: matrix
complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: u,vt
real(kind=dp), dimension(:), allocatable, intent(inout) :: sigma 
integer :: n,m
m=size(matrix,1)
n=size(matrix,2)

allocate(u(m,n))
allocate(vt( m, n ))
! need to check which is smallest
allocate(sigma(n))

allocate( work_svd( lwmax ))
allocate(rwork_svd(2*size(matrix,1)))

print*, 'Allocated temp work arrays for SVD'
end subroutine alloc_complex_svd

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


subroutine complex_eigenvects(a, w, vl, vr)
! matrix a in, eigenvals, eigenvects out
implicit none
! matrix in is a
!!!!!!!!!!!!!!!!!!!!!!! matrix a is overwritten WATCH OUT!
complex(kind=dp), dimension(:,:), allocatable :: a

! left & right vectors
complex(kind=dp), dimension(:,:), allocatable :: vl
complex(kind=dp), dimension(:,:), allocatable :: vr
! eigen values are w
complex(kind=dp), allocatable, dimension(:) :: w 

integer :: n, m

! use size of input matrix
n=size(a,1)
m=size(a,2)

lda=n
ldvl=n
ldvr=m

!     qUERY THE OPTIMAL WORKSPACE.
lwork = -1
call zgeev( 'V', 'v', n, a, lda, w, vl, ldvl, vr, ldvr, work_eigen, lwork, rwork_eigen, info )
lwork = min( lwmax, int( work_eigen( 1 ) ) )

!     sOLVE EIGENPROBLEM.
call zgeev( 'v', 'v', n, a, lda, w, vl, ldvl, vr, ldvr, work_eigen, lwork, rwork_eigen, info )

!     cHECK FOR CONVERGENCE.
if( info.gt.0 ) then
write(*,*)'tHE ALGORITHM FAILED TO COMPUTE EIGENVALUES.'
stop
end if

end subroutine complex_eigenvects

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine complex_svd(a, sigma, u, vt)
! matrix a in, eigenvals, eigenvects out
! A = U * sigma * V **H
implicit none 

integer :: n, m, i
integer :: lda, ldu, ldvt
integer, parameter   :: lwmax = 1000 

! matrix in is a
complex(kind=dp), dimension(:,:), allocatable, intent(inout) :: a
! singular values sigma
real(kind=dp), dimension(:), allocatable :: sigma
! left & right vectors
complex(kind=dp), dimension(:,:), allocatable :: u, vt
! eigen values are w
!complex(kind=dp), allocatable, dimension(:) ::  work 

! temp scalars
!integer :: info, lwork

!temp arrays
!real(kind=dp), allocatable, dimension(:) ::  rwork

! use size of input matrix
! might have got n & m the wrong way round
n=size(a,1)
m=size(a,2)
ldu=size(a,1)
ldvt=size(a,1)
lda=size(a,1)

! eigen vectors, eigen vals & temp arrays
!allocate(u( ldu, n ))
!allocate(vt( ldvt, n ))
!allocate(w( n ))
!allocate( work( lwmax ))
!allocate(rwork(2*n))

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
call zgesvd('S','S', m, n, a, lda, sigma, u, ldu, vt, ldvt, work_svd, lwork, rwork_svd, info )

! do svd
lwork = min( lwmax, int( work_svd( 1 ) ) )
call zgesvd('S','S', m, n, a, lda, sigma, u, ldu, vt, ldvt,  work_svd, lwork, rwork_svd, info )

!     cHECK FOR CONVERGENCE.
if( info.gt.0 ) then
write(*,*)'tHE ALGORITHM FAILED TO COMPUTE EIGENVALUES.'
stop
end if

! write singular vals back to matrix a
a=0.0_dp
do i =1, n
    a(i,i)=sigma(i)
end do
end subroutine complex_svd

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
end module olis_f90stdlib
