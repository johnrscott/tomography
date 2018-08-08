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

function complextrace(a,b)
        complex(kind=dp), dimension(:,:) :: a,b
        complex(kind=dp), dimension(:,:), allocatable :: c
        complex(kind=dp) :: complextrace 
        integer :: i

        allocate(c(size(a),size(a)))
        complextrace=0.0_dp
        c=matmul(a,b)
        do i=1, size(a)
                complextrace=complextrace+c(i,i)
        end do
end function complextrace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module olis_fstdlib
