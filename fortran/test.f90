program test
    use olis_f90stdlib
    implicit none
    
    integer :: n
    integer, parameter :: dp=selected_real_kind(15,300)
    complex(kind=dp), dimension(:,:), allocatable :: a 
    complex(kind=dp), dimension(:,:), allocatable :: b
    ! for svd and diagonalisation
    complex(kind=dp), dimension(:,:), allocatable ::u_diag, v_diag, u_svd, vt_svd
    complex(kind=dp), dimension(:), allocatable :: eigenvals
    real(kind=dp), dimension(:), allocatable :: sigma

    n=2
    allocate(a(n,n))
    
    ! do this at the start !
    call alloc_complex_eigenvects(a, eigenvals, u_diag, v_diag)
    call alloc_complex_svd(a, sigma, u_svd, vt_svd) 
    a=0.0_dp
    a(1,2)=(1.0_dp,0.0_dp)
    a(1,1)=(1.0_dp,0.0_dp)
    a(2,1)=(1.0_dp, 0.0_dp)
    b=a

    call printvectors(b, 'my matrix happens to be X')
    ! return diag matrix and u and v**H
    call complex_eigenvects(b,eigenvals,u_diag,v_diag)

    call printvectors(a, 'old matrix')
    call printvectors(b, 'new matrix')
    call printvectors(u_diag, 'u')
    call printvectors(v_diag, 'v')
    write(*,*) 'eigen vals', eigenvals

    call printvectors(matmul(u_diag,(matmul(b,v_diag))), 'matmul u, (a,v)')

    ! now test SVD
    
    print*, '--------------------------------------------'
    b=a
    !call printvectors(b, 'reset b')
    call complex_svd(b,sigma,u_svd,vt_svd)
   
    call printvectors(a, 'old matrix')
    call printvectors(b, 'matrix')
    !call printvectors(u_svd, 'u')
    !call printvectors(vt_svd, 'vt')
    write(*,*) 'singular vals', sigma

    ! vt * a * u = singular vals
    ! u * singular vals * vt = a
    print*,
    call printvectors(matmul(vt_svd, matmul(a,u_svd)), 'singular vals = matmul vt, (a,u)')
    call printvectors(matmul(u_svd, matmul(b,vt_svd)), 'a = matmul u, (singular vals, vt)')
end program test
