      SUBROUTINE Haar_U_N(n,QH,ROUTIN)
C
C     ########################################################################
C
C     Valerio Cappellini                                            April 2007
C
C     USES :: rgnf_lux, GinUE
C     COMPATIBILITY :: Fortran 90
C     VERSION :: 1.0.0
C
C     This Program produce N x N CPLX Random Unitary Matrices distributed according
C     to the Haar measure, making use of the QR decomposition of N x N CPLX Random (non 
C     Hermitian) Matrices from the Ginibre Ensemble. The QR decomposition is performed
C     by means of (N-1) Householder reflections. 
C
C     Algorithm:
C     ^^^^^^^^^
C      + we start producing an N x N matrix A(0) = a(0)_{ij} from the Ginibre Ensemble.
C
C                H
C      + we fix Q (0) = I_N ( the N x N identity; here and in the following the 
C        superscript H will denote the usual adjoint operation on CPLX matrices)
C
C         + we perform an iterated procedure on the index <k> running from 1 to N-1
C     
C             at each step a matrix H_k is produced and the matrices A(k-1) 
C
C                  H
C             and Q (k-1) are upgrated to 
C
C             A(k-1) -> A(k) = H_k * A(k-1) and 
C
C              H          H             H  
C             Q (k-1) -> Q (k) = H_k * Q (k-1)
C
C         + end of iterated procedures
C
C                                   H
C      + We divide the last raw of Q (N-1) by the unimodular Phase of A(N-1)_{NN}
C
C                    H 
C     Finally QH := Q (N-1) is the CPLX Random Unitary Matrix distributed according 
C     to the Haar measure and 
C          
C     R = QH * A(0) is an upper triangular matrix with REAL positive diagonal entries 
C     so that this upper positive triangularity could be used as a test for the subroutine, 
C     toghether with the Unitarity of QH.
C     
C     At each step the algorithm make use of the N-k+1 dimensional basis vector 
C
C                                    t
C     e_k : (1 , 0 , 0 , ... , 0 , 0)         , (the superscript t denotes the usual 
C                                                transposition) ,
C     of the N-k+1 dimensional CPLX vector
C                                                           t
C     v := (a_{k,k} , a_{k+1,k} , ... , a_{N-1,k} , a_{N,k})   ,
C
C                                               H
C     its Euclidean norm || v || = < v | v > = v * v,  its first entry v_1 = a_{kk}   ,
C
C     the phase of the latter exp(i theta) := v_1 / | v_1 | .
C     
C     Then we construct the vector   u := v + exp(i theta) * || v || * e_k and its 
C
C     positive rescaled-squared-norm c_k := || u ||^2 / 2  = || v ||*(|| v || + | v_1 |)
C     so that finally
C
C
C                /                   |                   \
C               |       I_{k-1}      |         0          |
C               |                    |                    |
C      H_k :=   | ___________________|___________________ |
C               |                    |                    |
C               |                    |                    |
C               |          0         |        M_k         |           
C                \                   |                   /
C     
C     with I_{k-1} being the (k-1) x (k-1) identity , and M_k the (N-k+1) x (N-k+1)
C
C                                                                          H
C                                                       /             u * u   \
C     dimensional CPLX matrix  M_k := -exp(-i theta) * | I_{N-k+1} - --------  |
C                                                       \              c_k    /
C
C     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     ROUTIN   : Name of the subroutine you are going to use for generating
C                the random numbers (suggested :: rgnf_lux.f)
C
C     IMPORTANT [1] !!! Remember to initialize the Subroutine rgnf_lux.f
C                       by means of "Call Ranset_()"
C     IMPORTANT [2] !!! Remember to put the file "i_seed_rng" into the 
C                       working directory.
C     IMPORTANT [3] !!! Remember to insert main_program: "EXTERNAL rgnf_lux"
C
C     ########################################################################
C
$     debug
C
      INTEGER n
      COMPLEX a(n,n),QH(n,n),phase,dk,sum_c,tau
      REAL ck
      REAL VARIANCE,sum_r
      INTEGER i,j,k

      EXTERNAL ROUTIN

      VARIANCE=1./SQRT(2.)
      call GinUE(n,VARIANCE,a,ROUTIN)  ! the matrix A(0) has been produced

      Do i=1,n
          Do j=1,n
              QH(i,j)=(0.,0.)
          Enddo
          QH(i,i)=(1.,0.)
      Enddo                   ! the identity matrix QH(0) has been produced

      sum_r=0.
      do 30 k=1,n-1

        do 11 i=k,n
          sum_r=max(sum_r,cabs(a(i,k)))     ! a(k:n,k) contains v
11      continue
        if(sum_r.ne.0.)then                 ! IF v =/= 0
          sum_r=0.
          do 13 i=k,n
            sum_r=sum_r+cabs(a(i,k))**2
13        continue                           ! sum_r = || v ||^2
      phase=sqrt(sum_r)*a(k,k)/cabs(a(k,k)) ! phase = exp(i theta) * || v ||
          a(k,k)=a(k,k)+phase              ! a(k:n,k) now contains u
          ck=real(conjg(phase)*a(k,k))    ! ck now contains c_k = || u ||^2 / 2
                                         !  = || v ||*( || v || + | v_1 | )
          dk=-conjg(phase)/sqrt(sum_r)  ! dk = - exp(-i theta)

          do 16 j=k+1,n                        !           _N_    _
            sum_c=(0.,0.)                      !          \   |   u_l * a_{l,j}
            do 14 i=k,n                        !   tau =   |     ---------------  
              sum_c=sum_c+conjg(a(i,k))*a(i,j) !          /___|        c_k
14          continue                           !           l=k
            tau=sum_c/ck                  
            do 15 i=k,n                      !
              a(i,j)=dk*(a(i,j)-tau*a(i,k))  !  a(k:n,j) -> M_k * a(k:n,j)
15          continue                         !
16        continue                           

          do 26 j=1,n                           !           _N_    _      H
            sum_c=(0.,0.)                       !          \   |   u_l * Q_{l,j}
            do 24 i=k,n                         !   tau =   |     ---------------  
              sum_c=sum_c+conjg(a(i,k))*QH(i,j) !          /___|        c_k
24          continue                            !           l=k
            tau=sum_c/ck                  
            do 25 i=k,n                        !   H                  H
              QH(i,j)=dk*(QH(i,j)-tau*a(i,k))  !  Q (k:n,j) -> M_k * Q (k:n,j)
25          continue                           !
26        continue                           

        endif
30    continue

      dk=-conjg(a(n,n))/cabs(a(n,n))   !
      Do i=1,n                         !  we multiply the last raw of QH(N-1) 
          QH(n,i)= -dk*QH(n,i)         !  by exp(-i theta)
      Enddo                            !
      
      return
      END
