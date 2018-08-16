      SUBROUTINE GinUE(N,VARIANCE,CMAT,ROUTIN)
C
C     ########################################################################
C
C     Valerio Cappellini                                            March 2007
C
C     USES :: rgnf_lux
C     COMPATIBILITY :: Fortran 90
C     VERSION :: 1.0.0
C
C     Generator of CPXL N x N non Hermitian matrices drawn according to the Ginibre ensemble
C     called <GinUE> \subset GL(N,C). The distribution of matrix Z and its matrix 
C     elements are given by
C                                      | z_{ij} |^2
C                        1          - --------------
C     P(z_{ij}) = --------------  e       2 s^2
C                    2\pi s^2
C
C      and
C
C                       1 
C     P(Z) = ------------------------  exp[ - Tr Z* Z / (2 s^2) ]
C             (2\pi)^(N^2) s^(2 N^2)
C
C     where s denotes the input variable <VARIANCE> and Z the output 
C     variable <CMAT>.
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
      INTEGER N, emme, enne, mu
      COMPLEX SW(0:1), CMAT(N,N)
      REAL U(2), S, T, A, B, R1, R2, V, X, Y, Q, DEVIAT, VARIANCE
      SAVE  S, T, A, B, R1, R2, SW
      EXTERNAL ROUTIN
      DATA  S, T, A, B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA  R1, R2 ,SW / 0.27597, 0.27846, (1. , 0.), (0. , 1.)/
C         generate pair of uniform deviates

      DO 200 emme = 1, N
          DO 200 mu = 2, 2*N+1
      IF(mu.LE.(N+1)) CMAT(emme,mu-1) = (0.,0.)
   50 CALL ROUTIN(U,2)
      V = 1.7156 * (U(2) - 0.5)
      X = U(1) - S
      Y = ABS(V) - T
      Q = X**2 + Y*(A*Y - B*X)
C           accept P if inside inner ellipse
      IF (Q .LT. R1)  GO TO 100
C           reject P if outside outer ellipse
      IF (Q .GT. R2)  GO TO 50
C           reject P if outside acceptance region
      IF (V**2 .GT. -4.0 *ALOG(U(1)) *U(1)**2)  GO TO 50
C           ratio of P's coordinates is normal deviate
  100 DEVIAT = V/U(1)*VARIANCE
  200 CMAT(emme,mu/2) = CMAT(emme,mu/2) + DEVIAT * SW(MOD(mu,2))

      RETURN
      END
