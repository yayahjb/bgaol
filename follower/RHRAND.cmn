
C**********************************************************************C
      SUBROUTINE INITRN (ISEED,INDX,JNDX,KNDX,BUFFER)
C   use a simple linear congruential random number generator to
C   fill the buffer array to start up a linear shift register
C   generator
C   set the integers indx,jndx,kndx to the "magic" initial
C   values described by Knuth, The Art of Computer Programming,
C   vol 2, Seminumerical Algorithms

C   set iseed negative to indicate that the initialization was done
      implicit none

      real*8 BUFFER(56)
      real*8 bufx
      INTEGER ISEED, INDX, JNDX, KNDX
C----------------------------------------------------------------------C

      JNDX = ISEED
      DO 1000 INDX= 1,56
         JNDX = MOD(JNDX*2349 + 14867,32767)
         bufx = JNDX
        BUFFER(INDX) = dABS(bufx/DBLE(32767))
 1000 CONTINUE
      ISEED = -ISEED
      IF (ISEED .EQ. 0) ISEED = -1
      INDX = 55
      KNDX = 54
      JNDX = 31
      END

C**********************************************************************C
      real*8 FUNCTION RANDF( ISEED )
C   real*8 function randf( iseed )
C        if seed .gt. 0 use seed to initialize
C        returns a uniformly distributed random number in (0. < 1.)
C
      implicit none

      INTEGER ISEED
C
      INTEGER INDX,JNDX,KNDX
      real*8    BUFFER(56)
      SAVE BUFFER, INDX, JNDX, KNDX
      DATA INDX /-1/
C----------------------------------------------------------------------C

      IF( ISEED .GE. 0 .OR. INDX .LT. 0)
     2   CALL INITRN (ISEED,INDX,JNDX,KNDX,BUFFER)

      INDX = MOD(INDX ,55)+1
      JNDX = MOD(JNDX ,55)+1
      KNDX = MOD(KNDX ,55)+1
      BUFFER(INDX) = DMOD( BUFFER(JNDX) + BUFFER(KNDX), 1.0D0)
      RANDF = BUFFER(INDX)
      END

C**********************************************************************C
      REAL*8 FUNCTION RANDFG( ISEED )
C   real function randfg( iseed )
C        if seed .gt. 0 use seed to initialize
C        returns a normally distributed random number with unit variance
C
C   EXACTLY FOLLOWS THE ALGORITHM OF KNUTH, THE ART OF COMPUTER
C   PROGRAMMING, V. 22D. ED., 1981, P. 117-118  (THE POLAR METHOD)
C   EXCEPT THAT LIMITS ARE INCLUDED BY MAKING THE MINIMUM VALUE OF "S"
C   BE EPS
C
      implicit none

      INTEGER ISEED
C
      INTEGER INDX,JNDX,KNDX
      REAL*8    BUFFER(56)
      REAL*8    EPS, R1, R2, S
      integer i
      SAVE BUFFER, INDX, JNDX, KNDX,EPS
      DATA INDX /-1/
C----------------------------------------------------------------------C

      IF( ISEED .GE. 0 .OR. INDX .LT. 0) THEN
         R1 = 1.0D0
         DO 100 I=1,100
             EPS = R1
             R1 = R1/2
             IF( R1 + 1.0D0 .EQ. 1.0D0) GOTO 101
100      CONTINUE
101      CONTINUE
         CALL INITRN (ISEED,INDX,JNDX,KNDX,BUFFER)
      ENDIF

 4000 CONTINUE
      INDX = MOD(INDX ,55)+1
      JNDX = MOD(JNDX ,55)+1
      KNDX = MOD(KNDX ,55)+1
      BUFFER(INDX) = DMOD( BUFFER(JNDX) + BUFFER(KNDX), 1.0D0)
      R1 = 2.0D0 * BUFFER(INDX) - 1.0D0
      INDX = MOD(INDX ,55)+1
      JNDX = MOD(JNDX ,55)+1
      KNDX = MOD(KNDX ,55)+1
      BUFFER(INDX) = DMOD( BUFFER(JNDX) + BUFFER(KNDX), 1.0D0)
      R2 = 2.0D0 * BUFFER(INDX) - 1.0D0
      S = R1*R1 + R2*R2
      IF (S .GE. 1.0D0) GO TO 4000
      IF (ABS(S) .LT. EPS) S = DSIGN(EPS,S)
      RANDFG = R1 * DSQRT(-2.0D0 * DLOG(S)/S)
      END

C**********************************************************************C
      subroutine fgbrand( iseed, x, nc )
      implicit none

C     This function returns a random direction 
C     vector v = (x_1,x_2,...,x_n) in n dimensions. 
C     The vector is NOT normalized to a unit vector.
C     The method uses the fact that a multivariate Gaussian 
C     distribution is spherically symmetric. 
C     Each component is generated to have a Gaussian distribution, 
C     The method is described 
C     by Knuth, v2, 3rd ed, p135,136, and attributed to G. W. Brown, 
C     Modern Mathematics for the Engineer (1956).

C     To obtain a random distribution on an S(n) sphere, normalize
C     the output to a unit vector.

       INTEGER ISEED,NC
       real*8 X(NC)
C
       INTEGER INDX,JNDX,KNDX
       real*8    BUFFER(55)
       SAVE BUFFER,INDX,JNDX,KNDX
       DATA INDX /-1/
       real*8 randfg
       integer i
C----------------------------------------------------------------------C
      DO 1000 I=1,NC
         x(i) = RANDFG( ISEED )
1000  CONTINUE
      END

C**********************************************************************C
       SUBROUTINE BRAND( ISEED, X, NC )
C subroutine brand( iseed,x,nc)
C  generate nc random numbers in real*8 array x
C
C derived from
C   real*8 function randf( iseed )
C        if seed .gt. 0 use seed to initialize
C        returns NC uniformly distributed random numbers (-1. < 1.)
C
      implicit none

       INTEGER ISEED,NC
       real*8 X(NC)
C
       INTEGER INDX,JNDX,KNDX
       real*8    BUFFER(56)
       SAVE BUFFER,INDX,JNDX,KNDX
       DATA INDX /-1/
       integer i
C----------------------------------------------------------------------C

      IF( ISEED .GT. 0 .OR. INDX .LT. 0)
     2     CALL INITRN (ISEED,INDX,JNDX,KNDX,BUFFER)

      DO 1000 I=1,NC
         INDX = MOD(INDX ,55)+1
         JNDX = MOD(JNDX ,55)+1
         KNDX = MOD(KNDX ,55)+1
         BUFFER(INDX) = DMOD( BUFFER(JNDX) + BUFFER(KNDX), 1.0D0)
         X(I) = 2.0D0*BUFFER(INDX)-1.0D0
1000  CONTINUE
      END


C**********************************************************************C
      LOGICAL FUNCTION ISUNIT(A)
      implicit none

      real*8 A(6,6)
      integer i, j
C----------------------------------------------------------------------C
      ISUNIT = .FALSE.
      DO 1000 I=1,6
      DO 1000 J=1,6
         IF (I.EQ.J) THEN
            IF (DABS(A(I,I)-1.0D0) .GT. 1.0D-5) RETURN
         ELSEIF(DABS(A(I,J)).GT. 1.0D-7) THEN
            RETURN
         ENDIF
 1000 CONTINUE
      ISUNIT = .TRUE.
      END
