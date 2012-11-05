C   program to investigate the environment around the 8-fold
C   corner in Niggli space. This is cF (face-center cubic),
C   and the environment is very complicated.

      include "BasicReduction.For"
      include "MKREFL.FOR"
      include "../NEAR.FOR"
      include "../E3TOG6.FOR"
      include "../MKGAOL.FOR"

C**********************************************************************C
      program Follower
      implicit none

      common/xdebug/xdebug
      logical xdebug

      real*8 v(6,42)
      character*3 roof(42)
      character*30 absymbols(42)
      integer nationaltables(42)

      character*8 chprj(215)
      character*8 level8

      real*8 stack(8),vp(8)

      real*8 prj(36,18)
      real*8 mred1(36)
      real*8 m1(36), m2(36), m3(36)
      real*8 prjout(36)

      real*8 v1(6), v2(6), v3(6), v4(6), v5(6), vran(6)
      real*8 c2(6)
      real*8 testCells(6,44)
      logical NextIndex
      real*8 base(6), rv1(6), rv2(6)
      real*8 ptrb6(6)
      real*8 temp(36)
      integer icount(20000), iprjComponent(2000)
      real*8 xforms(36,2000)
      integer nprjOut(10)
      logical reduce
      logical isunit
      logical logtest
      logical okcell
      logical quiet

      integer MXTREE
      integer i, j, iter, maxiter, nprevTrial, iseed, nprj, itemp
      integer nmatprv, nmat
      integer nprjTotal
      parameter (MXTREE=50*2000)
      real*8 TREE(MXTREE)
      integer FillPrjList
      real*8 eightfoldProjectors(36,1000)
      integer acceptedTrials
      integer nxform, ncycle
      integer xformList(2000)


      real*8 delta, anorm
      parameter (delta=0.1D0)
      parameter (maxiter=10000)

      data quiet /.TRUE./
      data iseed /137459/
      data xdebug/.false./

C----------------------------------------------------------------------
C The symbols for Bravais lattice types from Roof and Niggle
      Data roof /
     * "44A","44B","44C","45A","45B","45C","45D","45D",
     * "45E","48A","48B","49B","49C","49D","49E","50C",
     * "50A","50B","50D","50E","50F","51A","51B","52A",
     * "52B","52C","53A","53B","53C","54A","54B","54C",
     * "55A","55A","55B","55B","56A","56B","56C","57A",
     * "57B","57C" /

C the vector representations of each of the 42, non-triclinic
C bravais lattice types. The symbols such as r,r,r,s,s,s from
C Andrews and Bernstein, 1988 were converted to instantiated
C G6 vectors

      data (v(i, 1),i=1,6) /11,11,11,0,0,0/
      data (v(i, 2),i=1,6) /11,11,11,-7.33333,-7.33333,-7.33333/
      data (v(i, 3),i=1,6) /11,11,11,11,11,11/
      data (v(i, 4),i=1,6) /11,11,13,0,0,0/
      data (v(i, 5),i=1,6) /11,13,13,0,0,0/
      data (v(i, 6),i=1,6) /11,11,13,-11,-11,0/
      data (v(i, 7),i=1,6) /11,11,11,-24,-24,26/
      data (v(i, 8),i=1,6) /11,11,11,26,-24,-24/
      data (v(i, 9),i=1,6) /11,13,13,5.5,11,11/
      data (v(i,10),i=1,6) /11,11,13,0,0,-11/
      data (v(i,11),i=1,6) /11,13,13,-13,0,0/
      data (v(i,12),i=1,6) /11,11,13,11,11,11/
      data (v(i,13),i=1,6) /11,11,11,13,13,13/
      data (v(i,14),i=1,6) /11,11,11,13,13,13/
      data (v(i,15),i=1,6) /11,13,13,-9.33333,-7.33333,-7.33333/
      data (v(i,16),i=1,6) /11,13,17,0,0,0/
      data (v(i,17),i=1,6) /11,13,17,0,-11,0/
      data (v(i,18),i=1,6) /11,13,17,0,0,-11/
      data (v(i,19),i=1,6) /11,11,13,0,0,17/
      data (v(i,20),i=1,6) /11,13,13,17,0,0/
      data (v(i,21),i=1,6) /11,13,17,-13,0,0/
      data (v(i,22),i=1,6) /11,11,13,17,17,-56/
      data (v(i,23),i=1,6) /11,13,17,5.5,11,11/
      data (v(i,24),i=1,6) /11,11,11,13,17,-52/
      data (v(i,25),i=1,6) /11,13,13,17,11,11/
      data (v(i,26),i=1,6) /11,13,17,-13,-11,0/
      data (v(i,27),i=1,6) /11,13,17,0,19,0/
      data (v(i,28),i=1,6) /11,13,17,19,0,0/
      data (v(i,29),i=1,6) /11,13,17,0,0,19/
      data (v(i,30),i=1,6) /11,13,17,19,0,-11/
      data (v(i,31),i=1,6) /11,13,17,-13,19,0/
      data (v(i,32),i=1,6) /11,13,17,19,-11,0/
      data (v(i,33),i=1,6) /11,11,13,17,17,19/
      data (v(i,34),i=1,6) /11,11,13,17,17,19/
      data (v(i,35),i=1,6) /11,13,13,17,19,19/
      data (v(i,36),i=1,6) /11,13,13,17,19,19/
      data (v(i,37),i=1,6) /11,13,17,19,11,38/
      data (v(i,38),i=1,6) /11,13,17,13,19,38/
      data (v(i,39),i=1,6) /11,13,17,19,38,11/
C      data (v(i,40),i=1,6) /11,13,17,-32,-30,38/
      data (v(i,40),i=1,6) /10.0,19.9996,30.0003,-19.9996,0,0/
      data (v(i,41),i=1,6) /11,11,13,17,19,-58/
      data (v(i,42),i=1,6) /11,13,17,19,11,11/


C,F,cF,IT(01)
      data (testCells(i,01),i=1,6)/14.142,14.142,14.142,90,90,90/
C,R,hR,IT(02)
      data (testCells(i,02),i=1,6)/11.180,11.180,22.913,90,90,120/
C
C
C,I,cI,IT(05)
      data (testCells(i,05),i=1,6)/11.547,11.547,11.547,90,90,90/
C,I,tI,IT(06)
      data (testCells(i,06),i=1,6)/11.180,11.180,12.247,90,90,90/
C,I,tI,IT(07)
      data (testCells(i,07),i=1,6)/11.180,12.247,11.180,90,90,90/
C,I,oI,IT(08)
      data (testCells(i,08),i=1,6)/11.832,10.488,12.247,90,90,90/
C,R,hR,IT(09)
      data (testCells(i,09),i=1,6)/10,10,38.730,90,90,120/
C,C,mC,IT(10)
      data (testCells(i,10),i=1,6)/73.728,67.559,187.305,90,91.076,90/
C,P,tP,IT(11)
      data (testCells(i,11),i=1,6)/10,10,14.142,90,90,90/
C,P,hP,IT(12)
      data (testCells(i,12),i=1,6)/10,10,14.142,90,90,120/
C,C,oC,IT(13)
      data (testCells(i,13),i=1,6)/11.180,16.583,14.142,90,90,90/
C,C,mC,IT(14)
      data (testCells(i,14),i=1,6)/10.724,16.882,12.845,90,106.881,90/
C,I,tI,IT(15)
      data (testCells(i,15),i=1,6)/10,10,24.495,90,90,90/
C,F,oF,IT(16)
      data (testCells(i,16),i=1,6)/15.811,25.495,12.247,90,90,90/
C,I,mI,IT(17)
      data (testCells(i,17),i=1,6)/15.492,10.488,15.811,90,114.095,90/
C,I,tI,IT(18)
      data (testCells(i,18),i=1,6)/18.708,18.708,10,90,90,90/
C,I,oI,IT(19)
      data (testCells(i,19),i=1,6)/10,19.365,18.028,90,90,90/
C,C,mC,IT(20)
      data (testCells(i,20),i=1,6)/21.794,18.028,10,90,112.955,90/
C,P,tP,IT(21)
      data (testCells(i,21),i=1,6)/14.142,14.142,10,90,90,90/
C,P,hP,IT(22)
      data (testCells(i,22),i=1,6)/14.142,14.142,10,90,90,120/
C,C,oC,IT(23)
      data (testCells(i,23),i=1,6)/18.028,21.794,10,90,90,90/
C,R,hR,IT(24)
      data (testCells(i,24),i=1,6)/16.583,16.583,8.660,90,90,120/
C,R,hT,IT(24)
      data (testCells(i,24),i=1,6)/23.805,23.805,10,90,90,120/
C,C,mC,IT(25)
      data (testCells(i,25),i=1,6)/18.028,21.794,10,90,118.131,90/
C,F,oF,IT(26)
      data (testCells(i,26),i=1,6)/26.458,33.166,10,90,90,90/
C,I,mI,IT(27)
      data (testCells(i,27),i=1,6)/20.616,10,21.794,90,102.860,90/
C,C,mC,IT(28)
      data (testCells(i,28),i=1,6)/10,33.166,14.142,90,108.554,90/
C,C,mC,IT(29)
      data (testCells(i,29),i=1,6)/10,26.458,17.321,90,105.059,90/
C,C,mC,IT(30)
      data (testCells(i,30),i=1,6)/14.142,31.623,10,90,108.554,90/
C
C,P,oP,IT(32)
      data (testCells(i,32),i=1,6)/10,14.142,17.321,90,90,90/
C,P,mP,IT(33)
      data (testCells(i,33),i=1,6)/10,14.142,17.321,90,102.504,90/
C,P,mP,IT(34)
      data (testCells(i,34),i=1,6)/10,17.321,14.142,90,105.377,90/
C,P,mP,IT(35)
      data (testCells(i,35),i=1,6)/14.142,10,17.321,90,98.806,90/
C,C,oC,IT(36)
      data (testCells(i,36),i=1,6)/10,33.166,14.142,90,90,90/
C
C,C,oC,IT(38)
      data (testCells(i,38),i=1,6)/10,26.458,17.321,90,90,90/
C
C,C,oC,IT(40)
      data (testCells(i,40),i=1,6)/14.142,31.623,10,90,90,90/
C
C,I,oI,IT(42)
      data (testCells(i,42),i=1,6)/10,14.142,30,90,90,90/
C,I,mI,IT(43)
      data (testCells(i,43),i=1,6)/10,31.305,14.142,90,106.430,90/
C

C the boundary symbols of Andrews and Bernstein
C (probably 2012). There should be 215, up thru
C the 8-fold boundary
      Data chprj /
     * "1","2","3","4","5","6","7","8","9","A",
     * "B","C","D","E","F",
     * "12","13","14","15","16","17","18","19","1A","1B",
     * "1C","1D","1E","1F","23","24","25","26","27","28",
     * "29","2A","2B","2C","2D","2E","2F",
     * "34","35","3A","3B","3D","3E","45","47","48","4C",
     * "4E","56","58","59","5B","67","69","7C","8F","9A",
     * "9C","AD","BE","BF","CD","EF",
     * "123","124",
     * "125","126","127","128","129","12A","12B","12C",
     * "12D","12E","12F","134","135","13A","13B","13D",
     * "13E","145","147","148","14C","14E","156","158",
     * "159","15B","167","169","17C","18F","19A","1AD",
     * "1BF","1CD","1EF","234","235","23A","23B","23D",
     * "23E","245","247","248","24C","24E","256","258",
     * "259","25B","267","269","27C","28F","29A",
     * "29C","2AD","2BE","2BF","2CD","2EF","345","34E",
     * "35B","3AD","3BE","458","47C","569","BEF",
     * "1234","1235","123A","123B","123D","123E","1245","1247",
     * "1248","124C","124E","1256","1258","1259","125B",
     * "1267","1269","127C","128F","129A","12AD","12BF",
     * "12CD","12EF","1345","134E","135B","13AD","1458",
     * "147C","1569","2345","234E","235B","23AD","23BE",
     * "2458","247C","2569","2BEF","34CD","359A","4567",
     * "48EF","58BF","679C","9ACD",
     * "12345","1234E",
     * "1235B","123AD","12458","1247C","12569","134CD",
     * "1359A","13BEF","14567","148EF","158BF","234CD",
     * "2359A","24567","248EF","258BF","2679C","29ACD",
     * "1234CD","12359A","123BEF","124567","1248EF",
     * "1258BF",
     * "1679ACD",
     * "12679ACD" /

C the actual Andrews and Bernstein(1988) symbols for Bravais/Roof/Niggli lattice types
      Data absymbols /
     * "r,r,r,0,0,0","r,r,r,-2*r/3,-2*r/3,-2*r/3","r,r,r,r,r,r",
     * "r,r,s,0,0,0","r,s,s,0,0,0","r,r,s,-r,-r,0",
     * "r,r,r,-r-s,-r-s,2*s","r,r,r,2*s,-r-s,-r-s","r,s,s,r/2,r,r",
     * "r,r,s,0,0,-r","r,s,s,-s,0,0","r,r,s,r,r,r","r,r,r,s,s,s",
     * "r,r,r,s,s,s","r,s,s,-s+r/3,-2r/3,-2r/3","r,s,t,0,0,0",
     * "r,s,t,0,-r,0","r,s,t,0,0,-r","r,r,s,0,0,t","r,s,s,t,0,0",
     * "r,s,t,-s,0,0","r,r,s,t,t,-2*r-2*t","r,s,t,r/2,r,r",
     * "r,r,r,s,t,-2*r-s-t","r,s,s,t,r,r","r,s,t,-s,-r,0"
     * ,"r,s,t,0,u,0","r,s,t,u,0,0","r,s,t,0,0,u","r,s,t,u,0,-r",
     * "r,s,t,-s,u,0","r,s,t,u,-r,0","r,r,s,t,t,u","r,r,s,t,t,u",
     * "r,s,s,t,u,u","r,s,s,t,u,u","r,s,t,u,r,2*u","r,s,t,s,u,2*u",
     * "r,s,t,u,2*u,r","r,s,t,-s-u,-r-u,2*u","r,r,s,t,u,-2*r-t-u",
     * "r,s,t,u,r,r"  /

C the designations for the Bravais/Roof/Niggli lattice types
C as given in the International Tables
       Data nationaltables /
     * 3,5,1,11,21,15,6,7,18,12,22,9,2,4,24,32,
     * 36,38,13,23,40,16,26,8,19,42,33,35,34,
     * 39,41,37,10,14,20,25,28,30,29,43,17,27 /

      data base /10.0,10.0,10.0, 10.0,10.0,10.0/
C      data base /100,980.003025, 199.996164,  0, -79.99951599,  0/
      data level8 /"12679ACD"/

      integer*4 now(3)
C----------------------------------------------------------------------C

      open(  8, file="F_Summary.txt")
      open(  9, file="F_Step_&_NCDIST.txt")
      open( 10, file="F_LargeNCDIST.txt")
      open( 11, file="F_Step.txt")

      call itime( now )
      write(*,"('Time ',i2.2,':',i2.2,':',i2.2)") now

C get the list of all the projectors in the 8-fold corner
      nprjTotal = FillPrjList( eightfoldProjectors )

      nprevTrial = -1


C make a bunch of tests, creating a vector and looping over all the projectors

      nmat = 0
      acceptedTrials = 0
      write(*,*) " before set tree"
           !$OMP PARALLEL DO !I is private by default
      do i=1,mxtree
         tree(i) = -19191.0
      enddo
      write(*,*) " after set tree"

      do 9000 iter=1,maxiter
         nmatprv = nmat
C start from the same base vector each time
         call cpyvn( 6, base, v2 )
C         call cpyvn( 6, v(1,03), v2 )
         call cpyvn( 6, v2, v3 )
         if ( .not. OKCELL( "V", v3, quiet ) ) then
            write(*,*) "initial vector is invalid"
            write(*,*) v2
            write(*,*) v3
            read(*,*)
         endif

C nxform is the number of transfrom matrices in the transformation
C ncycle is the number of reduction cycles to reduce the cell
C delta is the fraction of the (near unit) perturbing vector to add
C nprj is the total number of projectors
C nmat is the number of transform (total) matrices found
C icount is the count of occurences of a transform
         nxform = 0
         ncycle = -1
         call prospect( v2, MXTREE, TREE, acceptedTrials, delta,
     *           nprj, nprevTrial, nmat, icount, xforms, iprjComponent,
     *  xformList, nxform, ncycle )

C if we got a new one, write it out
      if ( nmat .ne. nmatprv ) then
         write(*,"('nmat=',i4, ' ncycle=',i12,' xform',100i3)")
     *   nmat, ncycle, (xformList(j),j=1,nxform)
      endif

 9000 continue

      write(*,*) "nmat,acceptedTrials ", nmat,acceptedTrials
      write(1,*) "nmat,acceptedTrials ", nmat,acceptedTrials

      do 9200 i=1,nmat
         write(*,*) i, icount(i)
         write(1,*) i, icount(i)
 9200 continue

C sort the output data by occurrences
      do 9500 j=1,nmat
      do 9500 i=1,nmat-1
         if ( icount(i) .lt. icount(i+1) ) then

            call cpyvn( 36, xforms(1,i), temp )
            call cpyvn( 36, xforms(1,i+1), xforms(1,i) )
            call cpyvn( 36, temp, xforms(1,i+1) )

            itemp = icount(i)
            icount(i) = icount(i+1)
            icount(i+1) = itemp

            itemp = iprjComponent(i)
            iprjComponent(i) = iprjComponent(i+1)
            iprjComponent(i+1) = itemp
         endif

 9500 continue

      write(*,*)
      write(*,*)

      do 9600 i=1,nmat
         write(*,"(' # ',i4,', counts ',i8 )")
     *          i, icount(i)
         write(*,"(6F8.4)") (xforms(j,i),j=1,36)
         write(*,*)
 9600 continue

C sort the output data by found projector
      do 9700 j=1,nmat
      do 9700 i=1,nmat-1
         if ( iprjComponent(i) .gt. iprjComponent(i+1) ) then

            call cpyvn( 36, xforms(1,i), temp )
            call cpyvn( 36, xforms(1,i+1), xforms(1,i) )
            call cpyvn( 36, temp, xforms(1,i+1) )

            itemp = icount(i)
            icount(i) = icount(i+1)
            icount(i+1) = itemp

            itemp = iprjComponent(i)
            iprjComponent(i) = iprjComponent(i+1)
            iprjComponent(i+1) = itemp
         endif

 9700 continue

      write(*,*)
      write(*,*)

      do 9800 i=1,nmat
         write(*,"(' # ',i4,', counts ',i18, ', from projector ',i3 )")
     *          i, icount(i), iprjComponent(i)
         write(*,"(6F8.4)") (xforms(j,i),j=1,36)
         write(*,*)
 9800 continue

      end

C**********************************************************************C
      subroutine MakePerp( amat, v )
      real*8 amat(6,6), v(6)
      integer i,j
C----------------------------------------------------------------------C

      do 3000 i=1,6
      do 3000 j=1,6
         amat(i,j) = -(v(i)*v(j))
 3000 continue

      do 4000 i=1,6
 4000 amat(i,i) = 1.0 + amat(i,i)
      return
      end

C**********************************************************************C
      subroutine prospect( v2, MXTREE, TREE, acceptedTrials, delta,
     * nprj, nprevTrial, nmat, icount, xforms, iprjComponent,
     * xformList, nxform, ncycle )
      implicit none

      integer maxSteps
      parameter (maxSteps=2000)

      real*8 v2(6)
      real*8 anorm
      integer MXTREE, nprj, nprevTrial, nmat
      real*8 TREE(MXTREE)
      real*8 v3(6), vran(6), C2(6), rv2(6), rvtemp(6)
      integer imred1(36)
      real*8 mred1(36)
      real*8 delta
      integer icount(*), iprjComponent(*)
      real*8 xforms(36,*)
      integer acceptedTrials
      logical ISUNIT
      integer ip, id, iseed, mxtests, itests, k
      integer inear, isum, ncycle, nm
      real*8 sum
      real*8 perp(36),perp66(6,6), vperp(6)
      equivalence (perp,perp66)
      real*8 vtemp(6)
      real*8 distance(2000), NCDIST
      real*8 factor, distance2, distancesource, maxdelta, maxvalue
      integer i,j
      integer NEARST_N
      integer nsteps
      integer nxform
      integer xformList(2000)
      integer step1, step2
      integer im
      common /xdebug/xdebug
      logical xdebug
      logical OKCELL, REDUCE, bTest,quiet
      real*8 cell(6)
      integer icell
      integer nbuffer
      character*200 buffer(maxSteps)
      character*1 tell
      integer ibuff, dev

      real*8 deltaDist, absSum, absSumSQ, avg, stddev
      logical iprintStepDist, iprintNCDIST_Too_Large

C make the seed for starting the random number generator
C "SAVE" it so that it doesn't reinitialize every time this
C subroutine is called.
      data iseed /19191/
      data quiet /.true./
      SAVE iseed

C the following is from HJB emails
C
C      7 - a boundary manifold
C      7'- a manifold of the same order as 7, which intersects 7 to generate 7^
C      7^ - a sub-manifold of 7 consisting of special position points, i.e.
C      eigenvectors of eigenvalue 1 of M_7 in 7.
C      7-alt a manifold of the same order as 7, which intersects 7 along the
C      intersection of the closure of 6 with the closure of 7.
C    There are a lot of interesting distances.  For 7, the distances to
C
C    7 and 7^ are most important.  7^ is the intersection of 7 and 7'.
C    The there is 7 versus 7alt.  The intersection of 6 and 7 is
C    actually the intersection of 7 and 7alt.
C
C    On the distances from g to 7, 7' and 7^.  The distance from g to 7'
C    can be less than the distance from g to 7, but, if I did this right,
C    the distance from g to 7^ should always be greater than or equal
C    to the distance from g to 7.  I've posted a pdf with my first pass
C    at the necessary ^ projectors and the perps for those projectors
C    (makes the distance calculation easier) at
C
C       http://arcib.dowling.edu/~bernsteh/Embedding_G6.pdf
C
C    Note that that does not use 7alt, etc. ...
C
C    On the 67,9A,CD issue one of each pair should have a
C    valid othogonal distance, and the other should have
C    a larger distance that comes in at an angle, as the
C    distance, e.g. to distance to 7 intersect 7alt.  I'll
C    work up those three projectors later today and add them
C    to the PDF.  So the reporting would be done by
C    calacuting the orthogonal distance to 7, which could
C    either be the distance to 6 or the distance to 7, and
C    then computing the distance to 7 intersect 7alt, and
C    then deciding which of the two goes with 6 and which
C    with 7 on the basis of which side of 7alt it lies.
C
C----------------------------------------------------------------------C

C now perturb the projected vector (on the boundary) to search

      call fgbrand (iseed,vran,6)

      do 1000 i=1,6
 1000 v3(i) = v2(i)/anorm(6,v2)

      call MakePerp( perp, v3 )

C      write(*,*) "perp", perp
      call rmv6( vran, perp66, vperp )
C      write(*,*) "v2 ",v2
C      write(*,*) "vran ",vran
C      write(*,*) "vtemp ",vtemp

      do 1500 i=1,6
 1500    vran(i) = vperp(i)/anorm(6,vperp)
C      write(*,*) "vran ",vran


 1999 continue

        !$OMP PARALLEL DO !I is private by default
      do 3000 k=1,6
 3000 v3(k) = v2(k) + delta*(vran(k))
      call unitmx( 6, mred1 );
C      write(*,"(a,6F10.5)") "v3 ",v3

      ncycle = -1
      bTest = reduce (v3, imred1,rv2, 0.0D0, 'REDUCE  ')
      write(*,"(10x,6i3)" ) imred1
      write(*,*)
      do i=1,36
         mred1(i) = imred1(i)
      enddo

      call G6TOC(RV2,C2, quiet, 'G6TOC ')

C Determine if this is a new transform or one we've seen before

      write(*,"(6f4.0)" ) mred1
      write(*,*)
      iNear =
     2   NEARST_N(.FALSE.,MXTREE,36,mred1, 0.001D0,TREE,
     2   IP,ID,'NEARST_N')

      write(*,*) " iNear ", iNear

      IF ( iNear .LE. 0 ) THEN
C       new transformation, store it
         acceptedTrials = acceptedTrials + 1
         nmat = nmat + 1
         call cpyvn( 36, mred1, xforms(1,nmat) )
         iprjComponent(nmat) = nprj
         id = nmat
         nm = nmat
         CALL BLDTRE_N(.FALSE.,MXTREE,36,mred1,nm,TREE,'BLDTRE_N')
         ICOUNT(nmat) = 1
         write(1,*) "nprj, nmat, cycles ",
     *        nprj, nmat, ncycle
         write(1,"(6f8.3)" ) mred1
         write(1,"('('i4')'3x,6f9.4,5x,6f9.4)") ID, v3, rv2
         write(1,*)
      ELSEIF (ID .EQ. 0) THEN
         write(1,*) "failed reduce, nmat, cycles ", nmat, ncycle
         GO TO 9000
      ELSE
C seen this one before - count it
         acceptedTrials = acceptedTrials + 1
         ICOUNT(ID) = ICOUNT(ID) + 1
      ENDIF

      if ( nprevTrial .ne. nmat .or.
     *   mod(acceptedTrials,1000000) .eq. 0  ) then
         write(*,*) "trials, nmat ", acceptedTrials,nmat
         write(1,*) "trials, nmat ", acceptedTrials,nmat
         nprevTrial = nmat
         call FLUSH( )
      endif

C   write out some sample vectors
      if ( icount(id) < 31 ) then
         write(2,"('('i4')'3x,6f9.4,5x,6f9.4)") ID, v3, rv2
         call FLUSH( )


C the simple set of numbers is 101,1,100
         nsteps = 100
         step1 = 1
         step2 = 100

C         xdebug = ID .eq. 185

C            if (ID .eq. 185) then
         ibuff = 0
         iprintNCDIST_Too_Large = .false.
         iprintStepDist   = .false.
         do 8000 i=step1,step2
            factor = DBLE(i-1)/DBLE(nsteps-1)
            do 7000 j=1,6
               vtemp(j) = (1D0-factor)*v3(j) + factor*rv2(j)
 7000       continue

            do 7500 im=1,36
 7500       imred1(im) = 0.0D0
            do 7600 icell=1,6
               cell(icell) = vtemp(icell)
 7600       continue
            if ( OKCELL( "V", cell, .true. ) ) then
            if ( reduce (vtemp, imred1,rvtemp, 0.0D0, 'REDUCE  '))then

            distance(i) = NCDIST( rvtemp, rv2 )

            distance2 = 0.0
            distancesource = 0.0
            do 7800 k=1,6
               distance2 = distance2 + (rv2(k)-rvtemp(k))**2
C               distancesource = distancesource + (v3(k)-rvtemp(k))**2
 7800       continue
               ibuff = ibuff + 1
               tell = " "
               if ( distance(i)/1.01D0.gt.
     *             dsqrt(distance2) ) then
                  tell = "X"
                  iprintNCDIST_Too_Large = .true.
               endif
               write(buffer(ibuff),"('('i4')',i5,3x,f5.3,3x, 6f9.4,
     *            5x,6f9.4, 1x, 2F9.4,1x,A1)")
     *            ID, i, factor, rvtemp, vtemp,
     *            distance(i),
     *            dsqrt(distance2), tell
               if (xdebug)
     *            write(7,"('('i4')',i5,3x,f5.3,3x, 6f9.4,
     *            5x,6f9.4, 1x, 2F9.4,1x,A1)")
     *            ID, i, factor, rvtemp, vtemp,
     *            distance(i),
     *            dsqrt(distance2), tell
            else
               write(*,*) "bad reduce"
               write(*,*) i,vtemp
            endif
            else
               write(*,*) "bad cell"
            endif

 8000    continue

C compute the average and std dev. of the deltas
         absSum = 0.0D0
         absSumSQ = 0.0D0
         maxdelta = -1.0D0
         maxvalue = distance(1)
         do 8050 k=step1+1,step2
            deltaDist = dabs(distance(k-1)-distance(k))
            maxdelta = dmax1(maxdelta, deltaDist)
            absSum = absSum + deltaDist
            absSumSQ = absSumSQ + deltaDist**2
            maxvalue = dmax1(maxvalue,distance(k))
 8050    continue
         avg = absSum / (step2-step1+1)
         stddev = dsqrt( absSumSq/(step2-step1+1)**2 - avg**2/
     *    (step2-step1+1)**2 )

C look for large deltas
            do 8051 k=step1+1,step2
               deltaDist = abs(distance(k)-distance(k-1))

               if( deltaDist > 25.0 * stddev ) then
                  iprintStepDist = .true.
                  go to 8052
               endif
 8051       continue
 8052       continue

         if ( iprintStepDist .OR. iprintNCDIST_Too_Large ) then
            write(8,"(2h$ ,4f11.6, A)") 
     *         maxvalue, maxdelta, avg, stddev,
     *         " (maxvalue/maxdelta/avgdelta/stddevdelta)"       
            if ( iprintStepDist .and. iprintNCDIST_Too_Large ) then
               dev = 9
            elseif ( iprintNCDIST_Too_Large ) then
               dev = 10
            else
               dev = 11
            endif
            write(dev,"(2h$ ,4f11.6)") maxvalue, maxdelta, avg, stddev
            do 8100 k=1,ibuff
 8100           write(dev,*) buffer(k)

            write(3,
     *       "('('i4')'3x,6f9.4,5x,6f9.4, 2000F9.4)") ID, v3, rv2,
     *              (distance(j),j=step1,step2)
         endif
         call FLUSH( )
C      endif
      endif

 9000 continue
      end


C**********************************************************************C
      integer function initCleanupList( nc, c )
C make a list of a lot of rational numbers
C this will be used to "clean up" the iterated projectors to correct
C for floating point approximations (if any)
      implicit none

      integer nc, i5000, i, icount
      real*8 c(nc)
      integer maximumDenominator
      parameter (maximumDenominator=30)
C----------------------------------------------------------------------C
      icount = 1
      c(icount) = 0.0D0

      do 5000 i5000=1,maximumDenominator-1
         do 1000 i=i5000+1,maximumDenominator
            icount = icount + 1
            c(icount) = DBLE(i5000)/DBLE(i)
            icount = icount + 1
            c(icount) = -c(icount)
 1000    continue
 5000 continue

      initCleanupList = icount
      end

C**********************************************************************C
      subroutine cleanupProjector( p, nc, c )
C p is the projector to "clean up"
C nc is the number of rationals to compare to the values in p
C c is the array of previously computed rationals

      implicit none

      integer nc, i5000, i, icount
      real*8 p(36)
      real*8 c(nc)
C----------------------------------------------------------------------C
      do 5000 i5000=1,36
         do 1000 i=1,nc
            if ( abs(p(i5000)-c(i)) .lt. 0.0001D0 ) then
               p(i5000) = c(i)
               go to 4000
            endif
 1000    continue
 4000    continue
 5000 continue
      end

C**********************************************************************C
      subroutine InitStack( nlevels, stack )
C put the initial data into the stack of projector indices
      implicit none

      integer nlevels, i
      integer stack(nlevels)
C----------------------------------------------------------------------C
      do 1000 i=1,nlevels
         stack(i) = i
 1000 continue
      end

C**********************************************************************C
      subroutine InitPrj( nlevels, stack, prjout )
C create the initial product of nlevels projectors (for squaring later)
      implicit none

      integer nlevels, i, n, nprev
      common /mat36/ prj, mat, rslts
      real*8 mat(36,18)
      real*8 prj(36,18)
      character*1 rslts(155)
      real*8 prjout(36)
      integer stack(nlevels)
      real*8 unitmx(36)
      real*8 m1(36),m2(36),m3(36)
C----------------------------------------------------------------------C

      do 1000 i=1,36
 1000 m1(i) = 0.0D0
      do 1010 i=1,36,7
 1010 m1(i) = 1.0D0

      nprev = -1
      do 2000 i=1,nlevels
          n = stack(i)
          if (n.eq.7 .and. nprev.eq.6) n = 16
          if (n.eq.10 .and. nprev.eq.9) n = 17
          if (n.eq.13 .and. nprev.eq.12) n = 18
          call mm6( m1, prj(1,n), m2 )
          call cpyvn( 36, m2, m1 )
          nprev = n
 2000 continue

      call cpyvn( 36, m1, prjout )
      end

C**********************************************************************C
      real*8 FUNCTION DOTVN (N,V1,V2)
      implicit none

      integer i, N
      real*8 V1(N),V2(N)
C-----------------------------------------------------------------------
      DOTVN = 0.0D0
      DO 1000 I=1,N
         DOTVN = DOTVN + V1(I)*V2(I)
 1000 CONTINUE
      END

