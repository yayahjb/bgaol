#include <cmath>
#include <algorithm>

#include "BasicDistance.h"
#include "MakeGaol.h"
#include "Reducer.h"
#include "TNear.h"

/*  class MakeGaol
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A class to implement an embedding the G6 representation of Niggli space,
following the design of Andrews and Bernstein, 2012

   MakeGaol(void) == constructor, really does nothing --- for the present you need MakeGaolEntry to complete setup
   ~MakeGaol(void);

   void MakeGaolEntry( arma::vec6 v[], const arma::vec6& gred, const arma::vec6& g6ErrorBox, const double ratio,
      int& nv, double vdist[], int ivb[] );

   double NCDist( const arma::vec6& gv1, const arma::vec6& gv2 )
            == returns the distance in folded G6 Niggli space between two points
   static bool inncone( const arma::vec6& gvec )
            == returns whether a G6 point is in the Niggli cone
   void BDCoord( const arma::vec6& gvec, arma::vec6& xs, arma::vec6& ys )
            == returns the distances of gvec to related boundary "groups" -- this 
               function should probably not use vec6 as the return type
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/


const bool debug = true;
const double SMALLVDIST ( 0.0001 );
const double SMALLNEARREDDELTA( 1.0e-6 );

MakeGaol::MakeGaol(void)
   : m_prjList( )
{
}

MakeGaol::~MakeGaol(void)
{
}

double ms[36];


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



void MakeGaol::MakeGaolEntry( arma::vec6 v[], const arma::vec6& gred, 
   const arma::vec6& g6ErrorBox, const double ratio, int& nv, double vdist[], int ivb[] )
{
   bool goodseed[16];
   arma::vec6 seeds[16];
   double seeddist[16];
   arma::vec6 vtred;
   arma::vec6 newge;
   double dred;
   const int nvmax = 100;
   CNearTree<arma::vec6> tree;

   //     Make the perpendicular projectors
   //     Each of the PERP projectors takes a g6
   //     vector to the space othgonal to the
   //     corresponding PRJ.  The length of that
   //     vector is therefore the distance to that
   //     polytope
   //     Generate seed points on each of the 15 boundaries

   GenSeeds( gred, g6ErrorBox, seeds, seeddist, goodseed ); //call GENSEEDS(GRED,GE,SEEDS,SEEDDIST,

   tree.clear( );

   double gerr = ratio * arma::norm( gred, 2 );
   const double dmin = std::min( 0.05, 0.1*gerr );
   const double grmin = gred[0]*gred[0] + gred[1]*gred[1] + gred[2]*gred[2]; //  grmin = gred(1)**2+gred(2)**2+gred(3)**2

   if ( debug && dmin>999999999.9989 )
   {
      printf( "DMIN in MakeGaolEntry %g", dmin );
   }
   else
   {
      printf( "DMIN in MakeGaolEntry %f", dmin );
   }

   //     Populate the tree with the original cell and
   //     the nearby neighbors on the boundary

   nv = 0; //   ++nv; // start from zero for C-style indexing

   tree.insert( gred );
   v[nv] = gred;
   int nexamined = nv;
   ivb[nv] = 0;
   vdist[nv] = 0.0;
   //
   //     Add the good seeds that are not duplicates to
   //     the tree, and for each one apply the matrix
   //     for that boundary and if the image is nearly
   //     reduced, add that as well.
   //
   int nv1 = 0;
   //!!!!!!!!!!!!!!!!!!!!!! are the seeds in one-to-one with the boundaries - do we need to dimension them 16?
   for ( int i=1; i<=15; ++i ) //  do i = 1,15   
   {
      if( goodseed[i] )  //if (GOODSEED(i)) then
      {
         nv1 = nv + 1;
         arma::vec6 nearest;
         if ( ! tree.NearestNeighbor( dmin, nearest,  seeds[i] ) ) // IF (NEARST_N(.FALSE.,MXTREE,6,SEEDS(1,i),DMIN,TREE,NV1,ID,
         {
            ++nv; //   NV = NV + 1
            tree.insert( seeds[i] ); //  CALL BLDTRE_N (.FALSE., MXTREE,6, SEEDS(1,i),
            v[nv] = seeds[i];
            ivb[nv] = 1;
            vdist[nv] = seeddist[i];
            if ( vdist[nv] < SMALLVDIST ) vdist[nv] = 0.0;
         }

         const arma::vec6 vt = ms[i] * seeds[i]; //  call imv6(SEEDS(1,i),MS(1,i),vt)
         nv1 = nv + 1;
         int nred = Reducer::NearRed( vt, SMALLNEARREDDELTA );
         if ( nred > 0 )
         {
            const double disttemp = NCDist( vt, gred ); //  VDIST(NV1) = NCDIST(vt,GRED)

            arma::vec6 nearest1;
            if (disttemp<gerr && ! tree.NearestNeighbor( dmin, nearest1, vt ) ) //   if (VDIST(NV1).LT.gerr .and.
            {
               ++nv;
               tree.insert( vt ); //   CALL BLDTRE_N (
               v[nv] = vt;
               ivb[nv] = 1;
               vdist[nv] = disttemp;
               if ( vdist[nv]< SMALLVDIST ) vdist[nv] = 0.0;
            }
         }
         else // ELSE for if nred
         {
            nred = Reducer::Near2Red( vt, gerr, vtred, dred ); //  NRED=NEAR2RED(
            const double disttemp = NCDist( vtred, gred ); // VDIST(NV1) = NCDIST(vtred,GRED)
            if (disttemp<gerr && ! tree.NearestNeighbor( dmin, nearest, vtred ) )
            {
               ++nv; //   NV = NV+1
               tree.insert( vtred ); //  CALL BLDTRE_N(.FALSE.,MXTREE,6,vtred,
               v[nv] = vtred; // CALL CPYVN(6,vtred,V(1,NV))
               ivb[nv] = 1; // IVB(NV) = 1
               vdist[nv] = disttemp;
               if ( vdist[nv] < SMALLVDIST ) vdist[nv] = 0.0; // if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
            } // ENDIF
         } // ENDIF // end if nred
      } // end if GoodSeed
   } // end for

   for ( int i=0; i<6; ++i )  // DO I = 1,6
      newge[i] = std::max( g6ErrorBox[i], dmin ); //   NEWGE(I) = max(GE(I),DMIN)
   // enddo

   //     Now we have seeds to examine from nexamined+1
   //     through NV.  In the course of doing so, NV may
   //     increase.

   while ( true ) // 1000 nexamined = nexamined+1
   {
      ++nexamined;
      if ( nexamined > nv || nv > nvmax-31 )
      {
         printf(" NV in MKREFL %d \n", nv );
         nv = tree.size();
         return;
      } //  ENDIF

      if ( ivb[nexamined] < 10 ) continue; // go to 1000

      for ( int i=0; i<6; ++i )
      {
         newge[i] = std::max( pow( g6ErrorBox[i]/1.414, ivb[nexamined] ),dmin );
      }  // ENDDO

      gerr = arma::norm( newge, 2 );

      //     Generate seed points on each of the 15 boundaries
      GenSeeds( v[nexamined], newge, seeds, seeddist, goodseed );
      // GENSEEDS(V(1,nexamined),NEWGE,SEEDS,SEEDDIST,
      // GOODSEED, PRJ, PRJPERP)
      //
      //     Add the good seeds that are not duplicates to
      //     the tree, and for each one apply the matrix
      //     for that boundary and if the image is nearly
      //     reduced, add that as well.
      //
      //     Differs from the prior loop only in the distance
      //     calculation
      //
      for ( int i=1; i<=15; ++i ) // do i = 1,15
      {
         const double d = seeds[i][0]*seeds[i][0] + seeds[i][1]*seeds[i][1] + seeds[i][2]*seeds[i][2];
         if ( goodseed[i] && d<10.8*grmin ) //if (GOODSEED(i) .and.
         {
            nv1 = nv + 1; // NV1 = NV + 1
            arma::vec6 closest;
            const double dminTemp = dmin * pow(1.85,ivb[nexamined]);
            if ( ! tree.NearestNeighbor( dminTemp, closest, seeds[i] ) ) //IF (NEARST_N(.false.,MXTREE,6,SEEDS(1,i),DMIN
               //*       *(1.85D0**IVB(nexamined)),TREE,NV1,ID,
            { //*      'NEARST_N') .EQ. 0) THEN
               ++nv; // NV = NV + 1
               tree.insert( seeds[i] ); //CALL BLDTRE_N (.FALSE., MXTREE,6, SEEDS(1,i),
               v[nv] = seeds[i]; // CALL CPYVN(6,SEEDS(1,i),V(1,NV))
               ivb[nv] = ivb[nexamined] + 1; // IVB(NV) = IVB(nexamined)+1
               vdist[nv] = sqrt( seeddist[i]*seeddist[i] + vdist[nexamined]*vdist[nexamined] ); // VDIST(NV) = dsqrt(SEEDDIST(i)**2+
               //  *        VDIST(nexamined)**2)
               if ( vdist[nv] < SMALLVDIST ) vdist[nv] = 0.0; // if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
               //            write(*,*)"NV, SEED DIST,IB",NV,VDIST(NV),IB
               //            call printg6('Store SEED ',SEEDS(1,i))
            } //ENDIF
            const arma::vec6 vt = ms[i] * seeds[i]; // call imv6(SEEDS(1,i),MS(1,i),vt)
            nv1 = nv + 1; // NV1 = NV+1
            const bool nredtemp = Reducer::NearRed( vt, SMALLNEARREDDELTA ); // NRED=NEARRED(vt,1.0D-6,'NEARRED')
            if ( nredtemp ) //if (NRED) THEN
            {
               vdist[nv1] = NCDist( vt, gred ); // VDIST(NV1) = NCDIST(vt,GRED)
               if ( vdist[nv1]<gerr && ! tree.NearestNeighbor( dmin*pow(1.85,ivb[nexamined]), closest, vt ) )
               {
                  //     if (VDIST(NV1).LT.gerr .and.
                  //*       NEARST_N(.false.,MXTREE,6,VT,DMIN
                  //*       *(1.85D0**IVB(nexamined)),TREE,NV1,ID,
                  //*      'NEARST_N') .EQ. 0) THEN
                  ++nv; //NV = NV+1
                  tree.insert( vt ); //CALL BLDTRE_N (.FALSE.,MXTREE,6,VT,//*          NV*16+i,TREE,'BLDTRE_N')
                  v[nv] = vt; // CALL CPYVN(6,VT,V(1,NV))
                  ivb[nv] = ivb[nexamined] + 1; // IVB(NV) = IVB(nexamined)+1
                  if ( vdist[nv]<SMALLVDIST ) vdist[nv] = 0.0; //  if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
                  //            write(*,*)"NV, M*SEED DIST,IB",NV,VDIST(NV),IB
                  //            call printg6('Store M*SEED ',VT)
               } //ENDIF
            }
            else //ELSE
            {
               const bool nred = Reducer::Near2Red( vt, gerr, vtred, dred ); //NRED=NEAR2RED(vt,gerr,vtred,dred,'NEAR2RED')
               vdist[nv1] = NCDist( vtred, gred ); //VDIST(NV1) = NCDIST(vtred,GRED)
               if ( vdist[nv1]<gerr && ! tree.NearestNeighbor( dmin*pow(1.85,ivb[nexamined]), closest, vtred ) ) //if (VDIST(NV1).LT.gerr .and.
                  //       NEARST_N(.false.,MXTREE,6,vtred,DMIN
                  //       *(1.85D0**IVB(nexamined)),TREE,NV1,ID,
                  //      'NEARST_N') .EQ. 0) THEN
               {
                  ++nv; //NV = NV+1
                  tree.insert( vtred ); // CALL BLDTRE_N (.FALSE.,MXTREE,6,vtred,
                  //           NV*16+i,TREE,'BLDTRE_N')
                  v[nv] = vtred; //CALL CPYVN(6,vtred,V(1,NV))
                  ivb[nv] = 1; //IVB(NV) = 1
                  if ( vdist[nv]<SMALLVDIST ) vdist[nv] = 0.0; //  if (VDIST(NV).LT.1.D-4) VDIST(NV) = 0.D0
                  //            write(*,*)"NV, M*SEED ** DIST",NV,VDIST(NV)
                  //            call printg6('Store M*SEED ** ',VT)
               } // ENDIF
            } //ENDIF
            //        else
            //           call printg6('Reject SEED ',SEEDS(1,i))
         } //endif
      } // enddo  do i=1,15
   } //go to 1000 // end while 1000
   nv = tree.size();

}


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


      arma::mat66 MakeGaol::MKPerp( const arma::mat66& prj ) // SUBROUTINE MKPERP( PRJ, PRJPERP)
      {
         return( arma::mat66(arma::eye(6,6) - prj) );
      }


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


   void MakeGaol::PrintG6 ( const std::string& text, const arma::vec6& v ) //      SUBROUTINE PRINTG6(TEXT,VECTOR)
   {
      printf( "%s  %f %f %f %f %f %f", text.c_str(), v[0],v[1],v[2],v[3],v[4],v[5] ); //      WRITE(*,'(A,1x,6F9.2)')TEXT,VECTOR
   }

   void MakeGaol::PrintM66_INT ( const std::string& text, const arma::mat66& v ) //      SUBROUTINE PRINTG6(TEXT,VECTOR)
   {
      for ( int i=0; i<6; ++i )
      {
         printf( "%s ", text.c_str() );
         for ( int j=0; j<36; j+=6 )
         {
            printf( " %3d", int(v[i+j]) );
         }
         printf( "\n");
      }
   }

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//-----------------------------------------------------------------------------
// Name: GenSeeds()
//C     Generate G6 seed vectors for boundary definition
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void MakeGaol::GenSeeds( const arma::vec6& gvec, const arma::vec6& g6ErrorBox, arma::vec6 seeds[], double seeddist[], bool goodseed[] ) const
{
//      SUBROUTINE GENSEEDS( GVEC, GE, SEEDS, SEEDDIST,
//     * GOODSEED, PRJ, PRJPERP)
//C
//C     GENSEEDS TAKES AN ARBITRARY G6 VECTOR AND PROJECTS
//C     IT ONTO EACH OF THE 15 BOUNDARIES.  EVEN IF THE
//C     ORIGINAL VECTOR WAS REDUCED THE PROJECTED VECTOR
//C     MAY FAIL TO BE REDUCED, AND MAY EVEN FAIL TO BE
//C     NEARLY REDUCED.  IT WILL BE MARKED AS A GOOD SEED
//C     PROVIDED IT IS BOTH NEARLY REDUCED AND IS WITHIN
//C     3.5 TIMES THE ERROR BOX GE OF GVEC.
//C
//      implicit none
//      real*8 GVEC(6), GE(6)
//      real*8 PRJ(36,25),PRJPERP(36,25)
//      real*8 SEEDS(6,15), SEEDDIST(15)
//      logical GOODSEED(15)
//C      character*20 SEEDLAB
//      real*8 vtemp(6)
//      real*8 boundary67
//      real*8 boundary9A
//      real*8 boundaryCD
//      real*8 a,dred
//      integer i, ip
//      real*8 anorm, gerr
//      logical NEARRED,NEAR2RED
//C----------------------------------------------------------------------C
//
//C     Compute the seeds and the seed distances
//
//C      write(*,*) "GENSEEDS GVEC, GE"
//C      write(*,*) GVEC
//C      write(*,*) GE

   double gerr = arma::norm( g6ErrorBox, 2 ); //      gerr= anorm(6,GE)
   double a;
   //
   arma::vec6 vtemp = m_prjList[19].GetPerp() * gvec; //      call rmv6( GVEC, PRJPERP(1,19), vtemp )
   const double boundary67 = arma::norm( vtemp,2 ); //      boundary67 = anorm(6,vtemp)
   vtemp = m_prjList[20].GetPerp() * gvec; //      call rmv6( GVEC, PRJPERP(1,20), vtemp )
   const double boundary9A = arma::norm( vtemp, 2 ); //      boundary9A = anorm(6,vtemp)
   vtemp = m_prjList[21].GetPerp() * gvec; //      call rmv6( GVEC, PRJPERP(1,21), vtemp )
   const double boundaryCD = arma::norm( vtemp, 2 ); //      boundaryCD = anorm(6,vtemp)
   //
   //C      call printg6('GENSEEDS GVEC ',GVEC)
   //
   for ( int ip=1; ip<=15; ++ip )//      do ip=1,15
   {
      seeds[ip] = m_prjList[ip].GetProjector() * gvec; //        call rmv6( GVEC, PRJ(1,ip), SEEDS(1,ip))
      vtemp = m_prjList[ip].GetPerp() * gvec; //        call rmv6( GVEC, PRJPERP(1,ip), vtemp )
      //C        write(*,'(A,i2/,6(6F9.3/))'),"PRJ",
      //C     *  ip,(PRJ(i,ip),i=1,36)
      //C        write(*,'(A,i2/,6(6F9.3/))'),"PRJPERP",
      //C     *  ip,(PRJPERP(i,ip),i=1,36)
      a = arma::norm( vtemp, 2 ); //        a = anorm(6,vtemp)
      const double testg4 = gvec[3]-1.0E-6*sqrt(gvec[1]*gvec[2]);
      const double testg5 = gvec[4]-1.0E-6*sqrt(gvec[0]*gvec[2]);
      const double testg6 = gvec[5]-1.0E-6*sqrt(gvec[0]*gvec[2]);
      const bool ipTestA = ip==6 || ip==7 || ip==9 || ip==10 || ip==12 || ip==13;
      const bool ipTestB = ip==8 || ip==11 || ip==14 || ip==15;

      if ( ( testg4*testg5*testg6 <= 0.0 && ipTestA ) ||
           ( testg4*testg5*testg6 >  0.0 && ipTestB ) )
         //        if (((GVEC(4)-1.d-6*sqrt(GVEC(2)*GVEC(3)))*
         //     *     (GVEC(5)-1.d-6*sqrt(GVEC(1)*GVEC(3)))*
         //     *     (GVEC(6)-1.d-6*sqrt(GVEC(1)*GVEC(3)))
         //     *     .le. 0.D0
         //     *        .and. (ip .eq. 6 .or. ip .eq. 7
         //     *          .or. ip .eq. 9 .or. ip .eq. 10
         //     *          .or. ip .eq. 12.or. ip .eq. 13))
         //     *        .or. (
         //     *     (GVEC(4)-1.d-6*sqrt(GVEC(2)*GVEC(3)))*
         //     *     (GVEC(5)-1.d-6*sqrt(GVEC(1)*GVEC(3)))*
         //     *     (GVEC(6)-1.d-6*sqrt(GVEC(1)*GVEC(3)))
         //     *     .gt. 0.D0
         //     *        .and. (ip .eq. 8 .or. ip .eq. 11
         //     *          .or. ip .eq. 14. or. ip .eq. 15))) then
      {
         if ( ip>=6 && ip<=8 )   //            if ( ip .eq. 6 .or. ip .eq. 7 .or. ip .eq. 8 ) then
         {
            seeds[ip] = m_prjList[16].GetProjector() * gvec; //              call rmv6( GVEC, PRJ(1,16), SEEDS(1,ip) )
            vtemp = m_prjList[16].GetPerp() * gvec; //              call rmv6( GVEC, PRJPERP(1,16), vtemp )
         }   //            endif
         else if ( ip>=9 && ip<=15 )   //            if ( ip .g6ErrorBox. 9 .and. ip .le. 14) then
         {
            seeds[ip] = m_prjList[17].GetProjector() * gvec;             //              call rmv6( GVEC, PRJ(1,17), SEEDS(1,ip) )
            vtemp = m_prjList[17].GetPerp() * gvec;             //              call rmv6( GVEC, PRJPERP(1,17), vtemp )
            //            endif
         }
         else// if ( ip== 15 )  //            if ( ip .eq. 15) then
         {
            seeds[ip] = m_prjList[18].GetProjector() * gvec;             //              call rmv6( GVEC, PRJ(1,18), SEEDS(1,ip) )
            vtemp = m_prjList[18].GetPerp() * gvec;             //              call rmv6( GVEC, PRJPERP(1,18), vtemp )
         }   //            endif
         a = arma::norm( vtemp, 2 ); //            a = anorm(6,vtemp)
      }//        endif

           if ( ip== 6 && gvec[4]<=gvec[5] && gvec[4]>0.0 ) a = boundary67;
      else if ( ip== 7 && gvec[4]>=gvec[5] && gvec[5]>0.0 ) a = boundary67;
      else if ( ip== 9 && gvec[3]<=gvec[5] && gvec[3]>0.0 ) a = boundary9A;
      else if ( ip==10 && gvec[3]>=gvec[5] && gvec[5]>0.0 ) a = boundary9A;
      else if ( ip==12 && gvec[3]<=gvec[4] && gvec[3]>0.0 ) a = boundaryCD;
      else if ( ip==13 && gvec[3]>=gvec[4] && gvec[4]>0.0 ) a = boundaryCD;

      //        if ( ip .eq.  6 .and. GVEC(5) .le. GVEC(6)
      //     *     .and. GVEC(5) .gt. 0.D0) then
      //           a = boundary67
      //        elseif ( ip .eq.  7 .and. GVEC(5) .g6ErrorBox. GVEC(6)
      //     *     .and. GVEC(6) .gt. 0.D0) then
      //           a = boundary67
      //        elseif ( ip .eq.  9 .and. GVEC(4) .le. GVEC(6)
      //     *     .and. GVEC(4) .gt. 0.D0 ) then
      //           a = boundary9A
      //        elseif ( ip .eq. 10 .and. GVEC(4) .g6ErrorBox. GVEC(6)
      //     *     .and. GVEC(6) .gt. 0.D0 ) then
      //           a = boundary9A
      //        elseif ( ip .eq. 12 .and. GVEC(4) .le. GVEC(5)
      //     *     .and. GVEC(4) .gt. 0.D0 ) then
      //           a = boundaryCD
      //        elseif ( ip .eq. 13 .and. GVEC(4) .g6ErrorBox. GVEC(5)
      //     *     .and. GVEC(5) .gt. 0.D0 ) then
      //           a = boundaryCD
      //        endif
      //
      //C       get distance and mark as bad if outside
      //C       3.5 times the errorbox
      //
      seeddist[ip] = a;//        SEEDDIST(ip) = a
      goodseed[ip] = Reducer::NearRed( seeds[ip], SMALLNEARREDDELTA ); //        GOODSEED(ip) = NEARRED(SEEDS(1,ip),1.0D-6,'NEARRED')
      if ( ! goodseed[ip] ) //        if (.NOT.GOODSEED(ip)) then
      {
         double dred;
         goodseed[ip] = Reducer::Near2Red( seeds[ip], gerr, vtemp, dred );
         seeds[ip] = vtemp;
         //          GOODSEED(ip) = NEAR2RED(SEEDS(1,ip),gerr,
         //     *    vtemp,dred,'NEAR2RED')
         //          call CPYVN(6,vtemp,SEEDS(1,ip))
         //C         call printg6("Non-reduced projection",seeds(1,ip))
         //
      }//        endif
      for ( int i=0; i<6; ++i )//        do i = 1,6
      {
         goodseed[ip] = goodseed[ip] && abs(gvec[i])-seeds[ip][i] <= 3.5*g6ErrorBox[i];
         //          if (abs(GVEC(i)-SEEDS(i,ip)) .gt.
         //     *      3.5D0*GE(i)) GOODSEED(ip) = .false.
      }//        enddo
      //C        write(seedlab,'(''BDRY'', I2)')IP
      //C        call printg6(seedlab,seeds(1,ip))
   } // loop over 15 boundaries //      enddo

   //return
} //   end

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


//-----------------------------------------------------------------------------
// Name: BDCoord()
//C     Compute XS (distances to boundary sets) and
//C             YS (signed distances along boundaries)
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void MakeGaol::BDCoord( const arma::vec6& gvec, arma::vec6& xs, arma::vec6& ys )//      SUBROUTINE BDCOORD(GVEC,XS,YS)
{
   //      implicit none
   //      real*8 XS(6),YS(6)
   //      real*8 GVEC(6)
   //
   //     The boundaries are in the order
   //       1, 2, 678, 9AB, CDE, F
   const double sqrt2( sqrt(2.0) ); //      real*8 sqrt2, sqrt5
   const double sqrt5( sqrt(5.0) );
   xs[0] = (gvec[1] -     gvec[0]) /sqrt2;
   xs[1] = (gvec[2] -     gvec[1]) /sqrt2;
   xs[2] = (gvec[1] - abs(gvec[3]))/sqrt2;
   xs[3] = (gvec[0] - abs(gvec[4]))/sqrt2;
   xs[4] = (gvec[0] - abs(gvec[5]))/sqrt2;
   xs[5] = (gvec[0]+gvec[1] + gvec[3] + gvec[4] + gvec[5])/sqrt5;

   ys[0] = (abs(gvec[3]) - abs(gvec[4]))/sqrt2; // YS(1) = (ABS(GVEC(4)) - ABS(GVEC(5)))/sqrt2
   ys[1] = (abs(gvec[4]) - abs(gvec[5]))/sqrt2; // YS(2) = (ABS(GVEC(5)) - ABS(GVEC(6)))/sqrt2
   ys[2] = (    gvec[4]  - abs(gvec[5])/2.0)/sqrt2;
   ys[3] = (    gvec[3]  - abs(gvec[5])/2.0)/sqrt2;
   ys[4] = (    gvec[3]  - abs(gvec[4])/2.0)/sqrt2;
   ys[5] = (    gvec[1]  - gvec[0] + gvec[3] - gvec[4]) / 2.0;
} //      END


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



//-----------------------------------------------------------------------------
// Name: MKGReporter()
// Description: prints distance vector information for debugging, etc.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void MKGReporter( const std::string& str, const double ncdist, const arma::vec6& gv1, const arma::vec6& gv2 )
{
   //printf( "%s %g \n v1= %f %f %f %f %f %f \n v2= %f %f %f %f %f %f \n" , str.c_str(), ncdist, 
   //   gv1[0], gv1[1], gv1[2], gv1[3], gv1[4], gv1[5], gv2[0], gv2[1], gv2[2], gv2[3], gv2[4], gv2[5] );
}

//-----------------------------------------------------------------------------
// Name: NCDist()
// Description: Determines the minimal distance between to G6 vectors with folding
//              of the path to reside within the Niggli Cone (NC). This
//              function requires that both input vectors be for reduced cells.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double MakeGaol::NCDist( const arma::vec6& gv1, const arma::vec6& gv2 )     // REAL*8 FUNCTION NCDIST(gv1,gv2)
{
   // implicit none
   // real*8 gv1(6), gv2(6)
   // real*8 FOLDMDIST
   // real*8 g456dist
   // common/xdebug/xdebug
   // logical xdebug
   // real*8 prj(36,25), prjperp(36,25)
   // real*8 prjhat(36,15)
   // integer ms(36,18)
   // common /projectors/ prj,prjperp,prjhat,ms

   // if (xdebug) then
   //   write(7,'(a,1x,6f9.2)') "NCDIST gv1 ",gv1
   //   write(7,'(a,1x,6f9.2)') "NCDIST gv2 ",gv2
   // endif
   double ncdist( 9.0E38);
   double ncdistPrevious(ncdist);
   std::string str;
   ncdist = FoldMDist( gv1, gv2, 0, ncdist );      // NCDIST=foldmdist(gv1,gv2,0,9.D38)
   str = std::string("ncdist1"); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 7, ncdist );     // NCDIST=foldmdist(gv1,gv2,7,NCDIST)
   str = std::string("ncdist2") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 10, ncdist );     // NCDIST=foldmdist(gv1,gv2,10,NCDIST)
   str = std::string("ncdist3") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 13, ncdist );     // NCDIST=foldmdist(gv1,gv2,13,NCDIST)
   str = std::string("ncdist4") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 15, ncdist );     // NCDIST=foldmdist(gv1,gv2,15,NCDIST)
   str = std::string("ncdist5") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 23, ncdist );     // NCDIST=foldmdist(gv1,gv2,23,NCDIST)
   str = std::string("ncdist6") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 24, ncdist );     // NCDIST=foldmdist(gv1,gv2,24,NCDIST)
   str = std::string("ncdist7") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 25, ncdist );     // NCDIST=foldmdist(gv1,gv2,25,NCDIST)
   str = std::string("ncdist8") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif





   //ncdist = FoldMDist( gv2, gv1, 0, ncdist );      // NCDIST=foldmdist(gv1,gv2,0,9.D38)
   //str = std::string("ncdist1a") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   //MKGReporter( str, ncdist, gv2, gv1 );
   //// if (xdebug) then
   //// write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   ////* NCDIST
   //// endif
   //ncdist = FoldMDist( gv2, gv1, 7, ncdist );     // NCDIST=foldmdist(gv1,gv2,7,NCDIST)
   //str = std::string("ncdist2a") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   //MKGReporter( str, ncdist, gv2, gv1 );
   //// if (xdebug) then
   //// write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   ////* NCDIST
   //// endif
   //ncdist = FoldMDist( gv2, gv1, 10, ncdist );     // NCDIST=foldmdist(gv1,gv2,10,NCDIST)
   //str = std::string("ncdist3a") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   //MKGReporter( str, ncdist, gv2, gv1 );
   //// if (xdebug) then
   //// write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   ////* NCDIST
   //// endif
   //ncdist = FoldMDist( gv2, gv1, 13, ncdist );     // NCDIST=foldmdist(gv1,gv2,13,NCDIST)
   //str = std::string("ncdist4a") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   //MKGReporter( str, ncdist, gv2, gv1 );
   //// if (xdebug) then
   //// write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   ////* NCDIST
   //// endif
   //ncdist = FoldMDist( gv2, gv1, 15, ncdist );     // NCDIST=foldmdist(gv1,gv2,15,NCDIST)
   //str = std::string("ncdist5a") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   //MKGReporter( str, ncdist, gv2, gv1 );
   //// if (xdebug) then
   //// write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   ////* NCDIST
   //// endif
   //ncdist = FoldMDist( gv2, gv1, 23, ncdist );     // NCDIST=foldmdist(gv1,gv2,23,NCDIST)
   //str = std::string("ncdist6a") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   //MKGReporter( str, ncdist, gv2, gv1 );
   //// if (xdebug) then
   //// write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   ////* NCDIST
   //// endif
   //ncdist = FoldMDist( gv2, gv1, 24, ncdist );     // NCDIST=foldmdist(gv1,gv2,24,NCDIST)
   //str = std::string("ncdist7a") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   //MKGReporter( str, ncdist, gv2, gv1 );
   //// if (xdebug) then
   //// write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   ////* NCDIST
   //// endif
   //ncdist = FoldMDist( gv2, gv1, 25, ncdist );     // NCDIST=foldmdist(gv1,gv2,25,NCDIST)
   //str = std::string("ncdist8a") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   //MKGReporter( str, ncdist, gv2, gv1 );
   //// if (xdebug) then
   //// write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   ////* NCDIST
   //// endif



   return( ncdist );
}     // end

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



//-----------------------------------------------------------------------------
// Name: DistF()
//C     Compute the distance to the F boundary
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double MakeGaol::DistF( const arma::vec6& g ) const //      real*8 function DistF(g)
{
//      implicit none
//      real*8 g(6), x12
//      common/xdebug/xdebug
//      logical xdebug
//
//      x12 = dabs(g(1))+dabs(g(2))+dabs(g(3))
//     * -max(dabs(g(1)),dabs(g(2)),dabs(g(3)))
   const double x12( abs(g[0])+abs(g[1])+abs(g[2]) - std::max(abs(g[0]),std::max(abs(g[1]),abs(g[2]) ) ) );

   const double min1 = std::min( abs(x12+g[3]+g[4]+g[5]), //      DistF = min(dabs(x12+g(4)+g(5)+g(6)),
                                 abs(x12+g[3]-g[4]-g[5]));
   const double min2 = std::min( abs(x12-g[3]+g[4]-g[5]), //     *            dabs(x12-g(4)+g(5)-g(6)),
                                 abs(x12-g[3]-g[4]+g[5]));

   const double distftemp = std::min( min1, min2 ) / sqrt(5.0);//     *            dabs(x12+g(4)-g(5)-g(6)),
                      //     *            dabs(x12-g(4)-g(5)+g(6)))
//     *            /dsqrt(5.D0)
//       if (xdebug) write(7,*) "Distf ",DistF,g
      return( distftemp );
   }//      end


// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


//-----------------------------------------------------------------------------
// Name: inncone()
// Description: Determines whether a vector is within the Niggli cone
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool MakeGaol::inncone( const arma::vec6& gvec )    // logical function inncone(gvec)
{
     // implicit none
     // real*8 gvec(6),glow,gmid,ghigh,gmax
   int imax;  // integer ii,imax
   const double precn( 1.005 );  // data precn/1.005D0/
   bool innconetemp( false );
   const double glow( std::min( gvec[0], std::min(gvec[1], gvec[2]) ) );  // glow = min(gvec(1), gvec(2), gvec(3))
   if ( glow <= 0.0 ) return( false );
     // if (glow.le.0.D0) return
   const double ghigh( std::max( gvec[0], std::max(gvec[1], gvec[2]) ) );     // ghigh = max(gvec(1), gvec(2), gvec(3))
   const double gmid( gvec[0]+gvec[1]+gvec[2] -glow -ghigh );     // gmid = gvec(1)+gvec(2)+gvec(3)-glow-ghigh
   if ( std::max( abs(gvec[3]), std::max(abs(gvec[4]),abs(gvec[5]))) <= gmid*precn )  // if (max(dabs(gvec(4)),dabs(gvec(5)),
   {
     //*  dabs(gvec(6))).gt.gmid*precn) return
     double gmax = abs(gvec[3]); // gmax = dabs(gvec(4))
     imax = 3;
     if (abs(gvec[4]) > gmax) 
     {
        gmax = abs(gvec[4]);
        imax = 4;
     }
     if (abs(gvec[5]) > gmax) 
     {
        gmax = abs(gvec[5]);
        imax = 5;
     }

     if ( imax!=3 && abs(gvec[3])>glow*precn ) return( false );
     if ( imax!=4 && abs(gvec[4])>glow*precn ) return( false );
     if ( imax!=5 && abs(gvec[5])>glow*precn ) return( false );
     innconetemp = true; // inncone = .true.
   }
   return( innconetemp );
} //   end



// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


//C     FOLDMDIST(GVEC1,GVEC2,IP,CFOLDDIST)
//C       return the distance betweem two G6 vectors
//C       going through the given boundary
//C       (7, A, C, F) after taking the second
//C       vector through the 6 exchanges of A, B, C
//C

//-----------------------------------------------------------------------------
// Name: FoldMDist()
//       return the distance betweem two G6 vectors
//       going through the given boundary
//       (7, A, C, F) after taking the second
//       vector through the 6 exchanges of A, B, C
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double MakeGaol::FoldMDist( const arma::vec6& gvec1, const arma::vec6& gvec2, const int ip, const double cFoldDist ) const
{
//      real*8 function FOLDMDIST(GVEC1,GVEC2,IP,CFOLDDIST)
//      implicit none
//
//      REAL*8 GVEC1(6), GVEC2(6), CFOLDDIST
//      REAL*8 FoldxDist
//      integer ip
//      REAL*8 m1g1(6),m2g1(6),m1m2g1(6)
//      REAL*8 m2m1g1(6),m2m1m2g1(6)
//      REAL*8 m1g2(6),m2g2(6),m1m2g2(6)
//      REAL*8 m2m1g2(6),m2m1m2g2(6)
//      real*8 prj(36,25), prjperp(36,25)
//      real*8 prjhat(36,15)
//      integer ms(36,18)
//      common /projectors/ prj,prjperp,prjhat,ms
//      common/xdebug/xdebug
//      logical xdebug
//      real*8 h1,h1m1,h1m2,h1m1m2,h1m2m1,h1m2m1m2
//      real*8 h2,h2m1,h2m2,h2m1m2,h2m2m1,h2m2m1m2
//      real*8 distF
//      real*8 g123dist,g456dist
//    

   double h1 = DBL_MAX;
   double h2 = DBL_MAX;
   double h1m1 = DBL_MAX;
   double h2m1 = DBL_MAX;
   double h1m2 = DBL_MAX;
   double h2m2 = DBL_MAX;
   double h1m1m2 = DBL_MAX;
   double h2m1m2 = DBL_MAX;
   double h1m2m1 = DBL_MAX;
   double h2m2m1 = DBL_MAX;
   double h1m2m1m2 = DBL_MAX;
   double h2m2m1m2 = DBL_MAX;

   const arma::vec6 m1g1     = m_prjList[ 1].GetTransform( ) * gvec1; //      call imv6(gvec1,MS(1,1),m1g1)
   const arma::vec6 m2g1     = m_prjList[ 2].GetTransform( ) * gvec1; //      call imv6(gvec1,MS(1,2),m2g1)
   const arma::vec6 m1m2g1   = m_prjList[16].GetTransform( ) * gvec1; //      call imv6(gvec1,MS(1,16),m1m2g1)
   const arma::vec6 m2m1g1   = m_prjList[17].GetTransform( ) * gvec1; //      call imv6(gvec1,MS(1,17),m2m1g1)
   const arma::vec6 m2m1m2g1 = m_prjList[18].GetTransform( ) * gvec1; //      call imv6(gvec1,MS(1,18),m2m1m2g1)
   const arma::vec6 m1g2     = m_prjList[ 1].GetTransform( ) * gvec2; //      call imv6(gvec2,MS(1,1),m1g2)
   const arma::vec6 m2g2     = m_prjList[ 2].GetTransform( ) * gvec2; //      call imv6(gvec2,MS(1,2),m2g2)
   const arma::vec6 m1m2g2   = m_prjList[16].GetTransform( ) * gvec2; //      call imv6(gvec2,MS(1,16),m1m2g2)
   const arma::vec6 m2m1g2   = m_prjList[17].GetTransform( ) * gvec2; //      call imv6(gvec2,MS(1,17),m2m1g2)
   const arma::vec6 m2m1m2g2 = m_prjList[18].GetTransform( ) * gvec2; //      call imv6(gvec2,MS(1,18),m2m1m2g2)

   double foldmdist = cFoldDist; //      FOLDMDIST = CFOLDDIST
//      
//      if (ip.eq.0) go to 1000
//      if (ip.ge. 6 .and. ip.le. 8) go to 1007
//      if (ip.ge. 9 .and. ip.le.11) go to 1010
//      if (ip.ge.12 .and. ip.le.14) go to 1013
//      if (ip.eq.15) go to 1015
//      if (ip.eq.23) go to 1023
//      if (ip.eq.24) go to 1024
//      if (ip.eq.25) go to 1025
//      if (ip.eq.22) go to 1022
//
   if ( ip == 0 )// 1000 continue  //HJB this could go before all the above assignments !!!!!!!!!!!!!!!!!!!!!!!!!!!
   {
      foldmdist = std::min( foldmdist, BasicDistance::g123dist( gvec1, gvec2 ) ); //      FOLDMDIST = min(FOLDMDIST,
//     *  g123dist(gvec1,gvec2))
      return( foldmdist ); //       return
   }
//
   else if ( ip>=6 && ip<=8 ) // 1007 continue
   {
      h1       = abs(gvec1[1]    - abs(gvec1[3]))   /sqrt(2.0); //      h1 = abs(gvec1(2)-abs(gvec1(4)))/sqrt(2.D0)
      h2       = abs(gvec2[1]    - abs(gvec2[3]))   /sqrt(2.0); //      h2 = abs(gvec2(2)-abs(gvec2(4)))/sqrt(2.D0)
      h1m1     = abs(m1g1[1]     - abs(m1g1[3]))    /sqrt(2.0); //      h1m1 = abs(m1g1(2)-abs(m1g1(4)))/sqrt(2.D0)
      h2m1     = abs(m1g2[1]     - abs(m1g2[3]))    /sqrt(2.0); //      h2m1 = abs(m1g2(2)-abs(m1g2(4)))/sqrt(2.D0)
      h1m2     = abs(m2g1[1]     - abs(m2g1[3]))    /sqrt(2.0); //      h1m2 = abs(m2g1(2)-abs(m2g1(4)))/sqrt(2.D0)
      h2m2     = abs(m2g2[1]     - abs(m2g2[3]))    /sqrt(2.0); //      h2m2 = abs(m2g2(2)-abs(m2g2(4)))/sqrt(2.D0)
      h1m1m2   = abs(m1m2g1[1]   - abs(m1m2g1[3]))  /sqrt(2.0); //      h1m1m2 = abs(m1m2g1(2)-abs(m1m2g1(4)))/sqrt(2.D0)
      h2m1m2   = abs(m1m2g2[1]   - abs(m1m2g2[3]))  /sqrt(2.0); //      h2m1m2 = abs(m1m2g2(2)-abs(m1m2g2(4)))/sqrt(2.D0)
      h1m2m1   = abs(m2m1g1[1]   - abs(m2m1g1[3]))  /sqrt(2.0); //      h1m2m1 = abs(m2m1g1(2)-abs(m2m1g1(4)))/sqrt(2.D0)
      h2m2m1   = abs(m2m1g2[1]   - abs(m2m1g2[3]))  /sqrt(2.0); //      h2m2m1 = abs(m2m1g2(2)-abs(m2m1g2(4)))/sqrt(2.D0)
      h1m2m1m2 = abs(m2m1m2g1[1] - abs(m2m1m2g1[3]))/sqrt(2.0); //     * abs(m2m1m2g1(2)-abs(m2m1m2g1(4)))/sqrt(2.D0)
      h2m2m1m2 = abs(m2m1m2g2[1] - abs(m2m1m2g2[3]))/sqrt(2.0); //     * abs(m2m1m2g2(2)-abs(m2m1m2g2(4)))/sqrt(2.D0)
//      go to 1100
   }
//
   else if ( ip>=9 && ip<=11 ) // 1010 continue
   {
      h1       = abs(gvec1[0]    - abs(gvec1[4]))   /sqrt(2.0); //      h1 = abs(gvec1(1)-abs(gvec1(5)))/sqrt(2.D0)
      h2       = abs(gvec2[0]    - abs(gvec2[4]))   /sqrt(2.0); //      h2 = abs(gvec2(1)-abs(gvec2(5)))/sqrt(2.D0)
      h1m1     = abs(m1g1[0]     - abs(m1g1[4]))    /sqrt(2.0); //      h1m1 = abs(m1g1(1)-abs(m1g1(5)))/sqrt(2.D0)
      h2m1     = abs(m1g2[0]     - abs(m1g2[4]))    /sqrt(2.0); //      h2m1 = abs(m1g2(1)-abs(m1g2(5)))/sqrt(2.D0)
      h1m2     = abs(m2g1[0]     - abs(m2g1[4]))    /sqrt(2.0); //      h1m2 = abs(m2g1(1)-abs(m2g1(5)))/sqrt(2.D0)
      h2m2     = abs(m2g2[0]     - abs(m2g2[4]))    /sqrt(2.0); //      h2m2 = abs(m2g2(1)-abs(m2g2(5)))/sqrt(2.D0)
      h1m1m2   = abs(m1m2g1[0]   - abs(m1m2g1[4]))  /sqrt(2.0); //      h1m1m2 = abs(m1m2g1(1)-abs(m1m2g1(5)))/sqrt(2.D0)
      h2m1m2   = abs(m1m2g2[0]   - abs(m1m2g2[4]))  /sqrt(2.0); //      h2m1m2 = abs(m1m2g2(1)-abs(m1m2g2(5)))/sqrt(2.D0)
      h1m2m1   = abs(m2m1g1[0]   - abs(m2m1g1[4]))  /sqrt(2.0); //      h1m2m1 = abs(m2m1g1(1)-abs(m2m1g1(5)))/sqrt(2.D0)
      h2m2m1   = abs(m2m1g2[0]   - abs(m2m1g2[4]))  /sqrt(2.0); //      h2m2m1 = abs(m2m1g2(1)-abs(m2m1g2(5)))/sqrt(2.D0)
      h1m2m1m2 = abs(m2m1m2g1[0] - abs(m2m1m2g1[4]))/sqrt(2.0); //     * abs(m2m1m2g1(1)-abs(m2m1m2g1(5)))/sqrt(2.D0)
      h2m2m1m2 = abs(m2m1m2g2[0] - abs(m2m1m2g2[4]))/sqrt(2.0); //     * abs(m2m1m2g2(1)-abs(m2m1m2g2(5)))/sqrt(2.D0)
//      go to 1100
   }
//
   else if ( ip>=12 && ip<=13 )// 1013 continue
   {
      h1       = abs(gvec1[0]   - abs(gvec1[5]   ))/sqrt(2.0); //      h1 = abs(gvec1(1)-abs(gvec1(6)))/sqrt(2.0)
      h2       = abs(gvec2[0]   - abs(gvec2[5]   ))/sqrt(2.0); //      h2 = abs(gvec2(1)-abs(gvec2(6)))/sqrt(2.0)
      h1m1     = abs(m1g1[0]    - abs(m1g1[5]    ))/sqrt(2.0); //      h1m1 = abs(m1g1(1)-abs(m1g1(6)))/sqrt(2.0)
      h2m1     = abs(m1g2[0]    - abs(m1g2[5]    ))/sqrt(2.0); //      h2m1 = abs(m1g2(1)-abs(m1g2(6)))/sqrt(2.0)
      h1m2     = abs(m2g1[0]    - abs(m2g1[5]    ))/sqrt(2.0); //      h1m2 = abs(m2g1(1)-abs(m2g1(6)))/sqrt(2.0)
      h2m2     = abs(m2g2[0]    - abs(m2g2[5]    ))/sqrt(2.0); //      h2m2 = abs(m2g2(1)-abs(m2g2(6)))/sqrt(2.0)
      h1m1m2   = abs(m1m2g1[0]  - abs(m1m2g1[5]  ))/sqrt(2.0); //      h1m1m2 = abs(m1m2g1(1)-abs(m1m2g1(6)))/sqrt(2.0)
      h2m1m2   = abs(m1m2g2[0]  - abs(m1m2g2[5]  ))/sqrt(2.0); //      h2m1m2 = abs(m1m2g2(1)-abs(m1m2g2(6)))/sqrt(2.0)
      h1m2m1   = abs(m2m1g1[0]  - abs(m2m1g1[5]  ))/sqrt(2.0); //      h1m2m1 = abs(m2m1g1(1)-abs(m2m1g1(6)))/sqrt(2.0)
      h2m2m1   = abs(m2m1g2[0]  - abs(m2m1g2[5]  ))/sqrt(2.0); //      h2m2m1 = abs(m2m1g2(1)-abs(m2m1g2(6)))/sqrt(2.0)
      h1m2m1m2 = abs(m2m1m2g1[0]- abs(m2m1m2g1[5]))/sqrt(2.0); //      h1m2m1m2 =
      h2m2m1m2 = abs(m2m1m2g2[0]- abs(m2m1m2g2[5]))/sqrt(2.0); //      h2m2m1m2 =
      //      go to 1100
   }
//
//
   else if ( ip == 15 ) // 1015 continue
   {
      h1       = DistF(gvec1); //      h1 = DistF(gvec1)
      h2       = DistF(gvec2); //      h2 = DistF(gvec2)
      h1m1     = DistF(m1g1); //      h1m1 = DistF(m1g1)
      h2m1     = DistF(m1g2); //      h2m1 = DistF(m1g2)
      h1m2     = DistF(m2g1); //      h1m2 = DistF(m2g1)
      h2m2     = DistF(m2g2); //      h2m2 = DistF(m2g2)
      h1m1m2   = DistF(m1m2g1); //      h1m1m2 = DistF(m1m2g1)
      h2m1m2   = DistF(m1m2g2); //      h2m1m2 = DistF(m1m2g2)
      h1m2m1   = DistF(m2m1g1); //      h1m2m1 = DistF(m2m1g1)
      h2m2m1   = DistF(m2m1g2); //      h2m2m1 = DistF(m2m1g2)
      h1m2m1m2 = DistF(m2m1m2g1); //      h1m2m1m2 = DistF(m2m1m2g1)
      h2m2m1m2 = DistF(m2m1m2g2); //      h2m2m1m2 = DistF(m2m1m2g2)
//      goto 1100
   }
//
//C 8F boundary
   else if ( ip == 23 ) // 1023 continue
   {
      h1       = sqrt(pow(gvec1[1]   -abs(gvec1[3]),2)   /2.0 + pow(gvec1[0]   -abs(gvec1[4])   - abs(gvec1[5]),2)   /3.0);
      h2       = sqrt(pow(gvec2[1]   -abs(gvec2[3]),2)   /2.0 + pow(gvec2[0]   -abs(gvec2[4])   - abs(gvec2[5]),2)   /3.0);
      h1m1     = sqrt(pow(m1g1[1]    -abs(m1g1[3]),2)    /2.0 + pow(m1g1[0]    -abs(m1g1[4])    - abs(m1g1[5]),2)    /3.0);
      h2m1     = sqrt(pow(m1g2[1]    -abs(m1g2[3]),2)    /2.0 + pow(m1g2[0]    -abs(m1g2[4])    - abs(m1g2[5]),2)    /3.0);
      h1m2     = sqrt(pow(m2g1[1]    -abs(m2g1[3]),2)    /2.0 + pow(m2g1[0]    -abs(m2g1[4])    - abs(m2g1[5]),2)    /3.0);
      h2m2     = sqrt(pow(m2g2[1]    -abs(m2g2[3]),2)    /2.0 + pow(m2g2[0]    -abs(m2g2[4])    - abs(m2g2[5]),2)    /3.0);
      h1m1m2   = sqrt(pow(m1m2g1[1]  -abs(m1m2g1[3]),2)  /2.0 + pow(m1m2g1[0]  -abs(m1m2g1[4])  - abs(m1m2g1[5]),2)  /3.0);
      h2m1m2   = sqrt(pow(m1m2g2[1]  -abs(m1m2g2[3]),2)  /2.0 + pow(m1m2g2[0]  -abs(m1m2g2[4])  - abs(m1m2g2[5]),2)  /3.0);
      h1m2m1   = sqrt(pow(m2m1g1[1]  -abs(m2m1g1[3]),2)  /2.0 + pow(m2m1g1[0]  -abs(m2m1g1[4])  - abs(m2m1g1[5]),2)  /3.0);
      h2m2m1   = sqrt(pow(m2m1g2[1]  -abs(m2m1g2[3]),2)  /2.0 + pow(m2m1g2[0]  -abs(m2m1g2[4])  - abs(m2m1g2[5]),2)  /3.0);
      h1m2m1m2 = sqrt(pow(m2m1m2g1[1]-abs(m2m1m2g1[3]),2)/2.0 + pow(m2m1m2g1[0]-abs(m2m1m2g1[4])- abs(m2m1m2g1[5]),2)/3.0);
      h2m2m1m2 = sqrt(pow(m2m1m2g2[1]-abs(m2m1m2g2[3]),2)/2.0 + pow(m2m1m2g2[0]-abs(m2m1m2g2[4])- abs(m2m1m2g2[5]),2)/3.0);
//      go to 1100
   }
//
//C BF boundary
   else if ( ip == 24 ) // 1024 continue
   {
   h1       = sqrt( pow(gvec1[0]   -abs(gvec1[4]),2)   /2.0 + pow(gvec1[1]   -abs(gvec1[3])   -abs(gvec1[5])   ,2) /3.0);
   h2       = sqrt( pow(gvec2[0]   -abs(gvec2[4]),2)   /2.0 + pow(gvec2[1]   -abs(gvec2[3])   -abs(gvec2[5])   ,2) /3.0);
   h1m1     = sqrt( pow(m1g1[0]    -abs(m1g1[4]),2)    /2.0 + pow(m1g1[1]    -abs(m1g1[3])    -abs(m1g1[5])    ,2) /3.0);
   h2m1     = sqrt( pow(m1g2[0]    -abs(m1g2[4]),2)    /2.0 + pow(m1g2[1]    -abs(m1g2[3])    -abs(m1g2[5])    ,2) /3.0);
   h1m2     = sqrt( pow(m2g1[0]    -abs(m2g1[4]),2)    /2.0 + pow(m2g1[1]    -abs(m2g1[3])    -abs(m2g1[5])    ,2) /3.0);
   h2m2     = sqrt( pow(m2g2[0]    -abs(m2g2[4]),2)    /2.0 + pow(m2g2[1]    -abs(m2g2[3])    -abs(m2g2[5])    ,2) /3.0);
   h1m1m2   = sqrt( pow(m1m2g1[0]  -abs(m1m2g1[4]),2)  /2.0 + pow(m1m2g1[1]  -abs(m1m2g1[3])  -abs(m1m2g1[5])  ,2) /3.0);
   h2m1m2   = sqrt( pow(m1m2g2[0]  -abs(m1m2g2[4]),2)  /2.0 + pow(m1m2g2[1]  -abs(m1m2g2[3])  -abs(m1m2g2[5])  ,2) /3.0);
   h1m2m1   = sqrt( pow(m2m1g1[0]  -abs(m2m1g1[4]),2)  /2.0 + pow(m2m1g1[1]  -abs(m2m1g1[3])  -abs(m2m1g1[5])  ,2) /3.0);
   h2m2m1   = sqrt( pow(m2m1g2[0]  -abs(m2m1g2[4]),2)  /2.0 + pow(m2m1g2[1]  -abs(m2m1g2[3])  -abs(m2m1g2[5])  ,2) /3.0);
   h1m2m1m2 = sqrt( pow(m2m1m2g1[0]-abs(m2m1m2g1[4]),2)/2.0 + pow(m2m1m2g1[1]-abs(m2m1m2g1[3])-abs(m2m1m2g1[5]),2) /3.0);
   h2m2m1m2 = sqrt( pow(m2m1m2g2[0]-abs(m2m1m2g2[4]),2)/2.0 + pow(m2m1m2g2[1]-abs(m2m1m2g2[3])-abs(m2m1m2g2[5]),2) /3.0);
//   go to 1100
   }
//
//C EF boundary
   else if ( ip == 25 ) // 1025 continue
   {
      h1       = sqrt(pow((gvec1[0]   -abs(gvec1[5])),2)   /2.0 + pow((gvec1[1]-abs(gvec1[3])      -abs(gvec1[4])),2)   /3.0);
      h2       = sqrt(pow((gvec2[0]   -abs(gvec2[5])),2)   /2.0 + pow((gvec2[1]-abs(gvec2[3])      -abs(gvec2[4])),2)   /3.0);
      h1m1     = sqrt(pow((m1g1[0]    -abs(m1g1[5])),2)    /2.0 + pow((m1g1[1]-abs(m1g1[3])        -abs(m1g1[4])),2)    /3.0);
      h2m1     = sqrt(pow((m1g2[0]    -abs(m1g2[5])),2)    /2.0 + pow((m1g2[1]-abs(m1g2[3])        -abs(m1g2[4])),2)    /3.0);
      h1m2     = sqrt(pow((m2g1[0]    -abs(m2g1[5])),2)    /2.0 + pow((m2g1[1]-abs(m2g1[3])        -abs(m2g1[4])),2)    /3.0);
      h2m2     = sqrt(pow((m2g2[0]    -abs(m2g2[5])),2)    /2.0 + pow((m2g2[1]-abs(m2g2[3])        -abs(m2g2[4])),2)    /3.0);
      h1m1m2   = sqrt(pow((m1m2g1[0]  -abs(m1m2g1[5])),2)  /2.0 + pow((m1m2g1[1]-abs(m1m2g1[3])    -abs(m1m2g1[4])),2)  /3.0);
      h2m1m2   = sqrt(pow((m1m2g2[0]  -abs(m1m2g2[5])),2)  /2.0 + pow((m1m2g2[1]-abs(m1m2g2[3])    -abs(m1m2g2[4])),2)  /3.0);
      h1m2m1   = sqrt(pow((m2m1g1[0]  -abs(m2m1g1[5])),2)  /2.0 + pow((m2m1g1[1]-abs(m2m1g1[3])    -abs(m2m1g1[4])),2)  /3.0);
      h2m2m1   = sqrt(pow((m2m1g2[0]  -abs(m2m1g2[5])),2)  /2.0 + pow((m2m1g2[1]-abs(m2m1g2[3])    -abs(m2m1g2[4])),2)  /3.0);
      h1m2m1m2 = sqrt(pow((m2m1m2g1[0]-abs(m2m1m2g1[5])),2)/2.0 + pow((m2m1m2g1[1]-abs(m2m1m2g1[3])-abs(m2m1m2g1[4])),2)/3.0);
      h2m2m1m2 = sqrt(pow((m2m1m2g2[0]-abs(m2m1m2g2[5])),2)/2.0 + pow((m2m1m2g2[1]-abs(m2m1m2g2[3])-abs(m2m1m2g2[4])),2)/3.0);
//      go to 1100
   }
//
//
   else if ( ip == 22 ) // 1022 continue
   {
      h1 = sqrt(std::max(0.0,2.0*(gvec1[0]*gvec1[0]+gvec1[1]*gvec1[1]+gvec1[2]*gvec1[2]
           - gvec1[1]*gvec1[2]
           - gvec1[0]*gvec1[2]
           - gvec1[0]*gvec1[1])/3.0));
      h2 = sqrt(std::max(0.0,2.0*(gvec2[0]*gvec2[0]+gvec2[1]*gvec2[1]+gvec2[2]*gvec2[2]
           - gvec2[1]*gvec2[2]
           - gvec2[0]*gvec2[2]
           - gvec2[0]*gvec2[1])/3.0));
      h1m1 = h1;
      h2m1 = h2;
      h1m2 = h1;
      h2m2 = h2;
      h1m1m2 = h1;
      h2m1m2 = h2;
      h1m2m1 = h1;
      h2m2m1 = h2;
      h1m2m1m2 = h1;
      h2m2m1m2 = h2;
   }
   else
   {
      // The fortran had no else branch. This is a place to put an assertion or a breakpoint.
      const int i19191 = 19191;
   }

      //else ?????
//
//
// 1100 continue
//      if (xdebug) then
//        write(7,*) "FOLDMDIST ip ",
//     *  "h1, h2, h1m1, h2m1, h1m2, h2m2 ",
//     *  ip, h1, h2, h1m1, h2m1, h1m2, h2m2
//        write(7,*) "FOLDMDIST ip ",
//     *  "h1m1m2, h2m1m2, h1m2m1, ",
//     *  "h2m2m1, h1m2m1m2, h2m2m1m2 ",
//     *  ip, h1m1m2, h2m1m2, h1m2m1,
//     *  h2m2m1, h1m2m1m2, h2m2m1m2
//      end if
   //
   if (h1 < foldmdist)
   {
      if (h1+h2       < foldmdist) foldmdist = FoldxDist(gvec1,gvec2,   ip,foldmdist,1,1);
      if (h1+h2m1     < foldmdist) foldmdist = FoldxDist(gvec1,m1g2,    ip,foldmdist,1,2);
      if (h1+h2m2     < foldmdist) foldmdist = FoldxDist(gvec1,m2g2,    ip,foldmdist,1,3);
      if (h1+h2m1m2   < foldmdist) foldmdist = FoldxDist(gvec1,m1m2g2,  ip,foldmdist,1,4);
      if (h1+h2m2m1   < foldmdist) foldmdist = FoldxDist(gvec1,m2m1g2,  ip,foldmdist,1,5);
      if (h1+h2m2m1m2 < foldmdist) foldmdist = FoldxDist(gvec1,m2m1m2g2,ip,foldmdist,1,6);
   }

   if (h1m1 < foldmdist)
   {
      if (h1m1+h2       < foldmdist) foldmdist = FoldxDist(m1g1,gvec2,   ip,foldmdist,2,1);
      if (h1m1+h2m1     < foldmdist) foldmdist = FoldxDist(m1g1,m1g2,    ip,foldmdist,2,2);
      if (h1m1+h2m2     < foldmdist) foldmdist = FoldxDist(m1g1,m2g2,    ip,foldmdist,2,3);
      if (h1m1+h2m1m2   < foldmdist) foldmdist = FoldxDist(m1g1,m1m2g2,  ip,foldmdist,2,4);
      if (h1m1+h2m2m1   < foldmdist) foldmdist = FoldxDist(m1g1,m2m1g2,  ip,foldmdist,2,5);
      if (h1m1+h2m2m1m2 < foldmdist) foldmdist = FoldxDist(m1g1,m2m1m2g2,ip,foldmdist,2,6);
   }

   if (h1m2 < foldmdist)
   {
      if (h1m2+h2       < foldmdist) foldmdist = FoldxDist(m2g1,gvec2,   ip,foldmdist,3,1);
      if (h1m2+h2m1     < foldmdist) foldmdist = FoldxDist(m2g1,m1g2,    ip,foldmdist,3,2);
      if (h1m2+h2m2     < foldmdist) foldmdist = FoldxDist(m2g1,m2g2,    ip,foldmdist,3,3);
      if (h1m2+h2m1m2   < foldmdist) foldmdist = FoldxDist(m2g1,m1m2g2,  ip,foldmdist,3,4);
      if (h1m2+h2m2m1   < foldmdist) foldmdist = FoldxDist(m2g1,m2m1g2,  ip,foldmdist,3,5);
      if (h1m2+h2m2m1m2 < foldmdist) foldmdist = FoldxDist(m2g1,m2m1m2g2,ip,foldmdist,3,6);
   }

   if (h1m1m2 < foldmdist)
   {
      if (h1m1m2+h2       < foldmdist) foldmdist = FoldxDist(m1m2g1,gvec2,   ip,foldmdist,4,1);
      if (h1m1m2+h2m1     < foldmdist) foldmdist = FoldxDist(m1m2g1,m1g2,    ip,foldmdist,4,2);
      if (h1m1m2+h2m2     < foldmdist) foldmdist = FoldxDist(m1m2g1,m2g2,    ip,foldmdist,4,3);
      if (h1m1m2+h2m1m2   < foldmdist) foldmdist = FoldxDist(m1m2g1,m1m2g2,  ip,foldmdist,4,4);
      if (h1m1m2+h2m2m1   < foldmdist) foldmdist = FoldxDist(m1m2g1,m2m1g2,  ip,foldmdist,4,5);
      if (h1m1m2+h2m2m1m2 < foldmdist) foldmdist = FoldxDist(m1m2g1,m2m1m2g2,ip,foldmdist,4,6);
   }

   if (h1m2m1 < foldmdist)
   {
      if (h1m2m1+h2       < foldmdist) foldmdist = FoldxDist(m2m1g1,gvec2,   ip,foldmdist,5,1);
      if (h1m2m1+h2m1     < foldmdist) foldmdist = FoldxDist(m2m1g1,m1g2,    ip,foldmdist,5,2);
      if (h1m2m1+h2m2     < foldmdist) foldmdist = FoldxDist(m2m1g1,m2g2,    ip,foldmdist,5,3);
      if (h1m2m1+h2m1m2   < foldmdist) foldmdist = FoldxDist(m2m1g1,m1m2g2,  ip,foldmdist,5,4);
      if (h1m2m1+h2m2m1   < foldmdist) foldmdist = FoldxDist(m2m1g1,m2m1g2,  ip,foldmdist,5,5);
      if (h1m2m1+h2m2m1m2 < foldmdist) foldmdist = FoldxDist(m2m1g1,m2m1m2g2,ip,foldmdist,5,6);
   }

   if (h1m2m1m2 < foldmdist)
   {
      if (h1m2m1m2+h2       < foldmdist) foldmdist = FoldxDist(m2m1m2g1,gvec2,   ip,foldmdist,6,1);
      if (h1m2m1m2+h2m1     < foldmdist) foldmdist = FoldxDist(m2m1m2g1,m1g2,    ip,foldmdist,6,2);
      if (h1m2m1m2+h2m2     < foldmdist) foldmdist = FoldxDist(m2m1m2g1,m2g2,    ip,foldmdist,6,3);
      if (h1m2m1m2+h2m1m2   < foldmdist) foldmdist = FoldxDist(m2m1m2g1,m1m2g2,  ip,foldmdist,6,4);
      if (h1m2m1m2+h2m2m1   < foldmdist) foldmdist = FoldxDist(m2m1m2g1,m2m1g2,  ip,foldmdist,6,5);
      if (h1m2m1m2+h2m2m1m2 < foldmdist) foldmdist = FoldxDist(m2m1m2g1,m2m1m2g2,ip,foldmdist,6,6);
   }

   return( foldmdist );
} //      end FOLDMDIST

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


//C     FOLDXDIST -- return the distance between two
//C     G6 vectors going through the given boundary
//C     (7, A, C or F)
//C
//C     is1 and is2 are the bdy 1,2 symmetry op used
//C

//-----------------------------------------------------------------------------
// Name: FoldxDist()
//       return the distance between two
//       G6 vectors going through the given boundary
//       (7, A, C or F)
//  
//       is1 and is2 are the boundary symmetry ops used
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double MakeGaol::FoldxDist( const arma::vec6& gvec1, const arma::vec6& gvec2, const int ip, const double cFoldDist,
   const int is1DUMMY, const int is2DUMMY) const
{
   is1DUMMY, is2DUMMY; // just to keep the compiler happy
//      REAL*8 FUNCTION FOLDXDIST(GVEC1,GVEC2,IP,CFOLDDIST,
//     *  is1, is2)
//
//      implicit none
//
//      REAL*8 GVEC1(6), GVEC2(6), CFOLDDIST
//      integer ip, kp, lp
//      REAL*8 pg1(6),pg2(6),mpg1(6),mpg2(6)
//      real*8 pfcase1(6),ppfcase1(6),hpfc1,hpfmc1
//      real*8 pfcase2(6),ppfcase2(6),hpfc2,hpfmc2
//      real*8 pfmcase1(6),pfmcase2(6)
//      real*8 ppfmcase1(6),ppfmcase2(6)
//      real*8 fc1(6),fc2(6),fcmp1(6),fcmp2(6)
//      real*8 pfcm1(6),pfcm2(6)
//      real*8 mpfc1(6),mpfc2(6)
//      real*8 mg1(6,3),mg2(6,3)
//      real*8 bd1,bd2
//      real*8 prj(36,25), prjperp(36,25)
//      real*8 prjhat(36,15)
//      integer ms(36,18)
//      common /projectors/ prj,prjperp,prjhat,ms
//      integer i, is1, is2
//      real*8 ht1,ht2,g123dist,g456dist,dist12
//      common/xdebug/xdebug
//      logical xdebug
//      integer ihist(22,36)
//      common /hist/ihist
//      real*8 DistF, anorm
//      logical inncone,lpg1,lmpg1,lpg2,lmpg2
//      logical lpfc1,lpfc2,lfc1,lfc2,lpfmc1,lpfmc2
//
//
//      if (xdebug) then
//        write(7,*) "FoldxDist ip, CFOLDDIST: ",ip,CFOLDDIST
//        write(7,'(a,1x,6f9.2)') "FOLDXDIST gvec1: ",gvec1
//        write(7,'(a,1x,6f9.2)') "FOLDXDIST gvec2: ",gvec2
//      endif
//
   double ht1        (DBL_MAX);
   double ht2        (DBL_MAX);
   double hpfc1      (DBL_MAX);
   double hpfc2      (DBL_MAX);
   double hpfmc1     ( 0.0 );
   double hpfmc2     ( 0.0 );
   bool lpfc1        (DBL_MAX);
   bool lpfc2        (DBL_MAX);
   bool lpfmc1       (DBL_MAX);
   bool lpfmc2       (DBL_MAX);

   double foldxdist = cFoldDist; //      foldxdist = cfolddist
//
//      do i = 1,6
   arma::vec6 pg1( arma::abs(gvec1) ); //        pg1(i) = ABS(GVEC1(i))
   arma::vec6 pg2( arma::abs(gvec2) ); //        pg2(i) = ABS(GVEC2(i))
   arma::vec6 mpg1( pg1 );             //        mpg1(i) = ABS(GVEC1(i))
   arma::vec6 mpg2( pg2 );             //        mpg2(i) = ABS(GVEC2(i))
   arma::vec6 fc1( gvec1 );            //        fc1(i) = GVEC1(i)
   arma::vec6 fc2( gvec2 );            //        fc2(i) = GVEC2(i)

   arma::vec6 pfcase1,ppfcase1;
   arma::vec6 fcmp1, fcmp2;
   arma::vec6 pfcase2,ppfcase2;
   arma::vec6 pfmcase1,pfmcase2;
   arma::vec6 ppfmcase1,ppfmcase2;
//      enddo
//
   bool lfc1( true );
   bool lfc2( true );
   int lp( ip );
   int kp;

//      lfc1 = .true.
//      lfc2 = .true.
//      lp = ip
//
   for ( int i=3; i<6; ++i )             //      do i = 4,6
   {
      if ( fc1[i]<= 0.0 ) lfc1 = ! lfc1; //        if (fc1(i).le.0.D0) lfc1 = .not. lfc1
      if ( fc2[i]<= 0.0 ) lfc2 = ! lfc2; //        if (fc2(i).lt.0.D0) lfc2 = .not. lfc2
   }                                     //      enddo
   //
   //C     Process 8F boundary 
   if ( ip==23 )                         //      if (ip.eq.23) then
   {                                     
      if ( lfc1 )                        //        if (lfc1) then
      {                                  
         fc1[3] = -fc1[3];               //          fc1(4) = -fc1(4)
         if ( abs(fc1[4])>=abs(fc1[5]) ) //          if (dabs(fc1(5)).ge.dabs(fc1(6))) then
            fc1[4] = -fc1[4];            //            fc1(5) = -fc1(5)
         else                            //          else
            fc1[5] = -fc1[5];            //            fc1(6) = -fc1(6)
      }                                  //        endif
      if ( lfc2 )                        //        if (lfc2) then
      {
         fc2[3] = -fc2[3];               //          fc2(4) = -fc2(4)
         if ( abs(fc2[4])>=abs(fc2[5]) ) //          if (dabs(fc2(5)).ge.dabs(fc2(6))) then
            fc2[4] = -fc2[4];            //            fc2(5) = -fc2(5)
         else                            //          else
            fc2[5] = -fc2[5];            //            fc2(6) = -fc2(6)
      }                                  //        endif

      pg1[0]  = (2.0*fc1[0]-fc1[4]-fc1[5])/3.0;  //        pg1(1) = (2.0*fc1(1)-fc1(5)-fc1(6))/3.0
      pg2[0]  = (2.0*fc2[0]-fc2[4]-fc2[5])/3.0;  //        pg2(1) = (2.0*fc2(1)-fc2(5)-fc2(6))/3.0
      mpg1[0] = pg1[0];                         //        mpg1(1) = pg1(1)
      mpg2[0] = pg2[0];                         //        mpg2(1) = pg2(1)
      pg1[1]  = (fc1[1]-fc1[3])/2.0;             //        pg1(2) = (fc1(2)-fc1(4))/2.0
      pg2[1]  = (fc2[1]-fc2[3])/2.0;             //        pg2(2) = (fc2(2)-fc2(4))/2.0
      mpg1[1] = pg1[1];                         //        mpg1(2) = pg1(2)
      mpg2[1] = pg2[1];                         //        mpg2(2) = pg2(2)
      pg1[3]  = -pg1[1];                         //        pg1(4) = -pg1(2)
      pg2[3]  = -pg2[1];                         //        pg2(4) = -pg2(2)
      mpg1[3] = pg1[1];                         //        mpg1(4) = pg1(2)
      mpg2[3] = pg2[1];                         //        mpg2(4) = pg2(2)
      pg1[4]  = (2.0*fc1[4]-fc1[0]-    fc1[5])/3.0;  //        pg1(5) = (2.0*fc1(5)-fc1(1)-fc1(6))/3.0
      pg2[4]  = (2.0*fc2[4]-fc2[0]-    fc2[5])/3.0;  //        pg2(5) = (2.0*fc2(5)-fc2(1)-fc2(6))/3.0
      mpg1[4] = (2.0*fc1[0]-fc1[4]-    fc1[5])/3.0; //        mpg1(5) = (2.0*fc1(1)-fc1(5)-fc1(6))/3.0
      mpg2[4] = (2.0*fc2[0]-fc2[4]-    fc2[5])/3.0; //        mpg2(5) = (2.0*fc2(1)-fc2(5)-fc2(6))/3.0
      pg1[5]  = (2.0*fc1[5]-fc1[0]-    fc1[4])/3.0;  //        pg1(6) = (2.0*fc1(6)-fc1(1)-fc1(5))/3.0
      pg2[5]  = (2.0*fc2[5]-fc2[0]-    fc2[4])/3.0;  //        pg2(6) = (2.0*fc2(6)-fc2(1)-fc2(5))/3.0
      mpg1[5] =     (fc1[0]+fc1[4]-2.0*fc1[5])/3.0; //        mpg1(6) = (fc1(1)+fc1(5)-2.0*fc1(6))/3.0
      mpg2[5] =     (fc2[0]+fc2[4]-2.0*fc2[5])/3.0; //        mpg2(6) = (fc2(1)+fc2(5)-2.0*fc2(6))/3.0

      ht1     = sqrt(pow(fc1[1]+fc1[3],2)/2.0 + pow(fc1[0]+fc1[4]+fc1[5],2)/3.0);
      ht2     = sqrt(pow(fc2[1]+fc2[3],2)/2.0 + pow(fc2[0]+fc2[4]+fc2[5],2)/3.0);
      hpfc1   = 9.E38;
      hpfc2   = 9.E38;
      hpfmc1  = 9.E38;
      hpfmc2  = 9.E38;
      lpfc1   = false;
      lpfc2   = false;
      lpfmc1  = false;
      lpfmc2  = false;
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
      //      go to 1000
   }//      endif
   //
   //
   //C     Process BF boundary 
   else if ( ip==24 )                            //      if (ip.eq.24) then
   {
      if ( lfc1 )                                //        if (lfc1) then
      {
         fc1[4] = -fc1[4];                       //          fc1(5) = -fc1(5)
         if ( abs(fc1[3])>=abs(fc1[5]) )         //          if (dabs(fc1(4)).ge.dabs(fc1(6))) then
            fc1[3] = -fc1[3];                    //            fc1(4) = -fc1(4)
         else                                    //          else
            fc1[5] = -fc1[5];                    //            fc1(6) = -fc1(6)
      }                                          //        endif
      if ( lfc2 ) //        if (lfc2) then
      {
         fc2[4] = -fc2[4];                       //          fc2(5) = -fc2(5)
         if ( abs(fc2[3])>=abs(fc2[5]) )         //          if (dabs(fc2(4)).ge.dabs(fc2(6))) then
            fc2[3] = -fc2[3];                    //            fc2(4) = -fc2(4)
         else                                    //          else
            fc2[5] = -fc2[5];                    //            fc2(6) = -fc2(6)
      }                                          //        endif
        pg1[1]  = (2.0*fc1[1]-fc1[3]-fc1[5])/3.0;
        pg2[1]  = (2.0*fc2[1]-fc2[3]-fc2[5])/3.0;
        mpg1[1] = pg1[1];
        mpg2[1] = pg2[1];
        pg1[0]  = (fc1[0]-fc1[4])/2.0;
        pg2[0]  = (fc2[0]-fc2[4])/2.0;
        mpg1[0] = pg1[0];
        mpg2[0] = pg2[0];
        pg1[4]  = -pg1[0];
        pg2[4]  = -pg2[0];
        mpg1[4] = pg1[0];
        mpg2[4] = pg2[0];
        pg1[3]  =  (2.0*fc1[3]-fc1[1]-    fc1[5])/3.0;
        pg2[3]  =  (2.0*fc2[3]-fc2[1]-    fc2[5])/3.0;
        mpg1[3] =  (2.0*fc1[1]-fc1[3]-    fc1[5])/3.0;
        mpg2[3] =  (2.0*fc2[1]-fc2[3]-    fc2[5])/3.0;
        pg1[5]  =  (2.0*fc1[5]-fc1[1]-    fc1[3])/3.0;
        pg2[5]  =  (2.0*fc2[5]-fc2[1]-    fc2[3])/3.0;
        mpg1[5] =      (fc1[1]+fc1[3]-2.0*fc1[5])/3.0;
        mpg2[5] =      (fc2[1]+fc2[3]-2.0*fc2[5])/3.0;
       
        ht1     = sqrt( pow(fc1[0]+fc1[4],2)/2.0 + pow(fc1[1]+fc1[3]+fc1[5],2)/3.0);
        ht2     = sqrt( pow(fc2[0]+fc2[4],2)/2.0 + pow(fc2[1]+fc2[3]+fc2[5],2)/3.0);
        hpfc1   = 9.E38;
        hpfc2   = 9.E38;
        hpfmc1  = 9.E38;
        hpfmc2  = 9.E38;
        lpfc1   = false;
        lpfc2   = false;
        lpfmc1  = false;
        lpfmc2  = false;
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
//      go to 1000
   }//      endif
//
//C     Process EF boundary 
   else if ( ip==25 )                    //      if (ip.eq.25) then
   {
      if ( lfc1 )                        //        if (lfc1) then
      {
         fc1[5] = -fc1[5];               //          fc1(6) = -fc1(6)
         if ( abs(fc1[4])>=abs(fc2[3]) ) //          if (dabs(fc1(5)).ge.dabs(fc1(4))) then
            fc1[4] = -fc1[4];            //            fc1(5) = -fc1(5)
         else                            //          else
            fc1[3] = -fc1[3];            //            fc1(4) = -fc1(4)
      }                                  //        endif
      if( lfc2 )                         //        if (lfc2) then
      {
         fc2[5] = -fc2[5];               //          fc2(6) = -fc2(6)
         if ( abs(fc2[4])>=abs(fc2[3]) ) //          if (dabs(fc2(5)).ge.dabs(fc2(4))) then
            fc2[4] = -fc2[4];            //            fc2(5) = -fc2(5)
         else                            //          else
            fc2[3] = -fc2[3];            //            fc2(4) = -fc2(4)
      }                                  //        endif
      pg1[1]  = (2.0*fc1[1]-fc1[3]-fc1[4])/3.0;
      pg2[1]  = (2.0*fc2[1]-fc2[3]-fc2[4])/3.0;
      mpg1[1] = pg1[1];
      mpg2[1] = pg2[1];
      pg1[0]  = (fc1[0]-fc1[5])/2.0;
      pg2[0]  = (fc2[0]-fc2[5])/2.0;
      mpg1[0] = pg1[0];
      mpg2[0] = pg2[0];
      pg1[5]  = -pg1[0];
      pg2[5]  = -pg2[0];
      mpg1[5] = pg1[0];
      mpg2[5] = pg2[0];
      pg1[3]  =  (2.0*fc1[3]-fc1[1]-    fc1[4])/3.0;
      pg2[3]  =  (2.0*fc2[3]-fc2[1]-    fc2[4])/3.0;
      mpg1[3] =  (2.0*fc1[1]-fc1[3]-    fc1[4])/3.0;
      mpg2[3] =  (2.0*fc2[1]-fc2[3]-    fc2[4])/3.0;
      pg1[4]  =  (2.0*fc1[4]-fc1[1]-    fc1[3])/3.0;
      pg2[4]  =  (2.0*fc2[4]-fc2[1]-    fc2[3])/3.0;
      mpg1[4] =      (fc1[1]+fc1[3]-2.0*fc1[4])/3.0;
      mpg2[4] =      (fc2[1]+fc2[3]-2.0*fc2[4])/3.0;

      ht1 = sqrt(pow(fc1[0]+fc1[5],2)/2.0 + pow(fc1[1]+fc1[3]+fc1[4],2)/3.0);
      ht2 = sqrt(pow(fc2[0]+fc2[5],2)/2.0 + pow(fc2[1]+fc2[3]+fc2[4],2)/3.0);
      hpfc1  = 9.E38;
      hpfc2  = 9.E38;
      hpfmc1 = 9.E38;
      hpfmc2 = 9.E38;
      lpfc1  = false;
      lpfc2  = false;
      lpfmc1 = false;
      lpfmc2 = false;
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
//      go to 1000
   }//      endif
//
//
   else if ( ip==6 || ip==7 || ip==8 ) //      if (ip.eq.7.or.ip.eq.6.or.ip.eq.8) then
   {
        pg1[5]  = gvec1[5];
        pg2[5]  = gvec2[5];
        pg1[1]  = (pg1[1]+pg1[3])/2.0;
        pg1[3]  = pg1[1];
        pg2[1]  = (pg2[1]+pg2[3])/2.0;
        pg2[3]  = pg2[1];
        mpg1[1] = pg1[1];
        mpg1[3] = pg1[1];
        mpg2[1] = pg2[1];
        mpg2[3] = pg2[1];
        mpg1[5] = pg1[5];
        mpg2[5] = pg2[5];
        mpg1[4] = pg1[5]-pg1[4];
        mpg2[4] = pg2[5]-pg2[4];
        ht1 = abs(gvec1[1]-abs(gvec1[3]))/sqrt(2.0);
        ht2 = abs(gvec2[1]-abs(gvec2[3]))/sqrt(2.0);
        if ( lfc1 )                        //        if (lfc1) then
        {
           fc1[3] = -fc1[3];               //          fc1(4) = -fc1(4)
           if ( abs(fc1[4])>=abs(fc1[5]) ) //          if (dabs(fc1(5)).ge.dabs(fc1(6))) then
              fc1[4] = -fc1[4];            //            fc1(5) = -fc1(5)
           else                            //          else
              fc1[5] = -fc1[5];            //            fc1(6) = -fc1(6)
        }                                  //        endif
        if ( lfc2 )                        //        if (lfc2) then
        {
           fc2[3] = -fc2[3];               //          fc2(4) = -fc2(4)
           if ( abs(fc2[4])>=abs(fc2[5]) ) //          if (dabs(fc2(5)).ge.dabs(fc2(6))) then
              fc2[4] = -fc2[4];            //            fc2(5) = -fc2(5)
           else                            //          else
              fc2[5] = -fc2[5];            //            fc2(6) = -fc2(6)
        }                                  //        endif
        ppfcase1 = m_prjList[23].GetPerp() * fc1; //        call rmv6(fc1,prjperp(1,23),ppfcase1)
        hpfc1 = arma::norm(ppfcase1,2);//        hpfc1 = anorm(6,ppfcase1)
        ppfcase2 = m_prjList[23].GetPerp() * fc2; //        call rmv6(fc2,prjperp(1,23),ppfcase2)
        hpfc2 = arma::norm(ppfcase2,2); //        hpfc2 = anorm(6,ppfcase2)
        pfcase1 = m_prjList[23].GetProjector() * fc1; //        call rmv6(fc1,prj(1,23),pfcase1)
        pfcase2 = m_prjList[23].GetProjector() * fc2; //        call rmv6(fc2,prj(1,23),pfcase2)
        lpfc1=inncone(pfcase1); //        lpfc1=inncone(pfcase1)
        lpfc2=inncone(pfcase2); //        lpfc2=inncone(pfcase2)
        kp = 8; //        kp = 8
        lpfmc1 = true; //        lpfmc1 = .true.
        lpfmc2 = true; //        lpfmc2 = .true.
        for ( int i=0; i<6; ++i ) //        do i = 1,6
        {
           fcmp1[i] = mpg1[i]; //          fcmp1(i) = mpg1(i)
           if ( fcmp1[i] <= 0.0 )  lpfmc1 = ! lpfmc1; //          if (fcmp1(i).le.0.D0) lpfmc1=.not.lpfmc1
           fcmp2[i] = mpg2[i]; //          fcmp2(i) = mpg2(i)
           if ( fcmp2[i] <= 0.0 ) lpfmc2 = ! lpfmc2; //          if (fcmp2(i).le.0.D0) lpfmc2=.not.lpfmc2
        }//        enddo
        for ( int i=3; i<6; ++i ) //        do i = 4,6
        {
           fcmp1[i] = -abs(fcmp1[i]); //          fcmp1(i) = -dabs(fcmp1(i))
           fcmp2[i] = -abs(fcmp2[i]); //          fcmp2(i) = -dabs(fcmp2(i))
        }//        enddo
        if ( lpfmc1 ) //        if (lpfmc1) then
        {
           if ( fcmp1[4] < fcmp1[5] ) //          if (fcmp1(5) .lt. fcmp1(6)) then
              fcmp1[5] = -fcmp1[5]; //            fcmp1(6) = -fcmp1(6)
           else //          else
              fcmp1[4] = -fcmp1[4]; //            fcmp1(5) = -fcmp1(5)
        }//        endif
        if ( lpfmc2 ) //        if (lpfmc2) then
        {
           if ( fcmp2[4] < fcmp2[5] ) //          if (fcmp2(5) .lt. fcmp2(6)) then
              fcmp2[5] = -fcmp2[5]; //            fcmp2(6) = -fcmp2(6)
           else //          else
              fcmp2[4] = -fcmp2[4]; //            fcmp2(5) = -fcmp2(5)
        }//        endif
        lp = 23; //        lp = 23
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
   }//        go to 1000
   //
   else if ( ip==10 || ip==9 || ip==11 ) //      else if (ip.eq.10.or.ip.eq.9.or.ip.eq.11) then
   {
      pg1[5]  = gvec1[5];
      pg2[5]  = gvec2[5];
      pg1[0]  =(pg1[0]+pg1[4])/2.0;
      pg1[4]  = pg1[0];
      pg2[0]  =(pg2[0]+pg2[4])/2.0;
      pg2[4]  = pg2[0];
      mpg1[0] = pg1[0];
      mpg1[4] = pg1[0];
      mpg2[0] = pg2[0];
      mpg2[4] = pg2[0];
      mpg1[5] = pg1[5];
      mpg2[5] = pg2[5];
      mpg1[3] = pg1[5]-pg1[3];
      mpg2[3] = pg2[5]-pg2[3];
      ht1 = abs(gvec1[0]-abs(gvec1[4]))/sqrt(2.0);
      ht2 = abs(gvec2[0]-abs(gvec2[4]))/sqrt(2.0);
      if ( lfc1 )                        //        if (lfc1) then
      {
         fc1[4] = -fc1[4];               //          fc1(5) = -fc1(5)
         if ( abs(fc1[3])>=abs(fc1[5]) ) //          if (dabs(fc1(4)).ge.dabs(fc1(6))) then
            fc1[3] = -fc1[3];            //            fc1(4) = -fc1(4)
         else                            //          else
            fc1[5] = -fc1[5];            //            fc1(6) = -fc1(6)
      }                                  //        endif
      if ( lfc2 )                        //        if (lfc2) then
      {
         fc2[4] = -fc2[4];               //          fc2(5) = -fc2(5)
         if ( abs(fc2[3])>=abs(fc2[5]) ) //          if (dabs(fc2(4)).ge.dabs(fc2(6))) then
            fc2[3] = -fc2[3];            //            fc2(4) = -fc2(4)
         else                            //          else
            fc2[5] = -fc2[5];            //            fc2(6) = -fc2(6)
      }                                  //        endif
      ppfcase1 = this->m_prjList[24].GetPerp() * fc1; //        call rmv6(fc1,prjperp(1,24),ppfcase1)
      hpfc1 = arma::norm( ppfcase1,2 );               //        hpfc1 = anorm(6,ppfcase1)
      ppfcase2 = this->m_prjList[24].GetPerp() * fc2; //        call rmv6(fc2,prjperp(1,24),ppfcase2)
      hpfc2 = arma::norm( ppfcase2,2 );               //        hpfc2 = anorm(6,ppfcase2)
      pfcase1 = m_prjList[24].GetProjector() * fc1;   //        call rmv6(fc1,prj(1,24),pfcase1)
      pfcase2 = m_prjList[24].GetProjector() * fc2;   //        call rmv6(fc2,prj(1,24),pfcase2)
      lpfc1 = inncone(pfcase1);                       //        lpfc1=inncone(pfcase1)
      lpfc2 = inncone(pfcase2);                       //        lpfc2=inncone(pfcase2)
      kp = 11;                                        //        kp = 11
      lp = 24;                                        //        lp = 24
      lpfmc1 = true;                                  //        lpfmc1 = .true.
      lpfmc2 = true;                                  //        lpfmc2 = .true.
      for( int i=0; i<6; ++i )                     //        do i = 1,6
      {
         fcmp1[i] = mpg1[i];                       //          fcmp1(i) = mpg1(i)
         if ( fcmp1[i] <= 0.0 ) lpfmc1 = ! lpfmc1; //          if (fcmp1(i).le.0.D0) lpfmc1=.not.lpfmc1
         fcmp2[i] = mpg2[i];                       //          fcmp2(i) = mpg2(i)
         if ( fcmp2[i] <= 0.0 ) lpfmc2 = ! lpfmc2; //          if (fcmp2(i).le.0.D0) lpfmc2=.not.lpfmc2
      }                                            //        enddo
      for ( int i=3; i<6; ++i )                    //        do i = 4,6
      {
         fcmp1[i] = -abs(fcmp1[i]);                //          fcmp2(i) = -dabs(fcmp2(i))
         fcmp2[i] = -abs(fcmp2[i]);
      }                                            //        enddo
      if ( lpfmc1 )                                //        if (lpfmc1) then
      {
         if ( fcmp1[3] < fcmp1[5] )                //          if (fcmp1(4) .lt. fcmp1(6)) then
            fcmp1[5] = -fcmp1[5];                  //            fcmp1(6) = -fcmp1(6)
         else                                      //          else
            fcmp1[3] = -fcmp1[3];                  //            fcmp1(4) = -fcmp1(4)
      }                                            //        endif
      if ( lpfmc2 )                 //        if (lpfmc2) then
      {
         if ( fcmp2[3] < fcmp2[5] ) //          if (fcmp2(4) .lt. fcmp2(6)) then
            fcmp2[5] = -fcmp2[5];   //            fcmp2(6) = -fcmp2(6)
         else                       //          else
            fcmp2[3] = -fcmp2[3];   //            fcmp2(4) = -fcmp2(4)
      }                             //        endif
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
      //        go to 1000
}//
   else if ( ip==13 || ip==12 || ip==14 ) //      else if (ip.eq.13.or.ip.eq.12.or.ip.eq.14) then
   {
      pg1[4]  = gvec1[4];
      pg2[4]  = gvec2[4];
      pg1[0]  =(pg1[0]+pg1[5])/2.0;
      pg1[5]  = pg1[0];
      pg2[0]  =(pg2[0]+pg2[5])/2.0;
      pg2[5]  = pg2[0];
      mpg1[0] = pg1[0];
      mpg1[5] = pg1[0];
      mpg2[0] = pg2[0];
      mpg2[5] = pg2[0];
      mpg1[4] = pg1[4];
      mpg2[4] = pg2[4];
      mpg1[3] = pg1[4]-pg1[3];
      mpg2[3] = pg2[4]-pg2[3];
      ht1 = abs(gvec1[0]-abs(gvec1[5]))/sqrt(2.0);
      ht2 = abs(gvec2[0]-abs(gvec2[5]))/sqrt(2.0);
      if ( lfc1 )                          //        if (lfc1) then
      {
         fc1[5] = -fc1[5];                 //          fc1(6) = -fc1(6)
         if ( abs(fc1[4]) >= abs(fc1[3]) ) //          if (dabs(fc1(5)).ge.dabs(fc1(4))) then
            fc1[4] = -fc1[4];              //            fc1(5) = -fc1(5)
         else                              //          else
            fc1[3] = -fc1[3];              //            fc1(4) = -fc1(4)
      }                                    //        endif
      if ( lfc2 )                          //        if (lfc2) then
      {
         fc2[5] = -fc2[5];                 //          fc2(6) = -fc2(6)
         if ( abs(fc2[4]) >= abs(fc2[3]) ) //          if (dabs(fc2(5)).ge.dabs(fc2(4))) then
            fc2[4] = -fc2[4];              //            fc2(5) = -fc2(5)
         else                              //          else
            fc2[3] = -fc2[3];              //            fc2(4) = -fc2(4)
      }                                    //        endif
      ppfcase1 = m_prjList[25].GetPerp() * fc1;     //        call rmv6(fc1,prjperp(1,25),ppfcase1)
      hpfc1 = arma::norm( ppfcase1,2 );             //        hpfc1 = anorm(6,ppfcase1)
      ppfcase2 = m_prjList[25].GetPerp() * fc2;     //        call rmv6(fc2,prjperp(1,25),ppfcase2)
      hpfc2 = arma::norm( ppfcase2,2 );             //        hpfc2 = anorm(6,ppfcase2)
      pfcase1 = m_prjList[25].GetProjector() * fc1; //        call rmv6(fc1,prj(1,25),pfcase1)
      pfcase2 = m_prjList[25].GetProjector() * fc2; //        call rmv6(fc2,prj(1,25),pfcase2)
      lpfc1 = inncone( pfcase1 );                   //        lpfc1=inncone(pfcase1)
      lpfc2 = inncone( pfcase2 );                   //        lpfc2=inncone(pfcase2)
      kp = 14;                                      //        kp = 14
      lp = 25;                                      //        lp = 25
      lpfmc1 = true;                                //        lpfmc1 = .true.
      lpfmc2 = true;                                //        lpfmc2 = .true.
      for ( int i=0; i<6; ++i )                  //        do i = 1,6
      {
         fcmp1[i] = mpg1[i];                     //          fcmp1(i) = mpg1(i)
         if ( fcmp1[i]<=0.0 ) lpfmc1 = ! lpfmc1; //          if (fcmp1(i).le.0.D0) lpfmc1=.not.lpfmc1
         fcmp2[i] = mpg2[i];                     //          fcmp2(i) = mpg2(i)
         if ( fcmp2[i]<=0.0 ) lpfmc2 = ! lpfmc2; //          if (fcmp2(i).le.0.D0) lpfmc2=.not.lpfmc2
      }                                          //        enddo
      for ( int i=3; i<6; ++i )     //        do i = 4,6
      {
         fcmp1[i] = -abs(fcmp1[i]); //          fcmp1(i) = -dabs(fcmp1(i))
         fcmp2[i] = -abs(fcmp2[i]); //          fcmp2(i) = -dabs(fcmp2(i))
      }                             //        enddo
      if ( lpfmc1 )                 //        if (lpfmc1) then
      {
         if ( fcmp1[3] < fcmp1[4] ) //          if (fcmp1(4) .lt. fcmp1(5)) then
            fcmp1[4] = -fcmp1[4];   //            fcmp1(5) = -fcmp1(5)
         else//          else
            fcmp1[3] = -fcmp1[3];   //            fcmp1(4) = -fcmp1(4)
      }                             //        endif
      if ( lpfmc2 ) //        if (lpfmc2) then
      {
         if ( fcmp2[3] < fcmp2[4] ) //          if (fcmp2(4) .lt. fcmp2(5)) then
            fcmp2[4] = -fcmp2[4]; //            fcmp2(5) = -fcmp2(5)
         else //          else
            fcmp2[3] = -fcmp2[3]; //            fcmp2(4) = -fcmp2(4)
      }//        endif
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
      //        go to 1000
   }         //
   else if (ip==15) 
   {
      const double bd1 = gvec1[3]+gvec1[4]+gvec1[5];
      const double bd2 = gvec2[3]+gvec2[4]+gvec2[5];
      pg1[0]  = .60*gvec1[0]-0.20*bd1;
      pg2[0]  = .60*gvec2[0]-0.20*bd2;
      pg1[1]  = pg1[0];
      pg2[1]  = pg2[0];
      pg1[3]  = -0.40*gvec1[0]+gvec1[3]-0.20*bd1;
      pg2[3]  = -0.40*gvec2[0]+gvec2[3]-0.20*bd2;
      pg1[4]  = -0.40*gvec1[0]+gvec1[4]-0.20*bd1;
      pg2[4]  = -0.40*gvec2[0]+gvec2[4]-0.20*bd2;
      pg1[5]  = -0.40*gvec1[0]+gvec1[5]-0.20*bd1;
      pg2[5]  = -0.40*gvec2[0]+gvec2[5]-0.20*bd2;
      mpg1[0] = pg1[0];
      mpg2[0] = pg2[0];
      mpg1[1] = pg1[1];
      mpg2[1] = pg2[1];
      mpg1[3] = pg1[4];
      mpg2[3] = pg2[4];
      mpg1[4] = pg1[3];
      mpg2[4] = pg2[3];
      mpg1[5] = pg1[5];
      mpg2[5] = pg2[5];
      ht1     = DistF(gvec1);
      ht2     = DistF(gvec2);
      hpfc1   = 9.E38;
      hpfc2   = 9.E38;
      hpfmc1  = 9.E38;
      hpfmc2  = 9.E38;
      lpfc1   = false;
      lpfc2   = false;
      lpfmc1  = false;
      lpfmc2  = false;
   } //endif
//      go to 2000
//
// 1000 continue
//void  Set1000( pfmcase1, fcmp1, lp, hpfmc1, ht1, [output] hpfmc1, lpfmc1 )
//      call rmv6(fcmp1,prj(1,lp),pfmcase1)
//      call rmv6(fcmp2,prj(1,lp),pfmcase2)
//      hpfmc1 = g456dist(pfmcase1,fcmp1)
//      hpfmc1 = sqrt(hpfmc1**2+ht1**2)
//      hpfmc2 = g456dist(pfmcase2,fcmp2)
//      hpfmc1 = sqrt(hpfmc2**2+ht2**2)
//      lpfmc1=inncone(pfmcase1)
//      lpfmc2=inncone(pfmcase2)
//      if (xdebug) then
//        write(7,'(a,f12.4,f12.4,1x,6f9.2,1x,6f9.2)')
//     *   'hpfmc1,ht1,fcmp1,pfmcase1',
//     *   hpfmc1,ht1,fcmp1,pfmcase1
//        write(7,'(a,f12.4,f12.4,1x,6f9.2,1x,6f9.2)')
//     *   'hpfmc2,ht2,fcmp2,pfmcase2',
//     *   hpfmc2,ht2,fcmp2,pfmcase2
//      endif
//
// 2000 continue

   const bool lpg1  = inncone(pg1);
   const bool lpg2  = inncone(pg2);
   const bool lmpg1 = inncone(mpg1);
   const bool lmpg2 = inncone(mpg2);
   hpfc1  = 9.E38;
   hpfc2  = 9.E38;
   hpfmc1 = 9.E38;
   hpfmc2 = 9.E38;
   lpfc1  = false;
   lpfc2  = false;
   lpfmc1 = false;
   lpfmc2 = false;

   //
   double dist12;
   if ( ht1+ht2 < foldxdist ) //      if (ht1+ht2 .lt. FoldxDist) then
   {
      dist12 = foldxdist; //        dist12 = FoldxDist
      if ( lpg1 ) //        if (lpg1) then
      {
         if ( lpg2  ) dist12 = std::min(dist12,BasicDistance::g123dist(pg1, pg2)); //          if (lpg2)
         //     *      dist12=min(dist12,g123dist(pg1,pg2))
         if ( lmpg2 ) dist12 = std::min(dist12,BasicDistance::g123dist(pg1,mpg2)); //          if (lmpg2)
         //     *      dist12=min(dist12,g123dist(pg1,mpg2))
      } //        endif
      if ( lmpg1 )//        if (lmpg1) then
      {
         if ( lpg2  ) dist12 = std::min(dist12,BasicDistance::g123dist(mpg1, pg2)); //          if (lpg2)
         //     *      dist12=min(dist12,g123dist(mpg1,pg2))
         if ( lmpg2 ) dist12 = std::min(dist12,BasicDistance::g123dist(mpg1,mpg2));//          if (lmpg2)
         //     *      dist12=min(dist12,g123dist(mpg1,mpg2))
      } //        endif
      foldxdist = std::min(foldxdist,sqrt( pow(ht1+ht2,2) + dist12*dist12) ); //        FoldxDist = MIN(FoldxDist,
      //     *    DSQRT((ht1+ht2)**2 +
      //     *    dist12**2))
      //        if (xdebug)
      //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
   } //      endif
   if ( hpfc1+ht2 < foldxdist ) //      if(hpfc1+ht2 .lt.FOLDXDIST) then
   {
      if ( lpfc1 && (lpg2 || lmpg2) )//        if (lpfc1
      {
         //     *    .and.(lpg2.or.lmpg2)) then
         dist12 = foldxdist; //        dist12 = FoldxDist
         if ( lpg2  ) dist12 = std::min(dist12,BasicDistance::g123dist(pfcase1,pg2)); //        if (lpg2)
         //     *  dist12 = min(dist12,
         //     *    g123dist(pfcase1,pg2))
         if ( lmpg2 ) dist12 = std::min(dist12,BasicDistance::g123dist(pfcase1,mpg2)); //        if (lmpg2)
         //     *  dist12 = min(dist12,
         //     *    g123dist(pfcase1,mpg2))
         if( lpg2 || lmpg2 ) foldxdist = std::min( foldxdist, sqrt( pow(hpfc1+ht2,2) + dist12*dist12 ) ); //        if (lpg2.or.lmpg2)
         //     *    FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((hpfc1+ht2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,3f12.4,6f9.2)')
         //     *  "FoldxDist: hpfc1, ht2, dist12",
         //     *   hpfc1, ht2, dist12, FoldxDist
      }  //        endif
   }     //      endif
   if ( ht1+hpfc2 < foldxdist )//      if(ht1+hpfc2 .lt.FOLDXDIST) then
   {
      if ( lpfc2 && (lpg1 || lmpg1 ) ) //        if (lpfc2
      {
         //     *    .and.(lpg1.or.lmpg1)) then
         dist12 = foldxdist; //        dist12 = FoldxDist
         if ( lpg1  ) dist12 = std::min( dist12,BasicDistance::g123dist(pg1,pfcase2)); //        if (lpg1)
         //     *    dist12 = min(dist12,
         //     *      g123dist(pg1,pfcase2))
         if ( lmpg1 ) dist12 = std::min(dist12,BasicDistance::g123dist(mpg1,pfcase2));//        if (lmpg1)
         //     *    dist12 = min(dist12,
         //     *      g123dist(mpg1,pfcase2))
         if ( lpg1 || lmpg1 ) foldxdist = std::min( foldxdist, sqrt( pow(ht1+hpfc2,2) + dist12*dist12) ); //        if (lpg1.or.lmpg1)
         //     *    FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((ht1+hpfc2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,3f12.4,6f9.2)')
         //     *  "FoldxDist: ht1, hpfc2, dist12",
         //     *   ht1, hpfc2, dist12, FoldxDist
      } //        endif
   } //      endif

   arma::vec6 mpfc1, mpfc2;
   if ( hpfc1+hpfc2 < foldxdist ) //      if(hpfc1+hpfc2 .lt.FOLDXDIST) then
   {
      if ( lpfc1 && lpfc2 ) //        if (lpfc1.and.lpfc2)  then
      {
         dist12 = BasicDistance::g123dist(pfcase1,pfcase2); //        dist12 = g123dist(pfcase1,pfcase2)
         foldxdist = std::min( foldxdist, sqrt( pow(hpfc1+hpfc2,2) + dist12*dist12 ) );//        FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((hpfc1+hpfc2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      } //        endif
   } //      endif
   if( lpfc1 && (hpfc1+std::min(ht2,hpfc2)<foldxdist ) )//      if (lpfc1 .and.
   {
      //     *  (hpfc1+min(ht2,hpfc2).lt.FoldxDist)) then
      mpfc1 = ms[kp] * pfcase1;//        call imv6(pfcase1,MS(1,kp),mpfc1)
   } //      endif
   if ( lpfc2 && std::min(ht1,hpfc1)+hpfc2 < foldxdist )//      if (lpfc2 .and.
   {
      //     *  (min(ht1,hpfc1)+hpfc2.lt.FoldxDist)) then
      mpfc2 = ms[kp] * pfcase2; //        call imv6(pfcase2,MS(1,kp),mpfc2)
   } //      endif
   if ( hpfc1+ht2 < foldxdist ) //      if(hpfc1+ht2 .lt.FOLDXDIST) then
   {
      if ( lpfc1 && (lpg2 || lmpg2 ) ) //        if (lpfc1
      {
         //     *    .and.(lpg2.or.lmpg2)) then
         dist12 = foldxdist; //        dist12 = FoldxDist
         if ( lpg2  ) dist12 = std::min( dist12,BasicDistance::g123dist(mpfc1,pg2)); //        if (lpg2)
         //     *  dist12 = min(dist12,
         //     *    g123dist(mpfc1,pg2))
         if ( lmpg2 ) dist12 = std::min( dist12,BasicDistance::g123dist(mpfc1,mpg2) ); //        if (lmpg2)
         //     *  dist12 = min(dist12,
         //     *    g123dist(mpfc1,mpg2))
         if ( lpg2 || lmpg2 ) foldxdist = std::min(foldxdist, sqrt( pow(hpfc1+ht2,2) + dist12*dist12 ) ); //        if (lpg2.or.lmpg2)
         //     *    FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((hpfc1+ht2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }//        endif
   }//      endif
   if ( ht1+hpfc2 < foldxdist ) //      if(ht1+hpfc2 .lt.FOLDXDIST) then
   {
      if ( lpfc2 && (lpg1 || lmpg1 ) ) //        if (lpfc2
      {
         //     *    .and.(lpg1.or.lmpg1)) then
         //        if (xdebug) then
         //          write(7,'(a,1x,$)')'pg1,mpfc2'
         //        endif
         dist12 = foldxdist; //        dist12 = FoldxDist
         if ( lpg1  ) dist12 = std::min( dist12,BasicDistance::g123dist(pg1,mpfc2)); //        if (lpg1)
         //     *    dist12 = min(dist12,
         //     *      g123dist(pg1,mpfc2))
         if ( lmpg1 ) dist12 = std::min( dist12,BasicDistance::g123dist(mpg1,mpfc2)); //        if (lmpg1)
         //     *    dist12 = min(dist12,
         //     *      g123dist(mpg1,mpfc2))
         if ( lpg1 || lmpg1 ) foldxdist = std::min(foldxdist,sqrt( pow(ht1+hpfc2,2) + dist12*dist12 ) ); //        if (lpg1.or.lmpg1)
         //     *    FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((ht1+hpfc2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }//        endif
   }//      endif
   if ( hpfc1+hpfc2 < foldxdist ) //      if(hpfc1+hpfc2 .lt.FOLDXDIST) then
   {
      if ( lpfc1 && lpfc2 ) //        if (lpfc1.and.lpfc2)  then
      {
         //        if (xdebug) then
         //          write(7,'(a,1x,$)')'mpfc1,mpfc2'
         //        endif
         dist12 = BasicDistance::g123dist(mpfc1,mpfc2);//        dist12 = g123dist(mpfc1,mpfc2)
         foldxdist = std::min(foldxdist, sqrt(hpfc1+hpfc2) + dist12*dist12 );//        FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((hpfc1+hpfc2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }//        endif
   }//      endif
   if ( hpfmc1+ht2 < foldxdist ) //      if(hpfmc1+ht2 .lt.FOLDXDIST) then
   {
      if ( lpfmc1 && (lpg2 || lmpg2) ) //        if (lpfmc1
      {
         //     *    .and.(lpg2.or.lmpg2)) then
         dist12 = foldxdist; //        dist12 = FoldxDist
         if ( lpg2  ) dist12 = std::min(dist12, BasicDistance::g123dist(pfmcase1,pg2));//        if (lpg2)
         //     *  dist12 = min(dist12,
         //     *    g123dist(pfmcase1,pg2))
         if ( lmpg2 ) dist12 = std::min(dist12, BasicDistance::g123dist(pfmcase1,mpg2)); //        if (lmpg2)
         //     *  dist12 = min(dist12,
         //     *    g123dist(pfmcase1,mpg2))
         if ( lpg2 || lmpg2 ) foldxdist = std::min(foldxdist, sqrt( pow(hpfmc1+ht2,2) + dist12*dist12 ) ); //        if (lpg2.or.lmpg2)
         //     *    FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((hpfmc1+ht2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }//        endif
   }//      endif
   if ( ht1+hpfmc2 < foldxdist ) //      if(ht1+hpfmc2 .lt.FOLDXDIST) then
   {
      if ( lpfmc2 && (lpg1 || lmpg1 ) ) //        if (lpfmc2
      {
         //     *    .and.(lpg1.or.lmpg1)) then
         dist12 = foldxdist; //        dist12 = FoldxDist
         if ( lpg1  ) dist12 = std::min(dist12,BasicDistance::g123dist(pg1,pfmcase2));//        if (lpg1)
         //     *    dist12 = min(dist12,
         //     *      g123dist(pg1,pfmcase2))
         if ( lmpg1 ) dist12 = std::min(dist12,BasicDistance::g123dist(mpg1,pfmcase2));//        if (lmpg1)
         //     *    dist12 = min(dist12,
         //     *      g123dist(mpg1,pfmcase2))
         if ( lpg1 || lmpg1 ) foldxdist = std::min(foldxdist, sqrt( pow(ht1+hpfmc2,2) + dist12*dist12 ) ); //        if (lpg1.or.lmpg1)
         //     *    FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((ht1+hpfmc2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }//        endif
   }//      endif
   if ( hpfmc1+hpfmc2 < foldxdist ) //      if(hpfmc1+hpfmc2 .lt.FOLDXDIST) then
   {
      if ( lpfmc1 && lpfmc2 ) //        if (lpfmc1.and.lpfmc2)  then
      {
         //        if (xdebug) then
         //          write(7,'(a,1x,$)')'pfmcase1,pfmcase2'
         //        endif
         dist12 = BasicDistance::g123dist(pfmcase1,pfmcase2); //        dist12 = g123dist(pfmcase1,pfmcase2)
         foldxdist = std::min(dist12,sqrt( pow(hpfmc1+hpfmc2,2) + dist12*dist12) );
         //        FoldxDist = MIN(FoldxDist,
         //     *    DSQRT((hpfmc1+hpfmc2)**2 +
         //     *    dist12**2))
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }//        endif
   }//      endif
//
//      if (xdebug) then
//          write(7,'(a,1x,6f9.2,1x,l1)')
//     *      "FoldxDist pg1: ",pg1,inncone(pg1)
//          write(7,'(a,1x,6f9.2,1x,l1)')
//     *      "FoldxDist mpg1: ",mpg1,inncone(mpg1)
//          write(7,'(a,1x,6f9.2,1x,l1)')
//     *      "FoldxDist pg2: ",pg2,inncone(pg2)
//          write(7,'(a,1x,6f9.2,1x,l1)')
//     *      "FoldxDist mpg2: ",mpg2,inncone(mpg2)
//          if (hpfc1.lt.9.D38)
//     *    write(7,'(a,1x,f12.4,1x,6f9.2,1x,l1)')
//     *      "FoldxDist pfcase1: ",
//     *      hpfc1,pfcase1,lpfc1
//          if (hpfc2.lt.9.D38)
//     *    write(7,'(a,1x,f12.4,1x,6f9.2,1x,l1)')
//     *      "FoldxDist hpfc1,pfcase2: ",
//     *      hpfc2,pfcase2,lpfc2
//          if (hpfmc1.lt.9.D38)
//     *    write(7,'(a,1x,f12.4,1x,6f9.2,1x,l1)')
//     *      "FoldxDist hfpmc1,pfmcase1: ",
//     *      hpfmc1,pfmcase1,lpfmc1
//          if (hpfmc2.lt.9.D38)
//     *    write(7,'(a,1x,f12.4,1x,6f9.2,1x,l1)')
//     *      "FoldxDist hfpmc2, pfmcase2: ",
//     *      hpfmc2,pfmcase2,lpfmc2
//          write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
//      endif
   return( foldxdist ); //        RETURN
} //      end
//

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//-----------------------------------------------------------------------------
// Name: Set1000()
//       Utility function that performs actions that were needed repeated in SOME
//       parts of FoldxDist
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void MakeGaol::Set1000( const arma::vec6& fcmp, const BoundaryPolytope& bt, const double ht, 
   arma::vec6& pfmcase, double& hpfmc, bool& lpfmc )
{
   pfmcase = bt.GetProjector() * fcmp;              //      call rmv6(fcmp,prj(,lp),pfmcase1)
   hpfmc   = BasicDistance::g456dist(pfmcase,fcmp); //      hpfmc1 = g456dist(pfmcase,fcmp1)
   hpfmc   = sqrt( pow(hpfmc,2) + pow(ht,2) );      //      hpfmc1 = sqrt(hpfmc1**2+ht1**2)
   lpfmc   = inncone(pfmcase);                      //      lpfmc1=inncone(pfmcase1)
//      if (xdebug) then
//        write(7,'(a,f12.4,f12.4,1x,6f9.2,1x,6f9.2)')
//     *   'hpfmc1,ht1,fcmp1,pfmcase1',
//     *   hpfmc1,ht1,fcmp1,pfmcase1
//      endif
//
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ######################################################################################
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// ######################################################################################
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
