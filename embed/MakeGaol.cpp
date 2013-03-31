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

   GenSeeds( gred, g6ErrorBox, seeds, seeddist, goodseed );

   tree.clear( );

   double gerr = ratio * arma::norm( gred, 2 );
   const double dmin = std::min( 0.05, 0.1*gerr );
   const double grmin = gred[0]*gred[0] + gred[1]*gred[1] + gred[2]*gred[2]; 

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

   nv = 0; // start from zero for C-style indexing

   tree.insert( gred );
   v[nv] = gred;
   int nexamined = nv;
   ivb[nv] = 0;
   vdist[nv] = 0.0;

   //     Add the good seeds that are not duplicates to
   //     the tree, and for each one apply the matrix
   //     for that boundary and if the image is nearly
   //     reduced, add that as well.

   int nv1 = 0;
   //!!!!!!!!!!!!!!!!!!!!!! are the seeds in one-to-one with the boundaries - do we need to dimension them 16?
   for ( int i=1; i<=15; ++i ) //  do i = 1,15
   {
      if( goodseed[i] )
      {
         nv1 = nv + 1;
         arma::vec6 nearest;
         if ( ! tree.NearestNeighbor( dmin, nearest,  seeds[i] ) )
         {
            ++nv;
            tree.insert( seeds[i] ); 
            v[nv] = seeds[i];
            ivb[nv] = 1;
            vdist[nv] = seeddist[i];
            if ( vdist[nv] < SMALLVDIST ) vdist[nv] = 0.0;
         }

         const arma::vec6 vt = ms[i] * seeds[i];
         nv1 = nv + 1;
         int nred = Reducer::NearRed( vt, SMALLNEARREDDELTA );
         if ( nred > 0 )
         {
            const double disttemp = NCDist( vt, gred );

            arma::vec6 nearest1;
            if (disttemp<gerr && ! tree.NearestNeighbor( dmin, nearest1, vt ) )
            {
               ++nv;
               tree.insert( vt );
               v[nv] = vt;
               ivb[nv] = 1;
               vdist[nv] = disttemp;
               if ( vdist[nv]< SMALLVDIST ) vdist[nv] = 0.0;
            }
         }
         else // ELSE for if nred
         {
            nred = Reducer::Near2Red( vt, gerr, vtred, dred );
            const double disttemp = NCDist( vtred, gred );
            if (disttemp<gerr && ! tree.NearestNeighbor( dmin, nearest, vtred ) )
            {
               ++nv;
               tree.insert( vtred );
               v[nv] = vtred; 
               ivb[nv] = 1;
               vdist[nv] = disttemp;
               if ( vdist[nv] < SMALLVDIST ) vdist[nv] = 0.0;
            } 
         } // ENDIF // end if nred
      } // end if GoodSeed
   } // end for

   for ( int i=0; i<6; ++i ) 
      newge[i] = std::max( g6ErrorBox[i], dmin ); 

   //     Now we have seeds to examine from nexamined+1
   //     through NV.  In the course of doing so, NV may
   //     increase.

   while ( true )
   {
      ++nexamined;
      if ( nexamined > nv || nv > nvmax-31 )
      {
         printf(" NV in MKREFL %d \n", nv );
         nv = tree.size();
         return;
      }

      if ( ivb[nexamined] < 10 ) continue;

      for ( int i=0; i<6; ++i )
      {
         newge[i] = std::max( pow( g6ErrorBox[i]/1.414, ivb[nexamined] ),dmin );
      }

      gerr = arma::norm( newge, 2 );

      //     Generate seed points on each of the 15 boundaries
      GenSeeds( v[nexamined], newge, seeds, seeddist, goodseed );

      //     Add the good seeds that are not duplicates to
      //     the tree, and for each one apply the matrix
      //     for that boundary and if the image is nearly
      //     reduced, add that as well.

      //     Differs from the prior loop only in the distance
      //     calculation

      for ( int i=1; i<=15; ++i )
      {
         const double d = seeds[i][0]*seeds[i][0] + seeds[i][1]*seeds[i][1] + seeds[i][2]*seeds[i][2];
         if ( goodseed[i] && d<10.8*grmin ) 
         {
            nv1 = nv + 1; 
            arma::vec6 closest;
            const double dminTemp = dmin * pow(1.85,ivb[nexamined]);
            if ( ! tree.NearestNeighbor( dminTemp, closest, seeds[i] ) )
            {
               ++nv;
               tree.insert( seeds[i] );
               v[nv] = seeds[i]; 
               ivb[nv] = ivb[nexamined] + 1;
               vdist[nv] = sqrt( seeddist[i]*seeddist[i] + vdist[nexamined]*vdist[nexamined] );
               if ( vdist[nv] < SMALLVDIST ) vdist[nv] = 0.0;
               //            write(*,*)"NV, SEED DIST,IB",NV,VDIST(NV),IB
               //            call printg6('Store SEED ',SEEDS(1,i))
            }
            const arma::vec6 vt = ms[i] * seeds[i];
            nv1 = nv + 1;
            const bool nredtemp = Reducer::NearRed( vt, SMALLNEARREDDELTA );
            if ( nredtemp )
            {
               vdist[nv1] = NCDist( vt, gred );
               if ( vdist[nv1]<gerr && ! tree.NearestNeighbor( dmin*pow(1.85,ivb[nexamined]), closest, vt ) )
               {
                  ++nv;
                  tree.insert( vt );
                  v[nv] = vt;
                  ivb[nv] = ivb[nexamined] + 1;
                  if ( vdist[nv]<SMALLVDIST ) vdist[nv] = 0.0;
                  //            write(*,*)"NV, M*SEED DIST,IB",NV,VDIST(NV),IB
                  //            call printg6('Store M*SEED ',VT)
               }
            }
            else
            {
               const bool nred = Reducer::Near2Red( vt, gerr, vtred, dred );
               vdist[nv1] = NCDist( vtred, gred );
               if ( vdist[nv1]<gerr && ! tree.NearestNeighbor( dmin*pow(1.85,ivb[nexamined]), closest, vtred ) )
               {
                  ++nv;
                  tree.insert( vtred ); 
                  v[nv] = vtred;
                  ivb[nv] = 1;
                  if ( vdist[nv]<SMALLVDIST ) vdist[nv] = 0.0;
                  //            write(*,*)"NV, M*SEED ** DIST",NV,VDIST(NV)
                  //            call printg6('Store M*SEED ** ',VT)
               }
            }
            //           call printg6('Reject SEED ',SEEDS(1,i))
         }
      }
   }
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


      arma::mat66 MakeGaol::MKPerp( const arma::mat66& prj )
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


   void MakeGaol::PrintG6 ( const std::string& text, const arma::vec6& v )
   {
      printf( "%s  %f %f %f %f %f %f", text.c_str(), v[0],v[1],v[2],v[3],v[4],v[5] );
   }

   void MakeGaol::PrintM66_INT ( const std::string& text, const arma::mat66& v )
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
//C
//C     GENSEEDS TAKES AN ARBITRARY G6 VECTOR AND PROJECTS
//C     IT ONTO EACH OF THE 15 BOUNDARIES.  EVEN IF THE
//C     ORIGINAL VECTOR WAS REDUCED THE PROJECTED VECTOR
//C     MAY FAIL TO BE REDUCED, AND MAY EVEN FAIL TO BE
//C     NEARLY REDUCED.  IT WILL BE MARKED AS A GOOD SEED
//C     PROVIDED IT IS BOTH NEARLY REDUCED AND IS WITHIN
//C     3.5 TIMES THE ERROR BOX GE OF GVEC.
//C
//C----------------------------------------------------------------------C

//C     Compute the seeds and the seed distances

//C      write(*,*) "GENSEEDS GVEC, GE"
//C      write(*,*) GVEC
//C      write(*,*) GE

   double gerr = arma::norm( g6ErrorBox, 2 );
   double a;

   arma::vec6 vtemp = m_prjList[19].GetPerp() * gvec;
   const double boundary67 = arma::norm( vtemp,2 );
   vtemp = m_prjList[20].GetPerp() * gvec;
   const double boundary9A = arma::norm( vtemp, 2 );
   vtemp = m_prjList[21].GetPerp() * gvec;
   const double boundaryCD = arma::norm( vtemp, 2 );

   //C      call printg6('GENSEEDS GVEC ',GVEC)

   for ( int ip=1; ip<=15; ++ip )
   {
      seeds[ip] = m_prjList[ip].GetProjector() * gvec;
      vtemp = m_prjList[ip].GetPerp() * gvec;
      //C        write(*,'(A,i2/,6(6F9.3/))'),"PRJ",
      //C     *  ip,(PRJ(i,ip),i=1,36)
      //C        write(*,'(A,i2/,6(6F9.3/))'),"PRJPERP",
      //C     *  ip,(PRJPERP(i,ip),i=1,36)
      a = arma::norm( vtemp, 2 );
      const double testg4 = gvec[3]-1.0E-6*sqrt(gvec[1]*gvec[2]);
      const double testg5 = gvec[4]-1.0E-6*sqrt(gvec[0]*gvec[2]);
      const double testg6 = gvec[5]-1.0E-6*sqrt(gvec[0]*gvec[2]);
      const bool ipTestA = ip==6 || ip==7 || ip==9 || ip==10 || ip==12 || ip==13;
      const bool ipTestB = ip==8 || ip==11 || ip==14 || ip==15;

      if ( ( testg4*testg5*testg6 <= 0.0 && ipTestA ) ||
           ( testg4*testg5*testg6 >  0.0 && ipTestB ) )
      {
         if ( ip>=6 && ip<=8 )
         {
            seeds[ip] = m_prjList[16].GetProjector() * gvec;
            vtemp = m_prjList[16].GetPerp() * gvec;
         }
         else if ( ip>=9 && ip<=15 )
         {
            seeds[ip] = m_prjList[17].GetProjector() * gvec;
            vtemp = m_prjList[17].GetPerp() * gvec;

         }
         else
         {
            seeds[ip] = m_prjList[18].GetProjector() * gvec;
            vtemp = m_prjList[18].GetPerp() * gvec;
         }
         a = arma::norm( vtemp, 2 );
      }

           if ( ip== 6 && gvec[4]<=gvec[5] && gvec[4]>0.0 ) a = boundary67;
      else if ( ip== 7 && gvec[4]>=gvec[5] && gvec[5]>0.0 ) a = boundary67;
      else if ( ip== 9 && gvec[3]<=gvec[5] && gvec[3]>0.0 ) a = boundary9A;
      else if ( ip==10 && gvec[3]>=gvec[5] && gvec[5]>0.0 ) a = boundary9A;
      else if ( ip==12 && gvec[3]<=gvec[4] && gvec[3]>0.0 ) a = boundaryCD;
      else if ( ip==13 && gvec[3]>=gvec[4] && gvec[4]>0.0 ) a = boundaryCD;

      //C       get distance and mark as bad if outside
      //C       3.5 times the errorbox

      seeddist[ip] = a;
      goodseed[ip] = Reducer::NearRed( seeds[ip], SMALLNEARREDDELTA );
      if ( ! goodseed[ip] )
      {
         double dred;
         goodseed[ip] = Reducer::Near2Red( seeds[ip], gerr, vtemp, dred );
         seeds[ip] = vtemp;
      }
      for ( int i=0; i<6; ++i )
      {
         goodseed[ip] = goodseed[ip] && fabs(gvec[i])-seeds[ip][i] <= 3.5*g6ErrorBox[i];


      }
      //C        write(seedlab,'(''BDRY'', I2)')IP
      //C        call printg6(seedlab,seeds(1,ip))
   } // loop over 15 boundaries //      enddo

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
// Name: BDCoord()
//C     Compute XS (distances to boundary sets) and
//C             YS (signed distances along boundaries)
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void MakeGaol::BDCoord( const arma::vec6& gvec, arma::vec6& xs, arma::vec6& ys )
{

   //     The boundaries are in the order
   //       1, 2, 678, 9AB, CDE, F
   const double sqrt2( sqrt(2.0) );
   const double sqrt5( sqrt(5.0) );
   xs[0] = (gvec[1] -     gvec[0]) /sqrt2;
   xs[1] = (gvec[2] -     gvec[1]) /sqrt2;
   xs[2] = (gvec[1] - fabs(gvec[3]))/sqrt2;
   xs[3] = (gvec[0] - fabs(gvec[4]))/sqrt2;
   xs[4] = (gvec[0] - fabs(gvec[5]))/sqrt2;
   xs[5] = (gvec[0]+gvec[1] + gvec[3] + gvec[4] + gvec[5])/sqrt5;

   ys[0] = (fabs(gvec[3]) - fabs(gvec[4]))/sqrt2;
   ys[1] = (fabs(gvec[4]) - fabs(gvec[5]))/sqrt2;
   ys[2] = (    gvec[4]  - fabs(gvec[5])/2.0)/sqrt2;
   ys[3] = (    gvec[3]  - fabs(gvec[5])/2.0)/sqrt2;
   ys[4] = (    gvec[3]  - fabs(gvec[4])/2.0)/sqrt2;
   ys[5] = (    gvec[1]  - gvec[0] + gvec[3] - gvec[4]) / 2.0;
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
inline void MKGReporter( const std::string& str, const double ncdist, const arma::vec6& gv1, const arma::vec6& gv2 )
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
   // if (xdebug) then
   //   write(7,'(a,1x,6f9.2)') "NCDIST gv1 ",gv1
   //   write(7,'(a,1x,6f9.2)') "NCDIST gv2 ",gv2
   // endif
   double ncdist( 9.0E38);
   double ncdistPrevious(ncdist);
   std::string str;
   ncdist = FoldMDist( gv1, gv2, 0, ncdist );
   str = std::string("ncdist1"); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 7, ncdist );
   str = std::string("ncdist2") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 10, ncdist );
   str = std::string("ncdist3") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 13, ncdist );
   str = std::string("ncdist4") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 15, ncdist );
   str = std::string("ncdist5") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 23, ncdist );
   str = std::string("ncdist6") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 24, ncdist );
   str = std::string("ncdist7") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif
   ncdist = FoldMDist( gv1, gv2, 25, ncdist );
   str = std::string("ncdist8") + ((ncdist!=ncdistPrevious) ? std::string(" **************************") : std::string("")); ncdistPrevious = ncdist;
   MKGReporter( str, ncdist, gv1, gv2 );
   // if (xdebug) then
   // write(7,'(a,1x,6f13.4)') "NCDIST dist down to ",
   //* NCDIST
   // endif



   return( ncdist );
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
// Name: DistF()
//C     Compute the distance to the F boundary
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double MakeGaol::DistF( const arma::vec6& g ) const
{
   const double x12( fabs(g[0])+fabs(g[1])+fabs(g[2]) - std::max(fabs(g[0]),std::max(fabs(g[1]),fabs(g[2]) ) ) );

   const double min1 = std::min( fabs(x12+g[3]+g[4]+g[5]),
                                 fabs(x12+g[3]-g[4]-g[5]));
   const double min2 = std::min( fabs(x12-g[3]+g[4]-g[5]),
                                 fabs(x12-g[3]-g[4]+g[5]));

   const double distftemp = std::min( min1, min2 ) / sqrt(5.0);
      return( distftemp );
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
// Name: inncone()
// Description: Determines whether a vector is within the Niggli cone
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool MakeGaol::inncone( const arma::vec6& gvec )
{
   int imax;
   const double precn( 1.005 );
   bool innconetemp( false );
   const double glow( std::min( gvec[0], std::min(gvec[1], gvec[2]) ) );
   if ( glow <= 0.0 ) return( false );

   const double ghigh( std::max( gvec[0], std::max(gvec[1], gvec[2]) ) );
   const double gmid( gvec[0]+gvec[1]+gvec[2] -glow -ghigh );
   if ( std::max( fabs(gvec[3]), std::max(fabs(gvec[4]),fabs(gvec[5]))) <= gmid*precn )
   {
     double gmax = fabs(gvec[3]);
     imax = 3;
     if (fabs(gvec[4]) > gmax)
     {
        gmax = fabs(gvec[4]);
        imax = 4;
     }
     if (fabs(gvec[5]) > gmax)
     {
        gmax = fabs(gvec[5]);
        imax = 5;
     }

     if ( imax!=3 && fabs(gvec[3])>glow*precn ) return( false );
     if ( imax!=4 && fabs(gvec[4])>glow*precn ) return( false );
     if ( imax!=5 && fabs(gvec[5])>glow*precn ) return( false );
     innconetemp = true;
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


//-----------------------------------------------------------------------------
// Name: FoldMDist()
//       return the distance betweem two G6 vectors
//       going through the given boundary
//       (7, A, C, F) after taking the second
//       vector through the 6 exchanges of A, B, C
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double MakeGaol::FoldMDist( const arma::vec6& gvec1, const arma::vec6& gvec2, const int ip, const double cFoldDist ) const
{
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

   const arma::vec6 m1g1     = m_prjList[ 1].GetTransform( ) * gvec1;
   const arma::vec6 m2g1     = m_prjList[ 2].GetTransform( ) * gvec1;
   const arma::vec6 m1m2g1   = m_prjList[16].GetTransform( ) * gvec1;
   const arma::vec6 m2m1g1   = m_prjList[17].GetTransform( ) * gvec1;
   const arma::vec6 m2m1m2g1 = m_prjList[18].GetTransform( ) * gvec1;
   const arma::vec6 m1g2     = m_prjList[ 1].GetTransform( ) * gvec2;
   const arma::vec6 m2g2     = m_prjList[ 2].GetTransform( ) * gvec2;
   const arma::vec6 m1m2g2   = m_prjList[16].GetTransform( ) * gvec2;
   const arma::vec6 m2m1g2   = m_prjList[17].GetTransform( ) * gvec2;
   const arma::vec6 m2m1m2g2 = m_prjList[18].GetTransform( ) * gvec2;

   double foldmdist = cFoldDist;

   if ( ip == 0 )
   {
      foldmdist = std::min( foldmdist, BasicDistance::g123dist( gvec1, gvec2 ) );
      return( foldmdist );
   }

   else if ( ip>=6 && ip<=8 )
   {
      h1       = fabs(gvec1[1]    - fabs(gvec1[3]))   /sqrt(2.0);
      h2       = fabs(gvec2[1]    - fabs(gvec2[3]))   /sqrt(2.0);
      h1m1     = fabs(m1g1[1]     - fabs(m1g1[3]))    /sqrt(2.0);
      h2m1     = fabs(m1g2[1]     - fabs(m1g2[3]))    /sqrt(2.0);
      h1m2     = fabs(m2g1[1]     - fabs(m2g1[3]))    /sqrt(2.0);
      h2m2     = fabs(m2g2[1]     - fabs(m2g2[3]))    /sqrt(2.0);
      h1m1m2   = fabs(m1m2g1[1]   - fabs(m1m2g1[3]))  /sqrt(2.0);
      h2m1m2   = fabs(m1m2g2[1]   - fabs(m1m2g2[3]))  /sqrt(2.0);
      h1m2m1   = fabs(m2m1g1[1]   - fabs(m2m1g1[3]))  /sqrt(2.0);
      h2m2m1   = fabs(m2m1g2[1]   - fabs(m2m1g2[3]))  /sqrt(2.0);
      h1m2m1m2 = fabs(m2m1m2g1[1] - fabs(m2m1m2g1[3]))/sqrt(2.0);
      h2m2m1m2 = fabs(m2m1m2g2[1] - fabs(m2m1m2g2[3]))/sqrt(2.0);
   }

   else if ( ip>=9 && ip<=11 )
   {
      h1       = fabs(gvec1[0]    - fabs(gvec1[4]))   /sqrt(2.0);
      h2       = fabs(gvec2[0]    - fabs(gvec2[4]))   /sqrt(2.0);
      h1m1     = fabs(m1g1[0]     - fabs(m1g1[4]))    /sqrt(2.0);
      h2m1     = fabs(m1g2[0]     - fabs(m1g2[4]))    /sqrt(2.0);
      h1m2     = fabs(m2g1[0]     - fabs(m2g1[4]))    /sqrt(2.0);
      h2m2     = fabs(m2g2[0]     - fabs(m2g2[4]))    /sqrt(2.0);
      h1m1m2   = fabs(m1m2g1[0]   - fabs(m1m2g1[4]))  /sqrt(2.0);
      h2m1m2   = fabs(m1m2g2[0]   - fabs(m1m2g2[4]))  /sqrt(2.0);
      h1m2m1   = fabs(m2m1g1[0]   - fabs(m2m1g1[4]))  /sqrt(2.0);
      h2m2m1   = fabs(m2m1g2[0]   - fabs(m2m1g2[4]))  /sqrt(2.0);
      h1m2m1m2 = fabs(m2m1m2g1[0] - fabs(m2m1m2g1[4]))/sqrt(2.0);
      h2m2m1m2 = fabs(m2m1m2g2[0] - fabs(m2m1m2g2[4]))/sqrt(2.0);
   }

   else if ( ip>=12 && ip<=13 )
   {
      h1       = fabs(gvec1[0]   - fabs(gvec1[5]   ))/sqrt(2.0);
      h2       = fabs(gvec2[0]   - fabs(gvec2[5]   ))/sqrt(2.0);
      h1m1     = fabs(m1g1[0]    - fabs(m1g1[5]    ))/sqrt(2.0);
      h2m1     = fabs(m1g2[0]    - fabs(m1g2[5]    ))/sqrt(2.0);
      h1m2     = fabs(m2g1[0]    - fabs(m2g1[5]    ))/sqrt(2.0);
      h2m2     = fabs(m2g2[0]    - fabs(m2g2[5]    ))/sqrt(2.0);
      h1m1m2   = fabs(m1m2g1[0]  - fabs(m1m2g1[5]  ))/sqrt(2.0);
      h2m1m2   = fabs(m1m2g2[0]  - fabs(m1m2g2[5]  ))/sqrt(2.0);
      h1m2m1   = fabs(m2m1g1[0]  - fabs(m2m1g1[5]  ))/sqrt(2.0);
      h2m2m1   = fabs(m2m1g2[0]  - fabs(m2m1g2[5]  ))/sqrt(2.0);
      h1m2m1m2 = fabs(m2m1m2g1[0]- fabs(m2m1m2g1[5]))/sqrt(2.0);
      h2m2m1m2 = fabs(m2m1m2g2[0]- fabs(m2m1m2g2[5]))/sqrt(2.0);
   }


   else if ( ip == 15 )
   {
      h1       = DistF(gvec1);
      h2       = DistF(gvec2);
      h1m1     = DistF(m1g1);
      h2m1     = DistF(m1g2);
      h1m2     = DistF(m2g1);
      h2m2     = DistF(m2g2);
      h1m1m2   = DistF(m1m2g1);
      h2m1m2   = DistF(m1m2g2);
      h1m2m1   = DistF(m2m1g1);
      h2m2m1   = DistF(m2m1g2);
      h1m2m1m2 = DistF(m2m1m2g1);
      h2m2m1m2 = DistF(m2m1m2g2);
   }

//C 8F boundary
   else if ( ip == 23 )
   {
      h1       = sqrt(pow(gvec1[1]   -fabs(gvec1[3]),2)   /2.0 + pow(gvec1[0]   -fabs(gvec1[4])   - fabs(gvec1[5]),2)   /3.0);
      h2       = sqrt(pow(gvec2[1]   -fabs(gvec2[3]),2)   /2.0 + pow(gvec2[0]   -fabs(gvec2[4])   - fabs(gvec2[5]),2)   /3.0);
      h1m1     = sqrt(pow(m1g1[1]    -fabs(m1g1[3]),2)    /2.0 + pow(m1g1[0]    -fabs(m1g1[4])    - fabs(m1g1[5]),2)    /3.0);
      h2m1     = sqrt(pow(m1g2[1]    -fabs(m1g2[3]),2)    /2.0 + pow(m1g2[0]    -fabs(m1g2[4])    - fabs(m1g2[5]),2)    /3.0);
      h1m2     = sqrt(pow(m2g1[1]    -fabs(m2g1[3]),2)    /2.0 + pow(m2g1[0]    -fabs(m2g1[4])    - fabs(m2g1[5]),2)    /3.0);
      h2m2     = sqrt(pow(m2g2[1]    -fabs(m2g2[3]),2)    /2.0 + pow(m2g2[0]    -fabs(m2g2[4])    - fabs(m2g2[5]),2)    /3.0);
      h1m1m2   = sqrt(pow(m1m2g1[1]  -fabs(m1m2g1[3]),2)  /2.0 + pow(m1m2g1[0]  -fabs(m1m2g1[4])  - fabs(m1m2g1[5]),2)  /3.0);
      h2m1m2   = sqrt(pow(m1m2g2[1]  -fabs(m1m2g2[3]),2)  /2.0 + pow(m1m2g2[0]  -fabs(m1m2g2[4])  - fabs(m1m2g2[5]),2)  /3.0);
      h1m2m1   = sqrt(pow(m2m1g1[1]  -fabs(m2m1g1[3]),2)  /2.0 + pow(m2m1g1[0]  -fabs(m2m1g1[4])  - fabs(m2m1g1[5]),2)  /3.0);
      h2m2m1   = sqrt(pow(m2m1g2[1]  -fabs(m2m1g2[3]),2)  /2.0 + pow(m2m1g2[0]  -fabs(m2m1g2[4])  - fabs(m2m1g2[5]),2)  /3.0);
      h1m2m1m2 = sqrt(pow(m2m1m2g1[1]-fabs(m2m1m2g1[3]),2)/2.0 + pow(m2m1m2g1[0]-fabs(m2m1m2g1[4])- fabs(m2m1m2g1[5]),2)/3.0);
      h2m2m1m2 = sqrt(pow(m2m1m2g2[1]-fabs(m2m1m2g2[3]),2)/2.0 + pow(m2m1m2g2[0]-fabs(m2m1m2g2[4])- fabs(m2m1m2g2[5]),2)/3.0);
   }

//C BF boundary
   else if ( ip == 24 )
   {
   h1       = sqrt( pow(gvec1[0]   -fabs(gvec1[4]),2)   /2.0 + pow(gvec1[1]   -fabs(gvec1[3])   -fabs(gvec1[5])   ,2) /3.0);
   h2       = sqrt( pow(gvec2[0]   -fabs(gvec2[4]),2)   /2.0 + pow(gvec2[1]   -fabs(gvec2[3])   -fabs(gvec2[5])   ,2) /3.0);
   h1m1     = sqrt( pow(m1g1[0]    -fabs(m1g1[4]),2)    /2.0 + pow(m1g1[1]    -fabs(m1g1[3])    -fabs(m1g1[5])    ,2) /3.0);
   h2m1     = sqrt( pow(m1g2[0]    -fabs(m1g2[4]),2)    /2.0 + pow(m1g2[1]    -fabs(m1g2[3])    -fabs(m1g2[5])    ,2) /3.0);
   h1m2     = sqrt( pow(m2g1[0]    -fabs(m2g1[4]),2)    /2.0 + pow(m2g1[1]    -fabs(m2g1[3])    -fabs(m2g1[5])    ,2) /3.0);
   h2m2     = sqrt( pow(m2g2[0]    -fabs(m2g2[4]),2)    /2.0 + pow(m2g2[1]    -fabs(m2g2[3])    -fabs(m2g2[5])    ,2) /3.0);
   h1m1m2   = sqrt( pow(m1m2g1[0]  -fabs(m1m2g1[4]),2)  /2.0 + pow(m1m2g1[1]  -fabs(m1m2g1[3])  -fabs(m1m2g1[5])  ,2) /3.0);
   h2m1m2   = sqrt( pow(m1m2g2[0]  -fabs(m1m2g2[4]),2)  /2.0 + pow(m1m2g2[1]  -fabs(m1m2g2[3])  -fabs(m1m2g2[5])  ,2) /3.0);
   h1m2m1   = sqrt( pow(m2m1g1[0]  -fabs(m2m1g1[4]),2)  /2.0 + pow(m2m1g1[1]  -fabs(m2m1g1[3])  -fabs(m2m1g1[5])  ,2) /3.0);
   h2m2m1   = sqrt( pow(m2m1g2[0]  -fabs(m2m1g2[4]),2)  /2.0 + pow(m2m1g2[1]  -fabs(m2m1g2[3])  -fabs(m2m1g2[5])  ,2) /3.0);
   h1m2m1m2 = sqrt( pow(m2m1m2g1[0]-fabs(m2m1m2g1[4]),2)/2.0 + pow(m2m1m2g1[1]-fabs(m2m1m2g1[3])-fabs(m2m1m2g1[5]),2) /3.0);
   h2m2m1m2 = sqrt( pow(m2m1m2g2[0]-fabs(m2m1m2g2[4]),2)/2.0 + pow(m2m1m2g2[1]-fabs(m2m1m2g2[3])-fabs(m2m1m2g2[5]),2) /3.0);
   }

//C EF boundary
   else if ( ip == 25 )
   {
      h1       = sqrt(pow((gvec1[0]   -fabs(gvec1[5])),2)   /2.0 + pow((gvec1[1]-fabs(gvec1[3])      -fabs(gvec1[4])),2)   /3.0);
      h2       = sqrt(pow((gvec2[0]   -fabs(gvec2[5])),2)   /2.0 + pow((gvec2[1]-fabs(gvec2[3])      -fabs(gvec2[4])),2)   /3.0);
      h1m1     = sqrt(pow((m1g1[0]    -fabs(m1g1[5])),2)    /2.0 + pow((m1g1[1]-fabs(m1g1[3])        -fabs(m1g1[4])),2)    /3.0);
      h2m1     = sqrt(pow((m1g2[0]    -fabs(m1g2[5])),2)    /2.0 + pow((m1g2[1]-fabs(m1g2[3])        -fabs(m1g2[4])),2)    /3.0);
      h1m2     = sqrt(pow((m2g1[0]    -fabs(m2g1[5])),2)    /2.0 + pow((m2g1[1]-fabs(m2g1[3])        -fabs(m2g1[4])),2)    /3.0);
      h2m2     = sqrt(pow((m2g2[0]    -fabs(m2g2[5])),2)    /2.0 + pow((m2g2[1]-fabs(m2g2[3])        -fabs(m2g2[4])),2)    /3.0);
      h1m1m2   = sqrt(pow((m1m2g1[0]  -fabs(m1m2g1[5])),2)  /2.0 + pow((m1m2g1[1]-fabs(m1m2g1[3])    -fabs(m1m2g1[4])),2)  /3.0);
      h2m1m2   = sqrt(pow((m1m2g2[0]  -fabs(m1m2g2[5])),2)  /2.0 + pow((m1m2g2[1]-fabs(m1m2g2[3])    -fabs(m1m2g2[4])),2)  /3.0);
      h1m2m1   = sqrt(pow((m2m1g1[0]  -fabs(m2m1g1[5])),2)  /2.0 + pow((m2m1g1[1]-fabs(m2m1g1[3])    -fabs(m2m1g1[4])),2)  /3.0);
      h2m2m1   = sqrt(pow((m2m1g2[0]  -fabs(m2m1g2[5])),2)  /2.0 + pow((m2m1g2[1]-fabs(m2m1g2[3])    -fabs(m2m1g2[4])),2)  /3.0);
      h1m2m1m2 = sqrt(pow((m2m1m2g1[0]-fabs(m2m1m2g1[5])),2)/2.0 + pow((m2m1m2g1[1]-fabs(m2m1m2g1[3])-fabs(m2m1m2g1[4])),2)/3.0);
      h2m2m1m2 = sqrt(pow((m2m1m2g2[0]-fabs(m2m1m2g2[5])),2)/2.0 + pow((m2m1m2g2[1]-fabs(m2m1m2g2[3])-fabs(m2m1m2g2[4])),2)/3.0);
   }


   else if ( ip == 22 )
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

//-----------------------------------------------------------------------------
// Name: FoldxDist()
//       return the distance between two
//       G6 vectors going through the given boundary
//       (7, A, C or F)

//       is1 and is2 are the boundary symmetry ops used
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double MakeGaol::FoldxDist( const arma::vec6& gvec1, const arma::vec6& gvec2, const int ip, const double cFoldDist,
   const int is1DUMMY, const int is2DUMMY) const
{
   is1DUMMY, is2DUMMY; // just to keep the compiler happy
//      if (xdebug) then
//        write(7,*) "FoldxDist ip, CFOLDDIST: ",ip,CFOLDDIST
//        write(7,'(a,1x,6f9.2)') "FOLDXDIST gvec1: ",gvec1
//        write(7,'(a,1x,6f9.2)') "FOLDXDIST gvec2: ",gvec2
//      endif

   double ht1        (DBL_MAX);
   double ht2        (DBL_MAX);
   double hpfc1      (DBL_MAX);
   double hpfc2      (DBL_MAX);
   double hpfmc1     ( 0.0 );
   double hpfmc2     ( 0.0 );
   bool lpfc1        (false);
   bool lpfc2        (false);
   bool lpfmc1       (false);
   bool lpfmc2       (false);

   double foldxdist = cFoldDist;

   arma::vec6 pg1( arma::abs(gvec1) );
   arma::vec6 pg2( arma::abs(gvec2) );
   arma::vec6 mpg1( pg1 );
   arma::vec6 mpg2( pg2 );
   arma::vec6 fc1( gvec1 );
   arma::vec6 fc2( gvec2 );

   arma::vec6 pfcase1,ppfcase1;
   arma::vec6 fcmp1, fcmp2;
   arma::vec6 pfcase2,ppfcase2;
   arma::vec6 pfmcase1,pfmcase2;
   arma::vec6 ppfmcase1,ppfmcase2;

   bool lfc1( true );
   bool lfc2( true );
   int lp( ip );
   int kp;

   for ( int i=3; i<6; ++i )
   {
      if ( fc1[i]<= 0.0 ) lfc1 = ! lfc1;
      if ( fc2[i]<= 0.0 ) lfc2 = ! lfc2;
   }

   //C     Process 8F boundary
   if ( ip==23 )
   {
      if ( lfc1 )
      {
         fc1[3] = -fc1[3];
         if ( fabs(fc1[4])>=fabs(fc1[5]) )
            fc1[4] = -fc1[4];
         else
            fc1[5] = -fc1[5];
      }
      if ( lfc2 )
      {
         fc2[3] = -fc2[3];
         if ( fabs(fc2[4])>=fabs(fc2[5]) )
            fc2[4] = -fc2[4];
         else
            fc2[5] = -fc2[5];
      }

      pg1[0]  = (2.0*fc1[0]-fc1[4]-fc1[5])/3.0;
      pg2[0]  = (2.0*fc2[0]-fc2[4]-fc2[5])/3.0;
      mpg1[0] = pg1[0];
      mpg2[0] = pg2[0];
      pg1[1]  = (fc1[1]-fc1[3])/2.0;
      pg2[1]  = (fc2[1]-fc2[3])/2.0;
      mpg1[1] = pg1[1];
      mpg2[1] = pg2[1];
      pg1[3]  = -pg1[1];
      pg2[3]  = -pg2[1];
      mpg1[3] = pg1[1];
      mpg2[3] = pg2[1];
      pg1[4]  = (2.0*fc1[4]-fc1[0]-    fc1[5])/3.0;
      pg2[4]  = (2.0*fc2[4]-fc2[0]-    fc2[5])/3.0;
      mpg1[4] = (2.0*fc1[0]-fc1[4]-    fc1[5])/3.0;
      mpg2[4] = (2.0*fc2[0]-fc2[4]-    fc2[5])/3.0;
      pg1[5]  = (2.0*fc1[5]-fc1[0]-    fc1[4])/3.0;
      pg2[5]  = (2.0*fc2[5]-fc2[0]-    fc2[4])/3.0;
      mpg1[5] =     (fc1[0]+fc1[4]-2.0*fc1[5])/3.0;
      mpg2[5] =     (fc2[0]+fc2[4]-2.0*fc2[5])/3.0;

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
   }


   //C     Process BF boundary
   else if ( ip==24 )
   {
      if ( lfc1 )
      {
         fc1[4] = -fc1[4];
         if ( fabs(fc1[3])>=fabs(fc1[5]) )
            fc1[3] = -fc1[3];
         else
            fc1[5] = -fc1[5];
      }
      if ( lfc2 )
      {
         fc2[4] = -fc2[4];
         if ( fabs(fc2[3])>=fabs(fc2[5]) )
            fc2[3] = -fc2[3];
         else
            fc2[5] = -fc2[5];
      }
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
   }

//C     Process EF boundary
   else if ( ip==25 )
   {
      if ( lfc1 )
      {
         fc1[5] = -fc1[5];
         if ( fabs(fc1[4])>=fabs(fc2[3]) )
            fc1[4] = -fc1[4];
         else
            fc1[3] = -fc1[3];
      }
      if( lfc2 )
      {
         fc2[5] = -fc2[5];
         if ( fabs(fc2[4])>=fabs(fc2[3]) )
            fc2[4] = -fc2[4];
         else
            fc2[3] = -fc2[3];
      }
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
   }


   else if ( ip==6 || ip==7 || ip==8 )
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
        ht1 = fabs(gvec1[1]-fabs(gvec1[3]))/sqrt(2.0);
        ht2 = fabs(gvec2[1]-fabs(gvec2[3]))/sqrt(2.0);
        if ( lfc1 )
        {
           fc1[3] = -fc1[3];
           if ( fabs(fc1[4])>=fabs(fc1[5]) )
              fc1[4] = -fc1[4];
           else
              fc1[5] = -fc1[5];
        }
        if ( lfc2 )
        {
           fc2[3] = -fc2[3];
           if ( fabs(fc2[4])>=fabs(fc2[5]) )
              fc2[4] = -fc2[4];
           else
              fc2[5] = -fc2[5];
        }
        ppfcase1 = m_prjList[23].GetPerp() * fc1;
        hpfc1 = arma::norm(ppfcase1,2);
        ppfcase2 = m_prjList[23].GetPerp() * fc2;
        hpfc2 = arma::norm(ppfcase2,2);
        pfcase1 = m_prjList[23].GetProjector() * fc1;
        pfcase2 = m_prjList[23].GetProjector() * fc2;
        lpfc1=inncone(pfcase1);
        lpfc2=inncone(pfcase2);
        kp = 8;
        lpfmc1 = true;
        lpfmc2 = true;
        for ( int i=0; i<6; ++i )
        {
           fcmp1[i] = mpg1[i];
           if ( fcmp1[i] <= 0.0 )  lpfmc1 = ! lpfmc1;
           fcmp2[i] = mpg2[i];
           if ( fcmp2[i] <= 0.0 ) lpfmc2 = ! lpfmc2;
        }
        for ( int i=3; i<6; ++i )
        {
           fcmp1[i] = -fabs(fcmp1[i]);
           fcmp2[i] = -fabs(fcmp2[i]);
        }
        if ( lpfmc1 )
        {
           if ( fcmp1[4] < fcmp1[5] )
              fcmp1[5] = -fcmp1[5];
           else
              fcmp1[4] = -fcmp1[4];
        }
        if ( lpfmc2 )
        {
           if ( fcmp2[4] < fcmp2[5] )
              fcmp2[5] = -fcmp2[5];
           else
              fcmp2[4] = -fcmp2[4];
        }
        lp = 23;
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
   }

   else if ( ip==10 || ip==9 || ip==11 )
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
      ht1 = fabs(gvec1[0]-fabs(gvec1[4]))/sqrt(2.0);
      ht2 = fabs(gvec2[0]-fabs(gvec2[4]))/sqrt(2.0);
      if ( lfc1 )
      {
         fc1[4] = -fc1[4];
         if ( fabs(fc1[3])>=fabs(fc1[5]) )
            fc1[3] = -fc1[3];
         else
            fc1[5] = -fc1[5];
      }
      if ( lfc2 )
      {
         fc2[4] = -fc2[4];
         if ( fabs(fc2[3])>=fabs(fc2[5]) )
            fc2[3] = -fc2[3];
         else
            fc2[5] = -fc2[5];
      }
      ppfcase1 = this->m_prjList[24].GetPerp() * fc1;
      hpfc1 = arma::norm( ppfcase1,2 );
      ppfcase2 = this->m_prjList[24].GetPerp() * fc2;
      hpfc2 = arma::norm( ppfcase2,2 );
      pfcase1 = m_prjList[24].GetProjector() * fc1;
      pfcase2 = m_prjList[24].GetProjector() * fc2;
      lpfc1 = inncone(pfcase1);
      lpfc2 = inncone(pfcase2);
      kp = 11;
      lp = 24;
      lpfmc1 = true;
      lpfmc2 = true;
      for( int i=0; i<6; ++i )
      {
         fcmp1[i] = mpg1[i];
         if ( fcmp1[i] <= 0.0 ) lpfmc1 = ! lpfmc1;
         fcmp2[i] = mpg2[i];
         if ( fcmp2[i] <= 0.0 ) lpfmc2 = ! lpfmc2;
      }
      for ( int i=3; i<6; ++i )
      {
         fcmp1[i] = -fabs(fcmp1[i]);
         fcmp2[i] = -fabs(fcmp2[i]);
      }
      if ( lpfmc1 )
      {
         if ( fcmp1[3] < fcmp1[5] )
            fcmp1[5] = -fcmp1[5];
         else
            fcmp1[3] = -fcmp1[3];
      }
      if ( lpfmc2 )
      {
         if ( fcmp2[3] < fcmp2[5] )
            fcmp2[5] = -fcmp2[5];
         else
            fcmp2[3] = -fcmp2[3];
      }
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
}
   else if ( ip==13 || ip==12 || ip==14 )
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
      ht1 = fabs(gvec1[0]-fabs(gvec1[5]))/sqrt(2.0);
      ht2 = fabs(gvec2[0]-fabs(gvec2[5]))/sqrt(2.0);
      if ( lfc1 )
      {
         fc1[5] = -fc1[5];
         if ( fabs(fc1[4]) >= fabs(fc1[3]) )
            fc1[4] = -fc1[4];
         else
            fc1[3] = -fc1[3];
      }
      if ( lfc2 )
      {
         fc2[5] = -fc2[5];
         if ( fabs(fc2[4]) >= fabs(fc2[3]) )
            fc2[4] = -fc2[4];
         else
            fc2[3] = -fc2[3];
      }
      ppfcase1 = m_prjList[25].GetPerp() * fc1;
      hpfc1 = arma::norm( ppfcase1,2 );
      ppfcase2 = m_prjList[25].GetPerp() * fc2;
      hpfc2 = arma::norm( ppfcase2,2 );
      pfcase1 = m_prjList[25].GetProjector() * fc1;
      pfcase2 = m_prjList[25].GetProjector() * fc2;
      lpfc1 = inncone( pfcase1 );
      lpfc2 = inncone( pfcase2 );
      kp = 14;
      lp = 25;
      lpfmc1 = true;
      lpfmc2 = true;
      for ( int i=0; i<6; ++i )
      {
         fcmp1[i] = mpg1[i];
         if ( fcmp1[i]<=0.0 ) lpfmc1 = ! lpfmc1;
         fcmp2[i] = mpg2[i];
         if ( fcmp2[i]<=0.0 ) lpfmc2 = ! lpfmc2;
      }
      for ( int i=3; i<6; ++i )
      {
         fcmp1[i] = -fabs(fcmp1[i]);
         fcmp2[i] = -fabs(fcmp2[i]);
      }
      if ( lpfmc1 )
      {
         if ( fcmp1[3] < fcmp1[4] )
            fcmp1[4] = -fcmp1[4];
         else
            fcmp1[3] = -fcmp1[3];
      }
      if ( lpfmc2 )
      {
         if ( fcmp2[3] < fcmp2[4] )
            fcmp2[4] = -fcmp2[4];
         else
            fcmp2[3] = -fcmp2[3];
      }
      Set1000( fcmp1, m_prjList[lp], ht1, pfmcase1, hpfmc1, lpfmc1 );
      Set1000( fcmp2, m_prjList[lp], ht2, pfmcase2, hpfmc2, lpfmc2 );
   }
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
   }

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


   double dist12;
   if ( ht1+ht2 < foldxdist )
   {
      dist12 = foldxdist;
      if ( lpg1 )
      {
         if ( lpg2  ) dist12 = std::min(dist12,BasicDistance::g123dist(pg1, pg2));
         if ( lmpg2 ) dist12 = std::min(dist12,BasicDistance::g123dist(pg1,mpg2));
      }
      if ( lmpg1 )
      {
         if ( lpg2  ) dist12 = std::min(dist12,BasicDistance::g123dist(mpg1, pg2));
         if ( lmpg2 ) dist12 = std::min(dist12,BasicDistance::g123dist(mpg1,mpg2));
      }
      foldxdist = std::min(foldxdist,sqrt( pow(ht1+ht2,2) + dist12*dist12) );
      //        if (xdebug)
      //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
   }
   if ( hpfc1+ht2 < foldxdist )
   {
      if ( lpfc1 && (lpg2 || lmpg2) )
      {
         dist12 = foldxdist; //        dist12 = FoldxDist
         if ( lpg2  ) dist12 = std::min(dist12,BasicDistance::g123dist(pfcase1,pg2));
         if ( lmpg2 ) dist12 = std::min(dist12,BasicDistance::g123dist(pfcase1,mpg2));
         if( lpg2 || lmpg2 ) foldxdist = std::min( foldxdist, sqrt( pow(hpfc1+ht2,2) + dist12*dist12 ) );
         //        if (xdebug)
         //     * write(7,'(a,1x,3f12.4,6f9.2)')
         //     *  "FoldxDist: hpfc1, ht2, dist12",
         //     *   hpfc1, ht2, dist12, FoldxDist
      }
   }
   if ( ht1+hpfc2 < foldxdist )
   {
      if ( lpfc2 && (lpg1 || lmpg1 ) )
      {
         dist12 = foldxdist;
         if ( lpg1  ) dist12 = std::min( dist12,BasicDistance::g123dist(pg1,pfcase2));
         if ( lmpg1 ) dist12 = std::min(dist12,BasicDistance::g123dist(mpg1,pfcase2));
         if ( lpg1 || lmpg1 ) foldxdist = std::min( foldxdist, sqrt( pow(ht1+hpfc2,2) + dist12*dist12) );
         //        if (xdebug)
         //     * write(7,'(a,1x,3f12.4,6f9.2)')
         //     *  "FoldxDist: ht1, hpfc2, dist12",
         //     *   ht1, hpfc2, dist12, FoldxDist
      }
   }

   arma::vec6 mpfc1, mpfc2;
   if ( hpfc1+hpfc2 < foldxdist )
   {
      if ( lpfc1 && lpfc2 )
      {
         dist12 = BasicDistance::g123dist(pfcase1,pfcase2);
         foldxdist = std::min( foldxdist, sqrt( pow(hpfc1+hpfc2,2) + dist12*dist12 ) );
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }
   }
   if( lpfc1 && (hpfc1+std::min(ht2,hpfc2)<foldxdist ) )
   {
      mpfc1 = ms[kp] * pfcase1;
   }
   if ( lpfc2 && std::min(ht1,hpfc1)+hpfc2 < foldxdist )
   {
      mpfc2 = ms[kp] * pfcase2;
   }
   if ( hpfc1+ht2 < foldxdist )
   {
      if ( lpfc1 && (lpg2 || lmpg2 ) )
      {
         dist12 = foldxdist;
         if ( lpg2  ) dist12 = std::min( dist12,BasicDistance::g123dist(mpfc1,pg2));
         if ( lmpg2 ) dist12 = std::min( dist12,BasicDistance::g123dist(mpfc1,mpg2) );
         if ( lpg2 || lmpg2 ) foldxdist = std::min(foldxdist, sqrt( pow(hpfc1+ht2,2) + dist12*dist12 ) );
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }
   }
   if ( ht1+hpfc2 < foldxdist )
   {
      if ( lpfc2 && (lpg1 || lmpg1 ) )
      {
         //        if (xdebug) then
         //          write(7,'(a,1x,$)')'pg1,mpfc2'
         //        endif
         dist12 = foldxdist;
         if ( lpg1  ) dist12 = std::min( dist12,BasicDistance::g123dist(pg1,mpfc2));
         if ( lmpg1 ) dist12 = std::min( dist12,BasicDistance::g123dist(mpg1,mpfc2));
         if ( lpg1 || lmpg1 ) foldxdist = std::min(foldxdist,sqrt( pow(ht1+hpfc2,2) + dist12*dist12 ) );
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }
   }
   if ( hpfc1+hpfc2 < foldxdist )
   {
      if ( lpfc1 && lpfc2 )
      {
         //        if (xdebug) then
         //          write(7,'(a,1x,$)')'mpfc1,mpfc2'
         //        endif
         dist12 = BasicDistance::g123dist(mpfc1,mpfc2);
         foldxdist = std::min(foldxdist, sqrt(hpfc1+hpfc2) + dist12*dist12 );
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }
   }
   if ( hpfmc1+ht2 < foldxdist )
   {
      if ( lpfmc1 && (lpg2 || lmpg2) )
      {
         dist12 = foldxdist;
         if ( lpg2  ) dist12 = std::min(dist12, BasicDistance::g123dist(pfmcase1,pg2));
         if ( lmpg2 ) dist12 = std::min(dist12, BasicDistance::g123dist(pfmcase1,mpg2));
         if ( lpg2 || lmpg2 ) foldxdist = std::min(foldxdist, sqrt( pow(hpfmc1+ht2,2) + dist12*dist12 ) );
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }
   }
   if ( ht1+hpfmc2 < foldxdist ) //      if(ht1+hpfmc2 .lt.FOLDXDIST) then
   {
      if ( lpfmc2 && (lpg1 || lmpg1 ) ) //        if (lpfmc2
      {
         dist12 = foldxdist; //        dist12 = FoldxDist
         if ( lpg1  ) dist12 = std::min(dist12,BasicDistance::g123dist(pg1,pfmcase2));
         if ( lmpg1 ) dist12 = std::min(dist12,BasicDistance::g123dist(mpg1,pfmcase2));
         if ( lpg1 || lmpg1 ) foldxdist = std::min(foldxdist, sqrt( pow(ht1+hpfmc2,2) + dist12*dist12 ) );
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }
   }
   if ( hpfmc1+hpfmc2 < foldxdist )
   {
      if ( lpfmc1 && lpfmc2 )
      {
         //        if (xdebug) then
         //          write(7,'(a,1x,$)')'pfmcase1,pfmcase2'
         //        endif
         dist12 = BasicDistance::g123dist(pfmcase1,pfmcase2);
         foldxdist = std::min(dist12,sqrt( pow(hpfmc1+hpfmc2,2) + dist12*dist12) );
         //        if (xdebug)
         //     * write(7,'(a,1x,6f9.2)') "FoldxDist: ",FoldxDist
      }
   }

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
   return( foldxdist );
} //      end FoldxDust


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
   pfmcase = bt.GetProjector() * fcmp;
   hpfmc   = BasicDistance::g456dist(pfmcase,fcmp);
   hpfmc   = sqrt( pow(hpfmc,2) + pow(ht,2) );
   lpfmc   = inncone(pfmcase);
//      if (xdebug) then
//        write(7,'(a,f12.4,f12.4,1x,6f9.2,1x,6f9.2)')
//     *   'hpfmc1,ht1,fcmp1,pfmcase1',
//     *   hpfmc1,ht1,fcmp1,pfmcase1
//      endif

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
