/*
This is a main program for testing the C++ translation of some
of the cell reduction/lattice identification code.
*/

/*


a boundary polytope will need
   the projector
   the normal
   the boundary transform
   the "names:
    AB designation
    Niggli/Roof designation of exists
    sometimes list of included Niggli/Roof types
    same for IT types
    ? special positions
    text description of the boundary

Niggli/Roof classes
International table classes

Reduction matrices?
*/

/*lattice object members
    degrees of freedom
    projector (6-d) (including normalizer if integers)
    niggli notation
    IT notation
    lattice type (i.e. cI)
    



Lattice Tools
    matrix to remove centering (3-d and 6-d)
    matrix to re-center (3-d and 6-d)
    
    
Boundary info
    designation
    degrees of freedom
    projector (6-d) [integers in FillPrjList in BasicReduction.for, doubles inBLOCK DATA PRJS in MKGAOL.FOR]
    ?? Boundary transformations (ask hjb)
    

    */
#define USE_LOCAL_HEADERS

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>
#include <cassert>
#include <iostream>
#include <istream>
#include <ostream>
#include <string>

#include "BoundaryPolytope.h"
#include "BoundaryPolytopeList.h"
#include "NiggliPolytope.h"
#include "NiggliPolytopeList.h"
#include "Cell.h"
#include "Reducer.h"
#include "V7.h"
#include "Vector_3d.h"
#include "NCDist.h"
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void test1( void )
{
   //const arma::vec6 g1("10.1 10.2 10.3 10.4 10.5 10.6");
   const arma::vec6 g1("100.1 100.2 100.3 98 99 100");
   const Cell c1( g1 );
   const Cell c1i( c1.Inverse( ) );
   const Cell c1ii(c1i.Inverse( ) );
   const arma::vec6 v1ii( c1ii.Cell2V6() );
   const arma::vec6 vdiff( g1-v1ii );
   //const int i19191 = 19191;
}

static int iseed = 19191;

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
arma::vec6 GenRandG6( void )
{
   arma::vec6 g;
   arma::vec6 rands;

   for ( int i=0; i<6; ++i )
   {
      rands[i] = randfg( iseed );
   }

   rands /= arma::norm( rands, 2 );

   return( rands );
}

void NiggliClassIdentifier( void )
{
   arma::vec6 reducedBase;
   arma::mat66 m;
   m.eye( );
   // const arma::vec6 base("9.997321 10.003294 10.000253 0.001003 9.996918 9.997098");
   const arma::vec6 base("10 10 10   10 10 10");
   Reducer::Reduce( base, m, reducedBase, 0.0 );



   //     exit(0);

   const int nfollow=4;
   for ( int jfollow=0; jfollow<nfollow; ++jfollow )
   {
      // Generate the random step away from the universal starting point
      const arma::vec6 rands = GenRandG6( );
//      if( jfollow != 13 ) continue;

      // take the random step
      const arma::vec6 testBase( base + 0.01 * rands );
      Reducer::Reduce( testBase, m, reducedBase, 0.0 );

      // ERROR-CHECK -if the volume changed, the tranform was invalid !!!!!!!
      if ( fabs(Cell(testBase).Volume( ) - Cell(reducedBase).Volume( )) > 0.01 )
      {
         const double v1 = Cell(testBase).Volume( );
         const double v2 = Cell(reducedBase).Volume( );
         const int i19191 = 19191;
      }
      const V7 v7ReducedBase( reducedBase );
         printf( "\n\n\nV7Base %f %f %f %f %f %f %f\n", v7ReducedBase[0], v7ReducedBase[1], v7ReducedBase[2], v7ReducedBase[3], v7ReducedBase[4], v7ReducedBase[5], v7ReducedBase[6] );


      NiggliPolytope np;
      NiggliPolytopeList npl;
      double previousV7Dist  = 0.0;
      const int nstep(100);

      printf( "jfollow %d   \n", jfollow );
      arma::mat66 mtemp;
      arma::vec6 vprojReduced;
      arma::vec6 reducedTest;

      // Loop to run from the randomized starting point to the reduced
      // version of the randomized starting point
      for ( int i=0; i<=nstep; ++i )
      {
         const double di(i/double(nstep));
         // increment along the path
         const arma::vec6 test( di*reducedBase + (1.0-di)*testBase );
         // Reduce the current probe lattice
         Reducer::Reduce( test, m, reducedTest, 0.0 );

         // loop over the 44 Niggli lattice types
         for ( size_t iniggli=1; iniggli<npl.size(); ++ iniggli )
         {
            // project the reduced test lattice onto the current boundary
            const arma::vec6 vprojected( npl[iniggli].GetProjector() * reducedTest );
            // make sure the projection is reduced
            const bool bRed = Reducer::Reduce( vprojected, m, vprojReduced, 0.0 );
            if ( bRed )
            {
               int in;
               double vp[6], rT[6];
               for (in = 0; in < 6; in++) {
                 vp[in] = vprojected[in];
                 rT[in] = reducedTest[in];
               }
               const double dist = NCDist(vp,rT);
               printf( "%f ", dist );
            }
            else
            {
               // the projected G6 point is not a valid lattice in 3-D
               printf( "-1 " );
            }
         }
         printf( "\n" );
      }
      printf( "\n\njfollow= %d\n\n", jfollow );
   }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void Follower( void )
{
   const int N_FOLLOW_STEPS = 4;
   double distance[1000+20];
   std::string buffer[10000];
   char cBuffer[100];
   arma::vec6 reducedBase;
   arma::mat66 m;
   m.eye( );
   const arma::vec6 base("5 10 12  0 0 0"); // BF
   Reducer::Reduce( base, m, reducedBase, 0.0 );




   const int nfollow=N_FOLLOW_STEPS;
   int lines;
   for ( int jfollow=0; jfollow<nfollow; ++jfollow )
   {
      lines = 0;
      arma::vec6 rands = GenRandG6( );
      rands[0] = 0;
      rands[1] = 0;
      rands[2] = 0;

      const arma::vec6 testBase( base + 0.01 * rands );
      Reducer::Reduce( testBase, m, reducedBase, 0.0 );

      if ( fabs(Cell(testBase).Volume( ) - Cell(reducedBase).Volume( )) > 0.01 )
      {
         const double v1 = Cell(testBase).Volume( );
         const double v2 = Cell(reducedBase).Volume( );
         const int i19191 = 19191;
      }
      const V7 v7ReducedBase( reducedBase );
      if ( v7ReducedBase[0] < 1.0 ) 
      {
         const double a = sqrt( testBase[0] );
         const double b = sqrt( testBase[1] );
         const double c = sqrt( testBase[2] );
         const double alpha = acos(testBase[3]/2.0/b/c)*180.0/4.0/atan(1.0);
         const double beta  = acos(testBase[4]/2.0/a/c)*180.0/4.0/atan(1.0);
         const double gamma = acos(testBase[5]/2.0/a/b)*180.0/4.0/atan(1.0);
         const int i19191 = 19191;
      }
         sprintf( cBuffer, "\n\nV7Base %f %f %f %f %f %f %f\n", v7ReducedBase[0], v7ReducedBase[1], v7ReducedBase[2], v7ReducedBase[3], v7ReducedBase[4], v7ReducedBase[5], v7ReducedBase[6] );
         buffer[lines] = cBuffer;
         ++lines;

      sprintf( cBuffer, "jfollow= %d\n\n", jfollow );
      buffer[lines] = cBuffer;
      ++lines;

      double previousV7Dist  = 0.0;
      const int nstep(100);
      bool bShouldPrint = false;
      double ncdistMax = 0.0;

      for ( int i=0; i<=nstep; ++i )
      {
         const double di(i/double(nstep));
         const arma::vec6 test( di*reducedBase + (1.0-di)*testBase );
         arma::vec6 reducedTest;
         Reducer::Reduce( test, m, reducedTest, 0.0 );
         const V7         v7ReducedTest( reducedTest );
         const double     v7Dist       ( (v7ReducedTest-v7ReducedBase).Norm() );
         double    dvncdist;
         const arma::vec6 v6diff       ( reducedTest - reducedBase );
         const double     v6diffNorm   ( arma::norm( v6diff,2 ) );
         {
               int in;
               double rB[6], rT[6];
               for (in = 0; in < 6; in++) {
                 rB[in] = reducedBase[in];
                 rT[in] = reducedTest[in];
               }
               dvncdist = NCDist(rT,rB);
         }
         ncdistMax = std::max( ncdistMax, dvncdist );
         distance[i] = dvncdist;

         previousV7Dist = v7Dist;
         sprintf( cBuffer, "Distance(v7,nc) %g  %g\n", v7Dist, dvncdist );
         buffer[lines] = cBuffer;
         ++lines;
         //printf( "V7Base %f %f %f %f %f %f %f\n", v7ReducedBase[0], v7ReducedBase[1], v7ReducedBase[2], v7ReducedBase[3], v7ReducedBase[4], v7ReducedBase[5], v7ReducedBase[6] );
         //printf( "V7Test %f %f %f %f %f %f %f\n", v7ReducedTest[0], v7ReducedTest[1], v7ReducedTest[2], v7ReducedTest[3], v7ReducedTest[4], v7ReducedTest[5], v7ReducedTest[6] );
      }

      bShouldPrint = true;

      for ( int i=0; bShouldPrint && i<lines; ++i )
      {
         printf( "%s", buffer[i].c_str() );
      }
   }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int main(int argc, char* argv[])
{


   Follower( );
}

