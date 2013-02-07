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


#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>
#include <string>

#include "BoundaryPolytope.h"
#include "BoundaryPolytopeList.h"
#include "NiggliPolytope.h"
#include "NiggliPolytopeList.h"
#include "Cell.h"
#include "MakeGaol.h"
#include "Reducer.h"
#include "V7.h"
#include "Vector_3d.h"

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
      if ( abs(Cell(testBase).Volume( ) - Cell(reducedBase).Volume( )) > 0.01 )
      {
         const double v1 = Cell(testBase).Volume( );
         const double v2 = Cell(reducedBase).Volume( );
         const int i19191 = 19191;
      }
      const V7 v7ReducedBase( reducedBase );
      MakeGaol mkg;
         printf( "\n\n\nV7Base %f %f %f %f %f %f %f\n", v7ReducedBase[0], v7ReducedBase[1], v7ReducedBase[2], v7ReducedBase[3], v7ReducedBase[4], v7ReducedBase[5], v7ReducedBase[6] );

      //   exit(0);

      NiggliPolytope np;
      NiggliPolytopeList npl;
      double previousV7Dist  = 0.0;
      const int nstep(100);

      printf( "jfollow %d   ", jfollow );
      arma::mat66 mtemp;
      arma::vec6 vprojReduced;
      arma::vec6 reducedTest;

      // Loop to run from the randomized starting point to the reduced
      // version of the randomized starting point
      for ( int i=0; i<=nstep; ++i )
      {
//         if ( i <11 || i > 27 ) continue;
         const double di(i/double(nstep));
         // increment along the path
         const arma::vec6 test( di*reducedBase + (1.0-di)*testBase );
         // Reduce the current probe lattice
         Reducer::Reduce( test, m, reducedTest, 0.0 );

         // loop over the 44 Niggli lattice types
         for ( int iniggli=1; iniggli<npl.size(); ++ iniggli )
         {
            // project the reduced test lattice onto the current boundary
            const arma::vec6 vprojected( npl[iniggli].GetProjector() * reducedTest );
            // make sure the projection is reduced
            const bool bRed = Reducer::Reduce( vprojected, m, vprojReduced, 0.0 );
            if ( bRed )
            {
               // calculate how far the projected point is from the reduced probe point
               const double dist = mkg.NCDist( vprojReduced, reducedTest );
               const bool isInCone = MakeGaol::inncone( reducedTest );// to get the relative direction whether we are inside or outside of the cone
               if ( ! isInCone )
               {
                  // just a place to put a breakpoint for point that really wasn't reduced
                  const int i19191 = 19191;
               }
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
   arma::vec6 reducedBase;
   arma::mat66 m;
   m.eye( );
   // const arma::vec6 base("9.997321 10.003294 10.000253 0.001003 9.996918 9.997098");
   const arma::vec6 base("10 10 10   10 10 10");
   Reducer::Reduce( base, m, reducedBase, 0.0 );



   //     exit(0);

   const int nfollow=100;
   for ( int jfollow=0; jfollow<nfollow; ++jfollow )
   {
      const arma::vec6 rands = GenRandG6( );
//      if( jfollow != 13 ) continue;

      const arma::vec6 testBase( base + 0.01 * rands );
      Reducer::Reduce( testBase, m, reducedBase, 0.0 );

      if ( abs(Cell(testBase).Volume( ) - Cell(reducedBase).Volume( )) > 0.01 )
      {
         const double v1 = Cell(testBase).Volume( );
         const double v2 = Cell(reducedBase).Volume( );
         const int i19191 = 19191;
      }
      const V7 v7ReducedBase( reducedBase );
      MakeGaol mkg;
         printf( "\n\n\nV7Base %f %f %f %f %f %f %f\n", v7ReducedBase[0], v7ReducedBase[1], v7ReducedBase[2], v7ReducedBase[3], v7ReducedBase[4], v7ReducedBase[5], v7ReducedBase[6] );

      //   exit(0);

      double previousV7Dist  = 0.0;
      const int nstep(100);
      for ( int i=0; i<=nstep; ++i )
      {
//         if ( i <11 || i > 27 ) continue;
         const double di(i/double(nstep));
         const arma::vec6 test( di*reducedBase + (1.0-di)*testBase );
         arma::vec6 reducedTest;
         Reducer::Reduce( test, m, reducedTest, 0.0 );
         const V7         v7ReducedTest( reducedTest );
         const double     v7Dist       ( (v7ReducedTest-v7ReducedBase).Norm() );
         const double     dvncdist     ( mkg.NCDist(reducedTest, reducedBase) );
         const arma::vec6 v6diff       ( reducedTest - reducedBase );
         const double     v6diffNorm   ( arma::norm( v6diff,2 ) );

         if ( i              == 0 ) previousV7Dist = v7Dist;

         if ( abs(v7Dist-previousV7Dist) > 1.0e-3 )
         {
            const int i19191 = 19191; // place to put a breakpoint or assert to indicate large steps
         }
         previousV7Dist = v7Dist;
         printf( "\nDistance(v7,nc) %f  %f\n", v7Dist, dvncdist );
         //printf( "V7Base %f %f %f %f %f %f %f\n", v7ReducedBase[0], v7ReducedBase[1], v7ReducedBase[2], v7ReducedBase[3], v7ReducedBase[4], v7ReducedBase[5], v7ReducedBase[6] );
         //printf( "V7Test %f %f %f %f %f %f %f\n", v7ReducedTest[0], v7ReducedTest[1], v7ReducedTest[2], v7ReducedTest[3], v7ReducedTest[4], v7ReducedTest[5], v7ReducedTest[6] );
      }
      printf( "\njfollow= %d\n\n", jfollow );
   }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int main(int argc, char* argv[])
{

   //test1();

//   Follower( );
   NiggliClassIdentifier( );

//   BoundaryDistances();
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   //BoundaryPolytopeList btl;
   //BoundaryPolytope bp = btl[0x1];

   //CNearTree<V7> nt;
   //const arma::vec6 g1("10.1 10.2 10.3 10.4 10.5 10.6");
   //const arma::vec6 g2 = bp.GetProjector()*g1;
   //const V7 v1(g1);
   //const V7 v2(g2);

   //nt.insert( v1 );
   //nt.insert( v2 );

   //const CNearTree<V7>::iterator it = nt.NearestNeighbor( 10.0, v2 );
   //const double dvnt( (v1-v2).Norm( ) );
   //MakeGaol mkg;
   //const double dvncdist( mkg.NCDist(arma::vec6(g1), arma::vec6(g2)) );

   //const int i19191 = 19191;
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   //BoundaryPolytopeList btl;
   //BoundaryPolytope bp = btl[0x1];

   //const arma::vec6 g1("10.1 10.2 10.3 10.4 10.5 10.6");
   //const arma::vec6 g2 = bp.GetProjector()*g1;

   //arma::vec6 gout;
   //arma::mat66 m;
   //m.eye( );
   //Reducer::Reduce( g1, m, gout, 0.0 );
   //MakeGaol mkg;
   //const double ncdist = mkg.NCDist( g1, g2 );


   ////void MakeGaol::MakeGaolEntry( CNearTree<arma::vec6>& nt, arma::vec6 v[], arma::vec6& gred, 
   ////arma::vec6& ge, const double ratio, double vdist[], int ivb[] )


   //arma::vec6 returnList[1000];
   //double vdist[1000];
   //int ivb[1000];
   //arma::vec6 gred(g1), g6ErrorBox(g2);
   //int nv;
   //mkg.MakeGaolEntry( returnList, gred, g6ErrorBox, 1.0, nv, vdist, ivb );

   //const arma::vec6 gi(".1 .11 .12 .1 .1 .1");

   //Reducer::Reduce( gi, m, gout, 0.0 );

   //const V7 v7( g1 );


   return 0;
}

