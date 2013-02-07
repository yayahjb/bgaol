#include "BasicDistance.h"

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: g123dist()
// Description: Compute the best distance between 2 G6 vectors
//              allowing for permutations of g1, g2, g3 as
//              well as sign changes
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double BasicDistance::g123dist( const arma::vec6& v1, const arma::vec6& v2 )
{
   arma::vec6 vtemp(v2); // copy v2 into vtemp
   //C     123
   double g123temp = g456dist( v1, vtemp );
   //C     213
   vtemp[0]=v2[1];
   vtemp[1]=v2[0];
   vtemp[3]=v2[4];
   vtemp[4]=v2[3];
   g123temp = std::min(g123temp,g456dist(v1,vtemp));
   //C     231
   vtemp[1]=v2[2];
   vtemp[2]=v2[0];
   vtemp[4]=v2[5];
   vtemp[5]=v2[3];
   g123temp = std::min(g123temp,g456dist(v1,vtemp));
   //C     321
   vtemp[0]=v2[2];
   vtemp[1]=v2[1];
   vtemp[3]=v2[5];
   vtemp[4]=v2[4];
   g123temp = std::min(g123temp,g456dist(v1,vtemp));
   //C     312
   vtemp[2]=v2[1];
   vtemp[1]=v2[0];
   vtemp[5]=v2[4];
   vtemp[4]=v2[3];
   g123temp = std::min(g123temp,g456dist(v1,vtemp));
   //C     132
   vtemp[0]=v2[0];
   vtemp[1]=v2[2];
   vtemp[3]=v2[3];
   vtemp[4]=v2[5];
   g123temp = std::min(g123temp,g456dist(v1,vtemp));
   return( g123temp );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: g123dist()
// Description: Compute the best distance between 2 G6 vectors
//              allowing for cell-preserving sign changes in
//              g4,5,6
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double BasicDistance::g456dist( const arma::vec6& v1, const arma::vec6& v2 )
{
   const double xdot = arma::dot( v1-v2, v1-v2 ); 
   const double g456temp = sqrt( xdot 
      + 4.0*std::min( std::min(0.0, 
               v1[3]*v2[3]+v1[4]*v2[4]),
      std::min(v1[3]*v2[3]+v1[5]*v2[5],
               v1[4]*v2[4]+v1[5]*v2[5]))) ;
   return( g456temp );
}
