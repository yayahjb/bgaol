#ifndef MAKEGAOL_H
#define MAKEGAOL_H

#define USE_ARMADILLO_LIBRARY
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>

#include "BoundaryPolytopeList.h"
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

class MakeGaol
{
public:
   MakeGaol(void);
   ~MakeGaol(void);

   void MakeGaolEntry( arma::vec6 v[], const arma::vec6& gred, const arma::vec6& g6ErrorBox, const double ratio,
      int& nv, double vdist[], int ivb[] );

   double NCDist( const arma::vec6& gv1, const arma::vec6& gv2 );
   static void PrintG6 ( const std::string& text, const arma::vec6& v );
   static void PrintM66_INT ( const std::string& text, const arma::mat66& v );
   static bool inncone( const arma::vec6& gvec );
   void BDCoord( const arma::vec6& gvec, arma::vec6& xs, arma::vec6& ys );

private: // member data
   BoundaryPolytopeList m_prjList;

private: // functions
   bool GoodSeed( const int i );
   void GenSeeds( const arma::vec6& gred, const arma::vec6& g6ErrorBox, arma::vec6 seeds[], double seeddist[], bool goodseed[] ) const;
   static arma::mat66 MKPerp ( const arma::mat66& prj );
   double DistF( const arma::vec6& g ) const;
   double FoldMDist( const arma::vec6& gvec1, const arma::vec6& gvec2, const int ip, const double cFoldDist ) const;
   double FoldxDist( const arma::vec6& gvec1, const arma::vec6& gvec2, const int ip, const double cFoldDist,
      const int is1DUMMY, const int is2DUMMY) const;
   static void Set1000( const arma::vec6& fcmp1, const BoundaryPolytope& bt, const double ht1,
      arma::vec6& pfmcase1, double& hpfmc1, bool& lpfmc1 );
};

#endif //MAKEGAOL_H