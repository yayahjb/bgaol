#ifndef BASICDISTANCE_H
#define BASICDISTANCE_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>
class BasicDistance

/*
Functions for determining partial distances in G6. They are placed
here in a support class because they are used both in Reducer
and in MakeGaol.
*/
{
public:
   static double g456dist( const arma::vec6& v1, const arma::vec6& v2 );
   static double g123dist( const arma::vec6& v1, const arma::vec6& v2 );
};

#endif // BASICDISTANCE_H
