#ifndef REDUCER_H
#define REDUCER_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>

/*  class Reducer
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A class to implement Niggli reduction in the G6 space. Reduce returns a reduced cell
   and the matrix that converts the input cell to the reduced cell.

   Reduced(void) 
                 == constructor -- nothing to do here
   static void Reduce( const arma::vec6& vi, arma::mat66& m, arma::vec6& vout, const double delta )
                 == returns vout as the reduced cell of vi and m, the conversion matrix
   static bool NearRed( const arma::vec6& gvec, const double delta );
                 == determines whether gvec is reduced within delta
   static bool Near2Red( const arma::vec6& gvec, const double delta, arma::vec6& vout, double& dist )
                 == determines whether gvec is reduced within delta, and returns the reduced cell
                    and how far from reduced
   static void Reporter( const std::string& text, const arma::vec6& vin, const arma::vec6& vout, const arma::mat66& m )
                 == prints information about each step in reduction (including standard presentation)

private:
   static void MKnorm( const arma::vec6& vi, arma::mat66& m, arma::vec6& vout, const double delta  )
                 == internal function to convert vi to standard presentation - the matrix
                    to implement that change is returned
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

class Reducer
{
public:
   static bool Reduce( const arma::vec6& vi, arma::mat66& m, arma::vec6& vout, const double delta );
   static bool NearRed( const arma::vec6& gvec, const double delta );
   static bool Near2Red( const arma::vec6& gvec, const double delta, arma::vec6& vout, double& dist );

private:
   static void MKnorm( const arma::vec6& vi, arma::mat66& m, arma::vec6& vout, const double delta  );
   static void Reporter( const std::string& text, const arma::vec6& vin, const arma::vec6& vout, const arma::mat66& m );

   // at least for now, all functions are static, and there is no member data
   // forbid constructor and destructor
   Reducer(void);
   ~Reducer(void);

public:

};

#endif // REDUCER_H

