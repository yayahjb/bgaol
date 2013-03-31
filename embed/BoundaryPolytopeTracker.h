#ifndef BOUNDARYPOLYTOPETRACKER_H
#define BOUNDARYPOLYTOPETRACKER_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>
#include <string>
#include <vector>

#include "BoundaryPolytope.h"

class BoundaryPolytopeTracker : public BoundaryPolytope
{
public:
   BoundaryPolytopeTracker( void );
   BoundaryPolytopeTracker( const BoundaryPolytope& bp );
   ~BoundaryPolytopeTracker( void );

   BoundaryPolytopeTracker operator<< ( const BoundaryPolytope& bp );
   BoundaryPolytope operator[] ( const int n) const;
   size_t size( void ) const;
   double FirstDistance( void ) const { return( m_firstDistance ); }
   void FirstDistance( const double d ) { m_firstDistance = d; };
   double TotalDistance( void ) const { return( m_totalDistance ); }
   void TotalDistance( const double d ) { m_totalDistance = d; };
   arma::vec6 ProjectedVector( void ) const { return( m_projectedVector ); }
   void ProjectedVector( const arma::vec6& v ) { m_projectedVector = v; }

private:
   std::vector<BoundaryPolytope> m_list;
   double m_firstDistance;
   double m_totalDistance;
   arma::vec6 m_seedVector;
   arma::vec6 m_projectedVector;
};

#endif //BOUNDARYPOLYTOPETRACKER_H