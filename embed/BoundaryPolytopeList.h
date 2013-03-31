#ifndef BOUNDARYPOLYTOPELIST_H
#define BOUNDARYPOLYTOPELIST_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>
#include <vector>

#include "BoundaryPolytope.h"

/*
class BoundaryPolytopeList
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Container for the primary boundaries of the Niggli cone in G6, enumerated by
   Andrews and Bernstein, 2012. There are 15 5-D boundary polytopes. They can be
   determined by Monte-Carlo methods or by enumeration by examining the Niggli
   reduction conditions. The boundary polytopes of lower dimension can be derived
   from the intersections of the primary boundaries. 

   BoundaryPolytopeList(void)                         == constructor that initializes the list of Boundary polytopes
   size_t size( void ) const                             == returns the number of items in the List
   const BoundaryPolytope& operator[] ( const int n ) == returns the n-th item in the list
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

class BoundaryPolytopeList
{
public:
   BoundaryPolytopeList( void );
   ~BoundaryPolytopeList( void );

   const BoundaryPolytope& operator[] ( const int n ) const;
 //      BoundaryPolytope& operator[] ( const int n ); // do NOT implement assignment
   int DegreesofFreedom( const int ITDESG );
   size_t size( void ) const { return( m_list.size() ); }
   std::vector<BoundaryPolytope> Distances( const double filter, const arma::vec6& v ) const;

private: // member data
   std::vector<BoundaryPolytope> m_list;
};

#endif // BOUNDARYPOLYTOPELIST_H
