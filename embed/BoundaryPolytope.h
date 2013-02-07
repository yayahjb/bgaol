#ifndef BOUNDARYPOLYTOPE_H
#define BOUNDARYPOLYTOPE_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>
#include <string>

/*
class BoundaryPolytope
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Representation of the primary boundaries of the Niggli cone in G6, enumerated by
   Andrews and Bernstein, 2012. There are 15 5-D boundary polytopes. They can be
   determined by Monte-Carlo methods or by enumeration by examining the Niggli
   reduction conditions. The boundary polytopes of lower dimension can be derived
   from the intersections of the primary boundaries. 

   BoundaryPolytope( void ) == default constructor
   BoundaryPolytope::BoundaryPolytope( const std::string& proj, const std::string& prjhat, const std::string& transform,
                     const int degreesOfFreedom, const int ABindex, const std::string& ABname, const std::string& subspace,
                     const std::string& descr )
                     == constructor to entirely build a polytope
   DegreesOfFreedom( void ) == returns the number of degrees of freedom for a polytope (CAUTION: this is copied from the
                               fortran code and may be incorrect)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

class BoundaryPolytope
{

public:
   BoundaryPolytope( void );

   BoundaryPolytope( const std::string& proj, const std::string& prjhat, const std::string& tranform,
                     const int degreesOfFreedom, const int ABindex, const std::string& ABname, const std::string& subspace,
                     const std::string& descr );

   ~BoundaryPolytope( void );

   arma::mat66 GetProjector( void ) const { return( m_projector ); }
   arma::mat66 GetPerp     ( void ) const { return( m_perp ); }
   arma::mat66 GetTransform( void ) const { return( m_transform ); }

   int GetDOF                ( void ) const { return( m_degreesOfFreedom ); }
   int GetABindex            ( void ) const { return( m_ABindex ); }
   std::string GetABname     ( void ) const { return( m_ABname ); }
   std::string GetSubspace   ( void ) const { return( m_subspace ); }
   std::string GetDescription( void ) const { return( m_description ); }
   int DegreesofFreedom( const int ITDESG );

private: // member data
   arma::mat66 m_projector;   //    projector onto the boundary
   arma::mat66 m_perp;        //    projector onto the space orthogonal to the boundary (used of distance calculations)
   arma::mat66 m_prjhat;
   arma::mat66 m_transform;   //    the transformations applied to points as they pass thru the boundary
   int m_degreesOfFreedom;    // the number of independent degrees of freedom in the bounary

   // descriptions of the boundary
   int m_ABindex;             // the index number assigned by Andrews and Bernstein, 2012
   std::string m_ABname;      // the index as text
   std::string m_subspace;    // the symbolic representation of the subspace constraints
   std::string m_description; // text of the constraint equations for the boundary
};

#endif // BOUNDARYPOLYTOPE_H
