
#include <string>

#include "BoundaryPolytope.h"


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


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: BoundaryPolytope()
// Description: simple default constructor for building arrays
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
BoundaryPolytope::BoundaryPolytope( void )
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: BoundaryPolytope()
// Description: constructor to initialize all data in individual polytopes 
//              By design there are no setting functions
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
BoundaryPolytope::BoundaryPolytope( const std::string& proj, const std::string& prjhat, const std::string& transform,
                  const int degreesOfFreedom, const int ABindex, const std::string& ABname, const std::string& subspace,
                  const std::string& descr )
                  : m_projector( proj )
                  , m_prjhat( prjhat )
                  , m_perp( arma::eye(6,6) - m_projector )
                  , m_transform( transform )
                  , m_degreesOfFreedom( degreesOfFreedom )
                  , m_ABindex( ABindex )
                  , m_ABname( ABname )
                  , m_subspace( subspace )
                  , m_description( descr )
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: ~BoundaryPolytope()
// Description: destructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
BoundaryPolytope::~BoundaryPolytope( void )
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: DegreesofFreedom()
// Description: I BELIEVE THIS IS NOT CORRECT !!!!!!!!!!!!!!!!!!!!!!!!!
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int BoundaryPolytope::DegreesofFreedom( const int ITDESG )
{
   static const
      int idf[] = { 0,1,0,1,0,1,1,2,1,3,
                    1,1,2,3,1,2,3,1,2,3,
                    1,1,2,1,3,2,3,3,3,3,
                    5,2,3,3,3,2,3,2,3,2,
                    3,2,3,5 };
   return( idf[std::max(1,std::min(abs(ITDESG),44))] );                            
}

