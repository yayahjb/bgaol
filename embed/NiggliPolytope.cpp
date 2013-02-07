#include "NiggliPolytope.h"

/*  class NiggliPolytope
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A class initialize and contain the G6 representations of the 44 Niggli forms and
   their various properties.

   NiggliPolytope( void ) == default constructor

   NiggliPolytope( const std::string& lattice, const int IT_number, const std::string& niggliRoof, const std::string& subspace,
                   const std::string& descr, const double projectorNormalizer, const std::string& projector );
                 == constructor to set the properties of one Niggli form

   arma::mat66 GetProjector   ( void ) const == return the projector that takes a G6 vector and projects onto the polytope
   arma::mat66 GetPerp        ( void ) const == return the "perp" that takes a G6 vector and projects onto the normal to the polytope
   std::string GetLattice     ( void ) const == return the lattice type description (such as cF for face-centered cubic)
   int         GetITnumber    ( void ) const == return the ordinal assigned in the International Tables for Crystallography for the polytope
   std::string GetNiggliRoofID( void ) const == return the label assigned by Niggli, 1928, and Roof, 1969
   std::string GetSubspace    ( void ) const == return the subspace description of Andrews and Bernstein, 1988 (such as r,r,r,s,s,s)
   std::string GetDescription ( void ) const == some descriptive text
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/


NiggliPolytope::NiggliPolytope( void )
{
}

NiggliPolytope::NiggliPolytope( const std::string& lattice, const int IT_number, const std::string& niggliRoof, const std::string& subspace,
                                const std::string& descr, const double projectorNormalizer, const std::string& projector)
   : m_projectorNormalizer( projectorNormalizer )
   , m_projector          ( projector )
   , m_perp               ( arma::eye(6,6) - m_projector )
   , m_lattice            ( lattice ) // such as "cI"
   , m_IT_number          ( IT_number )
   , m_NiggliRoofID       ( niggliRoof )
   , m_subspace           ( subspace )
   , m_description        ( descr ) // free form description
{
   m_projector /= m_projectorNormalizer;
}

NiggliPolytope::~NiggliPolytope(  )
{
}
