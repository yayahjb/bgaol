#ifndef NIGGLIPOLYTOPE_H
#define NIGGLIPOLYTOPE_H

#include <string>

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>


class NiggliPolytope
{
public:
   NiggliPolytope( void );
   NiggliPolytope( const std::string& lattice, const int IT_number, const std::string& niggliRoof, const std::string& subspace,
                   const std::string& descr, const double projectorNormalizer, const std::string& projector );
   ~NiggliPolytope( void );

   arma::mat66 GetProjector   ( void ) const { return( m_projector ); }
   arma::mat66 GetPerp        ( void ) const { return( m_perp ); }
   std::string GetLattice     ( void ) const { return( m_lattice ); }
   int         GetITnumber    ( void ) const { return( m_IT_number ); }
   std::string GetNiggliRoofID( void ) const { return( m_NiggliRoofID ); }
   std::string GetSubspace    ( void ) const { return( m_subspace ); }
   std::string GetDescription ( void ) const { return( m_description ); }

private:
   double m_projectorNormalizer;
   arma::mat66 m_projector;
   arma::mat66 m_perp;

   std::string m_lattice; // such as "cI"
   int         m_IT_number;
   std::string m_NiggliRoofID;
   std::string m_subspace;
   std::string m_description; // free form description
};

#endif //NIGGLIPOLYTOPE_H