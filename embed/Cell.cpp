
#include <cmath>

#include "Cell.h"

/*  class Cell
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A class to implement some common operations for unit cells as usually 
   used with xray crystallography. Conversions from G6 to unit cell and 
   from unit to G6 are included. The angles are ALWAYS in RADIANS.

   Cell(void)                                  == default constructor
   Cell( const arma::vec6& v )                 == constructor to convert a G6 vector to unit cell
   double Volume(void)                         == return the volume of a unit cell
   arma::vec6 Cell::CellWithDegrees            == return the unit cell as a vector with the angles as DEGREES
   double Cell::operator[](const int& n) const == return the n-th element of a cell (zero-based)
   Cell Cell::Inverse( void ) const            == compute the reciprocal cell
   arma::vec6 Cell::Cell2V6( void ) const      == return the G6 vector corresponding to a unit cell
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/


const double PI( 4.0 * atan(1.0) );
const double RAD2ANGLE( 180.0/PI );

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: Cell()
// Description: default constructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
Cell::Cell(void)
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: Cell()
// Description: constructor to convert an input G6 vector (as a vector of
//              doubles to E3 lengths and angles. The angles are 
//              stored as RADIANS (only). For angles as degrees, use the
//              function CellWithDegrees to obtain a VECTOR (arma::vec6)
//              with lengths and angles as DEGREES. For consistency, a
//              "Cell" object will never contain degrees.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
Cell::Cell( const arma::vec6& v )
{
   m_cell[0] = sqrt( v[0] );
   m_cell[1] = sqrt( v[1] );
   m_cell[2] = sqrt( v[2] );

   const double cosalpha( 0.5*v[3]/(m_cell[1]*m_cell[2]) );
   const double cosbeta ( 0.5*v[4]/(m_cell[0]*m_cell[2]) );
   const double cosgamma( 0.5*v[5]/(m_cell[0]*m_cell[1]) );

// compute the cell angles in radians
   m_cell[3] = atan2( sqrt(1.0-pow(cosalpha,2)),cosalpha);
   m_cell[4] = atan2( sqrt(1.0-pow(cosbeta ,2)),cosbeta );
   m_cell[5] = atan2( sqrt(1.0-pow(cosgamma,2)),cosgamma);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: Cell()
// Description: destructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
Cell::~Cell(void)
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: Volume()
// Description: Return the E3 volume of a Cell
// follows the formula of Stout and Jensen page 33
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double Cell::Volume( void ) const
{
   const arma::vec6& c( m_cell );
   const double c3(cos(c[3]));
   const double c4(cos(c[4]));
   const double c5(cos(c[5]));
   const double volume( c[0]*c[1]*c[2] * sqrt( 1.0-pow(c3,2)-pow(c4,2)-pow(c5,2)
       + 2.0*c3*c4*c5 ) );
   return( volume );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: CellWithDegrees()
// Description: Return the E3 cell with the angles as DEGREES in an
//              arma::vec6 vector of doubles.
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
arma::vec6 Cell::CellWithDegrees( void ) const
{
   arma::vec6 v;
   v[0] = m_cell[0];
   v[1] = m_cell[1];
   v[2] = m_cell[2];
   v[3] = RAD2ANGLE * m_cell[3];
   v[4] = RAD2ANGLE * m_cell[4];
   v[5] = RAD2ANGLE * m_cell[5];

   return( v );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: operator[]()
// Description: access function for the values in a Cell object
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double Cell::operator[](const int& n) const
{
   const int nn( std::max(0,std::min(5,n)) );
   return( m_cell[n] );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: Inverse()
// Description: Compute the reciprocal cell (inverse) of a cell
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
Cell Cell::Inverse( void ) const
{
   const double& a((*this)[0]);
   const double& b((*this)[1]);
   const double& c((*this)[2]);
   const double& alpha((*this)[3]);
   const double& beta ((*this)[4]);
   const double& gamma((*this)[5]);
   const double cosAlpha(cos(alpha));
   const double cosBeta (cos(beta));
   const double cosGamma(cos(gamma));
   const double sinAlpha(sin(alpha));
   const double sinBeta (sin(beta));
   const double sinGamma(sin(gamma));

   const double v     = (*this).Volume( );
   const double astar = b*c*sin(alpha)/v;
   const double bstar = a*c*sin(beta )/v;
   const double cstar = a*b*sin(gamma)/v;

   const double cosAlphaStar = (cosBeta *cosGamma-cosAlpha)/abs(sinBeta*sinGamma);
   const double cosBetaStar  = (cosAlpha*cosGamma-cosBeta )/abs(sinAlpha*sinGamma);
   const double cosGammaStar = (cosAlpha*cosBeta -cosGamma)/abs(sinAlpha*sinBeta);

   Cell cell;
   cell.m_cell[0] = astar;
   cell.m_cell[1] = bstar;
   cell.m_cell[2] = cstar;
   cell.m_cell[3] = atan2( sqrt(1.0-pow(cosAlphaStar,2)), cosAlphaStar);
   cell.m_cell[4] = atan2( sqrt(1.0-pow(cosBetaStar ,2)), cosBetaStar );
   cell.m_cell[5] = atan2( sqrt(1.0-pow(cosGammaStar,2)), cosGammaStar);

   return( cell );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: Cell2V6()
// Description: return the G6 vector of an input cell
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
arma::vec6 Cell::Cell2V6( void ) const
{
   const Cell& c( *this );
   arma::vec6 v;
   v[0] = c[0]*c[0];
   v[1] = c[1]*c[1];
   v[2] = c[2]*c[2];
   v[3] = 2.0*c[1]*c[2]*cos(c[3]);
   v[4] = 2.0*c[0]*c[2]*cos(c[4]);
   v[5] = 2.0*c[0]*c[1]*cos(c[5]);

   return( v );
}
