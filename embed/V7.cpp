#include "V7.h"

#include "Reducer.h"

#include <armadillo>

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name:V7 ()
// Description: default constructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
V7::V7(void)
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: V7()
// Description: copy constructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
V7::V7( const arma::vec7& v7 )
   :m_v7( v7 )
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: V7()
// Description: constructor from a G6 vector
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
V7::V7( const arma::vec6& v )
{
    m_v7[0] = std::sqrt(v[0]);
    m_v7[1] = std::sqrt(v[1]);
    m_v7[2] = std::sqrt(v[2]);

   const Cell c( v );
   const Cell cI( c.Inverse( ) );

   const arma::vec6 vinverse( cI.Cell2V6( ) );

   arma::mat66 m;
   arma::vec6 reducedBase;
   Reducer::Reduce( vinverse, m, reducedBase, 0.0 );


    m_v7[3] = std::sqrt(1.0/reducedBase[0]);
    m_v7[4] = std::sqrt(1.0/reducedBase[1]);
    m_v7[5] = std::sqrt(1.0/reducedBase[2]);

   m_v7[6] = pow( c.Volume( ), 1.0/3.0 );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: V7()
// Description: constructor from a Cell
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
V7::V7( const Cell& c )
{
   const arma::vec6 v( c.Cell2V6() );
   m_v7 = v;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: ~V7()
// Description: 
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
V7::~V7(void)
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: norm()
// Description: 
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double V7::Norm(void ) const
{
   return( arma::norm( m_v7, 2 ) );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: operator-()
// Description: 
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
const V7 V7::operator- ( const V7& v7 ) const
{
   const arma::vec7 vtemp( m_v7-v7.m_v7 );
   return( V7(vtemp) );
}
