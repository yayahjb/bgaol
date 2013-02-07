#ifndef V7_H
#define V7_H

#include <armadillo>
#include "Cell.h"

class V7
{
public:
   V7(void);
   explicit V7( const arma::vec7& v7 );
   explicit V7( const arma::vec6& v );
   explicit V7( const Cell& c );
   ~V7(void);

   double Norm( void ) const;
   const V7 operator- ( const V7& v7 ) const;
   double operator[](const int n ) const { return( m_v7[n] ); }

private:
   arma::vec7 m_v7;
};

#endif //V7_H
