#ifndef CELL_H
#define CELL_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>

class Cell
{
public:
   Cell(void);
   Cell( const arma::vec6& v6 );
   ~Cell(void);

   double Volume( void ) const;
   arma::vec6 CellWithDegrees( void ) const;
   double operator[](const int& n) const;
   arma::vec6 Cell2V6( void ) const;

   Cell Inverse( void ) const;

private:
   arma::vec6 m_cell;
   double operator[](const int& n);
};

#endif // CELL_H
