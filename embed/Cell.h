#ifndef CELL_H
#define CELL_H

#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>
#include <cmath>

class Cell
{
public:
    Cell(void);
    Cell( const double a, const double b, const double c,
         const double alpha, const double beta, const double gamma);
    Cell( const arma::vec6& v6 );
    ~Cell(void);
    
    double Volume( void ) const;
    arma::vec6 CellWithDegrees( void ) const;
    double operator[](const int& n) const;
    arma::vec6 Cell2V6( void ) const;
    
    Cell Inverse( void ) const;
    arma::mat66 LatSymMat66(const std::string& latsym ) const;
    arma::mat66 mat33tomat66(const arma::mat33& m3);
    arma::mat33 mat66tomat33(const arma::mat66& m6);
    
private:
    arma::vec6 m_cell;
    double operator[](const int& n);
    void vvtorow(const double m,
                    const int n1, const int n2,
                    arma::mat33& v1, arma::mat33& v2,
                    const int n3,
                              arma::mat66& m6);
    void btos(int bits, arma::ivec& s);
};

#endif // CELL_H
