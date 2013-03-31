
#include <cmath>
#include <cstdio>

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
const int fi[37] = { 0, 0, 6,12,18,24,30,
                        1, 7,13,19,25,31,
                        2, 8,14,20,26,32,
                        3, 9,15,21,27,33,
                        4,10,16,22,28,34,
                        5,11,17,23,29,35};

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: Cell()
// Description: default constructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
Cell::Cell(void)
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: Cell()
// Description: constructor to convert load a Cell from doubles
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
Cell::Cell( const double a, const double b, const double c,
           const double alpha, const double beta, const double gamma)
{
    m_cell[0] = a;
    m_cell[1] = b;
    m_cell[2] = c;
    m_cell[3] = alpha/RAD2ANGLE;
    m_cell[4] = beta/RAD2ANGLE;
    m_cell[5] = gamma/RAD2ANGLE;
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

   const double cosAlphaStar = (cosBeta *cosGamma-cosAlpha)/fabs(sinBeta*sinGamma);
   const double cosBetaStar  = (cosAlpha*cosGamma-cosBeta )/fabs(sinAlpha*sinGamma);
   const double cosGammaStar = (cosAlpha*cosBeta -cosGamma)/fabs(sinAlpha*sinBeta);

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


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: LatSymMat66(const std::string& latsym)
// Description: return the mat66 matrix for a lattice symbol
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
arma::mat66 Cell::LatSymMat66( const std::string& latsym ) const
{
    arma::mat66 M;
    arma::vec6 g;
    double ge;
    char LATSYM;
    int i;
    
    M.zeros();
    if (latsym.size()!=1) return(M);
    LATSYM = latsym[0];
    
    switch(LATSYM){
        case 'P':
        case 'p':
        case 'V':
        case 'v':
            M[fi[1]]=M[fi[8]]=M[fi[15]]=M[fi[22]]=M[fi[29]]=M[fi[36]]=1.;
            break;
            
        case 'I':
        case 'i':
            M[fi[1]] = 1.;
            M[fi[8]] = 1.;
            for (i=13;i<=18;i++) {
                M[fi[i]] = 0.25;
            }
            M[fi[20]] = 1.;
            M[fi[22]] = 0.5;
            M[fi[24]] = 0.5;
            M[fi[25]] = 1.;
            M[fi[29]] = 0.5;
            M[fi[30]] = 0.5;
            M[fi[36]] = 1.;
            break;

        case 'F':
        case 'f':
            M[fi[1]] = 0.25;
            M[fi[2]] = 0.25;
            M[fi[6]] = 0.25;
            M[fi[7]] = 0.25;
            M[fi[9]] = 0.25;
            M[fi[11]] = 0.25;
            M[fi[14]] = 0.25;
            M[fi[15]] = 0.25;
            M[fi[16]] = 0.25;
            M[fi[21]] = 0.5;
            for (i=22;i<=24;i++) {
                M[fi[i]] = 0.25;
            }
            M[fi[26]] = 0.5;
            for (i=28;i<=30;i++) {
                M[fi[i]] = 0.25;
            }
            M[fi[31]] = 0.5;
            for (i=34;i<=36;i++) {
                M[fi[i]] = 0.25;
            }
            break;
            
        case 'A':
        case 'a':
            M[fi[1]] = 1.;
            M[fi[8]] = 1.;
            M[fi[14]] = 0.25;
            M[fi[15]] = 0.25;
            M[fi[16]] = 0.25;
            M[fi[20]] = 1.;
            M[fi[22]] = 0.5;
            M[fi[29]] = 0.5;
            M[fi[30]] = 0.5;
            M[fi[36]] = 1.;
            break;

            
        case 'B':
        case 'b':
            M[fi[1]] = 1.;
            M[fi[8]] = 1.;
            M[fi[13]] = 0.25;
            M[fi[15]] = 0.25;
            M[fi[17]] = 0.25;
            M[fi[22]] = 0.5;
            M[fi[24]] = 0.5;
            M[fi[25]] = 1.;
            M[fi[29]] = 0.5;
            M[fi[36]] = 1.;
            break;

            
        case 'C':
        case 'c':
            M[fi[1]] = 1.;
            M[fi[7]] = 0.25;
            M[fi[8]] = 0.25;
            M[fi[12]] = 0.25;
            M[fi[15]] = 1.;
            M[fi[22]] = 0.5;
            M[fi[23]] = 0.5;
            M[fi[29]] = 1.;
            M[fi[31]] = 1.;
            M[fi[36]] = 0.5;
            break;

        case 'R':
        case 'r':
        {   arma::vec6 cwd=(*this).CellWithDegrees();
            g = (*this).Cell2V6();
            ge = arma::norm(g,2)*.005;
            // For R, as distinct from H, detect primitive R cases
            if (std::abs(g[0]-g[1])<ge) {
                // (r,r,?,?,?,?)
                if (std::abs(g[3]-g[4]) < ge
                   &&  std::abs(g[4]-g[5]) < ge) {
                // (r,r,?,s,s,s), (r,r,?,-s,-s,-s)
                    if (std::abs(g[3]-g[0]) < ge
                        || std::abs(g[1]-g[2]) < ge ) {
                        arma::vec6 cwd=(*this).CellWithDegrees();
                        // (r,r,s,r,r,r), (r,r,r,s,s,s), (r,r,r,-s,-s,-s)
                        M[fi[1]]=M[fi[8]]=M[fi[15]]=M[fi[22]]=M[fi[29]]=M[fi[36]]=1.;
                        //fprintf(stderr,"Treated non-hexagonal R as P: %g %g %g %g %g %g\n",
                        //        cwd[0],cwd[1],cwd[2],cwd[3],cwd[4],cwd[5]);
                        break;
                    }
                }
            }
            if (std::abs(g[1]-g[2]) < ge ) {
                if (std::abs(3.*g[4]/2.+g[0]) < ge
                    && std::abs(3.*g[5]/2.+g[0]) < ge
                    && std::abs(3.*(g[3]+g[1])-g[0]) < ge){
                    arma::vec6 cwd=(*this).CellWithDegrees();
                    M[fi[1]]=M[fi[8]]=M[fi[15]]=M[fi[22]]=M[fi[29]]=M[fi[36]]=1.;
                    // fprintf(stderr,"Treated non-hexagonal R as P: %g %g %g %g %g %g\n",
                    //        cwd[0],cwd[1],cwd[2],cwd[3],cwd[4],cwd[5]);
                    break;
                }
            }
            fprintf(stderr,"Treated non-rhombodedral R as H: %g %g %g %g %g %g\n",
                   cwd[0],cwd[1],cwd[2],cwd[3],cwd[4],cwd[5]);
        }

        case 'H':
        case 'h':
            // from Andrews and Bernstein, 1988
            M[fi[1]] = 1./9.;
            M[fi[2]] = 1./9.;
            M[fi[3]] = 1./9.;
            M[fi[4]] = 1./9.;
            M[fi[5]] = -1./9.;
            M[fi[6]] = -1./9.;

            M[fi[7]] = 4./9.;
            M[fi[8]] = 1./9.;
            M[fi[9]] = 1./9.;
            M[fi[10]] = 1./9.;
            M[fi[11]] = 2./9.;
            M[fi[12]] = 2./9.;

            M[fi[13]] = 1./9.;
            M[fi[14]] = 4./9.;
            M[fi[15]] = 1./9.;
            M[fi[16]] = -2./9.;
            M[fi[17]] = -1./9.;
            M[fi[18]] = 2./9.;

            M[fi[19]] = -4./9.;
            M[fi[20]] = -4./9.;
            M[fi[21]] = 2./9.;
            M[fi[22]] = -1./9.;
            M[fi[23]] = 1./9.;
            M[fi[24]] = -5./9.;

            M[fi[25]] = 2./9.;
            M[fi[26]] = -4./9.;
            M[fi[27]] = 2./9.;
            M[fi[28]] = -1./9.;
            M[fi[29]] = -2./9.;
            M[fi[30]] = 1./9.;

            M[fi[31]] = -4./9.;
            M[fi[32]] = 2./9.;
            M[fi[33]] = 2./9.;
            M[fi[34]] = 2./9.;
            M[fi[35]] = 1./9.;
            M[fi[36]] = 1./9.;
            break;
    }
   
    return( M );
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: vvtorow(const double m,
//               const int n1, const int n2,
//               arma::mat33& v1, arma::mat33& v2,
//               const int n3,
//               arma::mat66& m6)
// Description: Compute one row of a 6x6 matrix from
//              2 rows of 2 3x3 matrices
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void vvtorow(const double m,
                    const int n1, const int n2,
                    const arma::mat33& v1, const arma::mat33& v2,
                    const int n3,
                    arma::mat66& m6) {
    m6(n3,0) = m*v1(n1,0)*v2(n2,0);
    m6(n3,1) = m*v1(n1,1)*v2(n2,1);
    m6(n3,2) = m*v1(n1,2)*v2(n2,2);
    m6(n3,3) = (m*v1(n1,1)*v2(n2,2) + m*v1(n1,2)*v2(n2,1)) / 2.0;
    m6(n3,4) = (m*v1(n1,0)*v2(n2,2) + m*v1(n1,2)*v2(n2,0)) / 2.0;
    m6(n3,5) = (m*v1(n1,0)*v2(n2,1) + m*v1(n1,1)*v2(n2,0)) / 2.0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: mat33tomat66(arma::mat33& m3)
// Description: Compute a 6x6 matrix from a 3x3 matrix
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

arma::mat66 mat33tomat66(const arma::mat33& m3) {
    arma::mat66 m6;
    vvtorow(1.0,0,0,m3,m3,0,m6);
    vvtorow(1.0,1,1,m3,m3,1,m6);
    vvtorow(1.0,2,2,m3,m3,2,m6);
    vvtorow(2.0,1,2,m3,m3,3,m6);
    vvtorow(2.0,0,2,m3,m3,4,m6);
    vvtorow(2.0,0,1,m3,m3,5,m6);
    return m6;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: btos(int bits, arma::ivec s)
// Description: return vector of signs of bits
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void btos(int bits, arma::ivec& s) {
    int i;
    for (i=0; i < s.n_rows; i++) {
        if ((bits >> i)&1 == 1) s[i] = -1;
        else s[i] = 1;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: mat66tomat33(arma::mat66& m6)
// Description: return a 3x3 matrix for a 6x6 matrix
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
arma::mat33 mat66tomat33(const arma::mat66& m6) {
    arma::mat33 m;
    arma::mat33 m3;
    arma::ivec9 signs;
    arma::mat66 m6t;
    double t1, t1min;
    int bits, bitsmin;
    int i, j, l;
    
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            m3(i,j) = std::sqrt(m6(i,j));
        }
    }
    
    for (bits = 511; bits >=0; bits--) {
        btos(bits,signs);
        for (l = 0; l < 9; l++) {
            i = l/3;
            j = l - 3*i;
            m(i,j) = signs(l)*m3(i,j);
        }
        m6t = mat33tomat66(m);
        t1 = arma::norm((m6-m6t),2);
        if (bits == 511) t1min=t1;
        if (t1 <= t1min) {
            t1min = t1;
            bitsmin = bits;
        }
        if (t1 == 0.0) break;
    }
    if (t1 != 0.0) {
        btos(bitsmin,signs);
        for (l=0; l<9; l++) {
            i = l/3;
            j = l - 3*i;
            m(i,j) = signs(l)*m(i,j);
        }
    }
    if (arma::det(m,true)) {
        for (l=0; l<9; l++) {
            i = l/3;
            j = l - 3*i;
            m(i,j) = -m(i,j);
        }        
    }
    for (l=0; l<9; l++) {
        i = l/3;
        j = l - 3*i;
        if (std::abs(m(i,j))< 1.e-8) m(i,j) = 0.;
    }
    return m;
}
