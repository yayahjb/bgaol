#ifndef PROJECTORTOOL_H
#define PROJECTORTOOL_H

#define USE_ARMADILLO_LIBRARY
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>

#include <vector>

class ProjectorTools
{
public:
   ProjectorTools(void);
   ~ProjectorTools(void);

   int InitCleanupList( void );
   void CleanupProjector( arma::mat66& projector );
   static arma::mat66 Squaring( const int n, const arma::mat66& prjin );

private:
   std::vector<double> m_vCleanupList;

};

#endif // PROJECTORTOOL_H
//InitCleanupList
//   CleanupProjector
//   Squaring