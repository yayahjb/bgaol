#include "ProjectorTools.h"


ProjectorTools::ProjectorTools(void)
{
   InitCleanupList( );
}


ProjectorTools::~ProjectorTools(void)
{
   m_vCleanupList.clear();
}

//-----------------------------------------------------------------------------
// Name: InitCleanupList()
// Description: make a list of a lot of rational numbers this will be used 
//              to "clean up" the iterated projectors to correct for 
//              floating point approximations (if any)
//
//              NOTE: if you expand the size of this list, be careful, because
//              you may need to decrease the hard-coded constant 1.0E-4 in
//              function CleanupProjector
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int ProjectorTools::InitCleanupList( void )
{
   const long MAXDENOMINATOR = 30L;
   m_vCleanupList.push_back( 0.0 );

   for ( int inumerator=1; inumerator<=MAXDENOMINATOR-1; ++ inumerator )
   {
      for ( int idenominator=inumerator+1; idenominator<=MAXDENOMINATOR; ++idenominator )
      {
         m_vCleanupList.push_back( double(inumerator) / double(idenominator) );
         m_vCleanupList.push_back( -double(inumerator) / double(idenominator) );
      }
   }

   return( m_vCleanupList.size( ) );
}

//-----------------------------------------------------------------------------
// Name: Squaring()
// Description: Performs n squarings of an input projector to find the limiting
//              projector (Andrews and Bernstein, 1986)
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
arma::mat66 ProjectorTools::Squaring( const int n, const arma::mat66& prjin )
{
   arma::mat66 m( prjin );
   for ( int i=0; i<n; ++i )
   {
      m *=m;
   }

   return( m );
}

//-----------------------------------------------------------------------------
// Name: CleanupProjector()
// Description: make a list of a lot of rational numbers this will be used 
//              to "clean up" the iterated projectors to correct for 
//              floating point approximations (if any)
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void ProjectorTools::CleanupProjector( arma::mat66& projector )
{
   const size_t n = m_vCleanupList.size( );
   for ( int iproj=0; iproj<36; ++iproj )
   {
      for ( int i=0; i<n; ++i )
      {
         if ( fabs( projector[iproj]-m_vCleanupList[i] ) < 1.0E-4 )
         {
            projector[iproj] = m_vCleanupList[i];
            break;
         }
      }
   }
}
