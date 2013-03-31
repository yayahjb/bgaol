#include "BoundaryPolytopeTracker.h"


BoundaryPolytopeTracker::BoundaryPolytopeTracker(void)
{
}

BoundaryPolytopeTracker::~BoundaryPolytopeTracker(void)
{
}


BoundaryPolytopeTracker BoundaryPolytopeTracker::operator<< ( BoundaryPolytope const& bp )
{
   m_list.push_back( bp );
   return bp;
}


BoundaryPolytope BoundaryPolytopeTracker::operator[] ( const int n) const
{
   return( m_list[n] );
}

size_t BoundaryPolytopeTracker::size( void ) const
{
   return( m_list.size() );
}
