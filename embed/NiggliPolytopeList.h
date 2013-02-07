#ifndef NIGGLIPOLYTOPELIST_H
#define NIGGLIPOLYTOPELIST_H

#include <vector>

#include "NiggliPolytope.h"


/* class NiggliPolytopeList
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Description: initializes the 44 Niggli lattice types. There are 45 entries
             so that the index of each is the same as the number assigned
             in the International Table of Crystallography. Also included
             are the lattice type (cI, oP, etc.) and the numbering from
             the figures in Niggli's original publication which was 
             updated/corrected by Roof (44C, etc.). The subspace description
             from Andrews and Bernstein, 1988 (r,r,r,s,s,s and so on) are input
             The projectors are from Paciorek and Bonin (1992). The projectors
             in not normalized until divided by the normalizer, so that all
             the entries are integer values.

   NiggliPolytopeList(void)                         == constructor that initializes the list of Niggli polytopes
   int size( void ) const                           == returns the number of items in the List
   const NiggliPolytope& operator[] ( const int n ) == returns the n-th item in the list
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

class NiggliPolytopeList
{
public:
   NiggliPolytopeList(void);
   ~NiggliPolytopeList(void);
   int size( void ) const { return( m_list.size( ) ); };
   const NiggliPolytope& operator[] ( const int n ) const;

private: // member data
   std::vector<NiggliPolytope> m_list;

};

#endif // NIGGLIPOLYTOPELIST_H
