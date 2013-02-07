
#include "BoundaryPolytopeList.h"
/*
class BoundaryPolytopeList
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Container for the primary boundaries of the Niggli cone in G6, enumerated by
   Andrews and Bernstein, 2012. There are 15 5-D boundary polytopes. They can be
   determined by Monte-Carlo methods or by enumeration by examining the Niggli
   reduction conditions. The boundary polytopes of lower dimension can be derived
   from the intersections of the primary boundaries. 

   BoundaryPolytopeList(void)                         == constructor that initializes the list of Boundary polytopes
   int size( void ) const                             == returns the number of items in the List
   const BoundaryPolytope& operator[] ( const int n ) == returns the n-th item in the list
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: BoundaryPolytopeList()
// Description: initializes the 15 boundary manifolds and some extras.
//              The list begins with a null (meaningless) entry so that the 
//              index of each is the same as the number assigned
//              in Andrews and Bernstein, 2012. 
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
BoundaryPolytopeList::BoundaryPolytopeList( void )
   : m_list(26)
{
   m_list[0x0] = BoundaryPolytope(  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // projector
                                    "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // prjhat
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1", // transform
                                    0, 0, "0", "", "NULL BOUNDARY" ); // degrees of freedom, index, index as string, subspace, description
 
   m_list[0x1] = BoundaryPolytope(  ".5 .5 0 0 0 0; .5 .5 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1",
                                    ".5 .5 0 0 0 0; .5 .5 0 0 0 0; 0 0 1 0 0 0; 0 0 0 .5 .5 0; 0 0 0 .5 .5 0; 0 0 0 0 0 1",
                                    "0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1",
                                    5, 1, "1", "(r, r, s, t, u, v)", "Case 1: g1=g2, special position g4=g5" );
   
   m_list[0x2] = BoundaryPolytope(  "0 .5  .5 0 0 0; 0 .5 .5 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1",
                                    "1 0 0 0 0 0; 0 .5 .5 0 0 0; 0 .5 .5 0 0 0; 0 0 0 1 0 0; 0 0 0 0 .5 .5; 0 0 0 0 .5 .5",
                                    "1 0 0 0 0 0; 0 0 1 0 0 0; 0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1",
                                    5, 2, "2", "(r, s, s, t, u, v)", "Case 2: g2=g3, special position g5=g6" );

   m_list[0x3] = BoundaryPolytope(  "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0  1 0; 0 0 0 0 0  1",
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0",
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 -1 0; 0 0 0 0 0 -1",
                                    5, 3, "3", "(r, s, t, 0, -u, -v)", "Case 3: g4=0" );
   
   m_list[0x4] = BoundaryPolytope(  "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 0; 0 0 0 0 0 1",
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0",
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 -1",
                                     5, 4, "4", "(r, s, t, -u, 0, -v)", "Case 4: g5=0" );

   m_list[0x5] = BoundaryPolytope(  "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 0",
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 0",
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0; 0 0 0 0 0 1",
                                    5, 5, "5", "(r, s, t, -u, -v, 0)", "Case 5: g6=0" );

   m_list[0x6] = BoundaryPolytope(  "1 0 0 0 0 0; 0 .5 0 .5 0 0; 0 0 1 0 0 0; 0 .5 0 .5 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1",
                                    "1 0 0 0 0 0; 0 .5 0 .5 0 0; 0 0 1 0 0 0; 0 .5 0 .5 0 0; 0 0 0 0 .2 .4; 0 0 0 0 .4 .8",
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 1 1 -1 0 0; 0 -2 0 1 0 0; 0 0 0 0 -1 1; 0 0 0 0 0 -1",
                                    5, 6, "6", "(r, s, t, s, u + v, v)", "Case 6: g2 = g4, g5 > g6" );

   m_list[0x7] = BoundaryPolytope(  "1 0 0 0 0 0; 0 .5 0 .5 0 0; 0 0 1 0 0 0; 0 .5 0 .5 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1",
                                    "1 0 0 0 0 0; 0 .5 0 .5 0 0; 0 0 1 0 0 0; 0 .5 0 .5 0 0; 0 0 0 0 .2 .4; 0 0 0 0 .4 .8",
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 1 1 -1 0 0; 0 2 0 -1 0 0; 0 0 0 0 -1 1; 0 0 0 0 0 1",
                                    5, 7, "7", "(r, s, t, s, u, u + v)", "Case 7: g2 = g4, g5 < g6" );

   m_list[0x8] = BoundaryPolytope(  "1 0 0 0 0 0; 0 .5 0 -.5 0 0; 0 0 1 0 0 0; 0 -.5 0 .5 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1",
                                    "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                    "1 0 0 0 0 0; 0 1 0 0 0 0; 0 1 1 1 0 0; 0 2 0 1 0 0; 0 0 0 0 -1 -1; 0 0 0 0 0 -1",
                                    5, 8, "8", "(r, s, t, -s, -u, -v)", "Case 8: g2 = -g4" );

   m_list[0x9] = BoundaryPolytope( ".5 0 0 0 .5 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;  .5 0 0 0 .5 0; 0 0 0 0 0 0",
                                   ".5 0 0 0 .5 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 .2 0 .4; .5 0 0 0 .5 0; 0 0 0 .4 0 .8",
                                   "1 0 0 0 0 0; 0 1 0 0 0 0; 1 0 1 0 -1 0; 0 0 0 -1 0 1; -2 0 0 0 1 0; 0 0 0 0 0 -1",
                                    5, 9, "9", "(r, s, t, u + v, r, u)", "Case 9: g1 = g5, g4 > g6" );

   m_list[0xA] = BoundaryPolytope( ".5 0 0 0 .5 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0;  .5 0 0 0 .5 0; 0 0 0 0 0 0",
                                   ".5 0 0 0 .5 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 .2 0 .4; .5 0 0 0 .5 0; 0 0 0 .4 0 .8",
                                   "1 0 0 0 0 0; 0 1 0 0 0 0; 1 0 1 0 -1 0; 0 0 0 -1 0 1; 2 0 0 0 -1 0; 0 0 0 0 0 1",
                                    5, 10, "A", "(r, s, t, u, r, u + v)", "Case A: g1 = g5, g4 < g6" );

   m_list[0xB] = BoundaryPolytope( ".5 0 0 0 -.5 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; -.5 0 0 0 .5 0; 0 0 0 0 0 1",
                                   "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                   "1 0 0 0 0 0; 0 1 0 0 0 0; 1 0 1 0 1 0; 0 0 0 -1 0 -1; 2 0 0 0 1 0; 0 0 0 0 0 -1",
                                    5, 11, "B", "(r, s, t,-u,-r,-v)", "Case B: g1 = -g5" );

   m_list[0xC] = BoundaryPolytope( ".5 0 0 0 0 .5; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; .5 0 0 0 0 .5",
                                   ".5 0 0 0 0 .5; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 .2 .4 0; 0 0 0 .4 .8 0; .5 0 0 0 0 .5",
                                   "1 0 0 0 0 0; 1 1 0 0 0 -1; 0 0 1 0 0 0; 0 0 0 -1 1 0; 0 0 0 0 -1 0; -2 0 0 0 0 1",
                                    5, 12, "C", "(r, s, t, u + v, v, r)", "Case C: g1 = g6, g4 > g5" );   

   m_list[0xD] = BoundaryPolytope( ".5 0 0 0 0 .5; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; .5 0 0 0 0 .5",
                                   ".5 0 0 0 0 .5; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 .2 .4 0; 0 0 0 .4 .8 0; .5 0 0 0 0 .5",
                                   "1 0 0 0 0 0; 1 1 0 0 0 -1; 0 0 1 0 0 0; 0 0 0 -1 1 0; 0 0 0 0 1 0; 2 0 0 0 0 -1",
                                    5, 13, "D", "(r, s, t, u, u + v, r)", "Case D: g1 = g6, g4 < g5" );
   
   m_list[0xE] = BoundaryPolytope( ".5 0 0 0 0 -.5; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; -.5 0 0 0 0 .5",
                                   "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                   "1 0 0 0 0 0; 1 1 0 0 0 1; 0 0 1 0 0 0; 0 0 0 -1 -1 0; 0 0 0 0 -1 0; 2 0 0 0 0 1",
                                   5, 14, "E", "(r, s, t,-u,-v,-r)", "Case D: g1 = -g6" );

   m_list[0xF] = BoundaryPolytope( ".8 -.2 0 -.2 -.2 -.2; -.2 .8 0 -.2 -.2 -.2; 0 0 1 0 0 0; -.2 -.2 0 .8 -.2 -.2; -.2 -.2 0 -.2 .8 -.2; -.2 -.2 0 -.2 -.2 .8",
                                   ".55 .05 0 .05 -.45 -.2; .05 .55 0 -.45 .05 -.2; 0 0 1 0 0 0; .05 -.45 0 .55 .05 -.2; -.45 .05 0 .05 .55 -.2; -.2 -.2 0 -.2 -.2 .8",
                                   "1 0 0 0 0 0; 0 1 0 0 0 0; 1 1 1 1 1 1; 0 -2 0 -1 0 -1; -2 0 0 0 -1 -1; 0 0 0 0 0 1",
                                   5, 15, "F", "(r, s, t,-u,-v,-r - s + u + v)", "Case F: g1 + g2 + g3 + g4 + g5 + g6 = g3" );

   // Case678X g2=g4=g5=g6=0
   m_list[16] = BoundaryPolytope( "1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0.",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 1 0 0 0; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 0 0 0 1; 0 0 0 1 0 0; 0 0 0 0 1 0",
                                  1, 16, "16", "", "Case678X g2=g4=g5=g6=0" );

   // Case9ABCDEX g1=g4=g5=g6=0
   m_list[17] = BoundaryPolytope( "0 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0.",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 1 0 0 0 0; 0 0 1 0 0 0; 1 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 1 0 0",
                                  1, 17, "17", "", "Case9ABCDEX g1=g4=g5=g6=0" );

   // CaseFX g1=g2=g4=g5=g6=0
   m_list[18] = BoundaryPolytope( "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0.",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 1 0 0 0; 0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0; 0 0 0 1 0 0",
                                  1, 18, "18", "", "CaseFX g1=g2=g4=g5=g6=0" );

   // DATA boundary67
   m_list[19] = BoundaryPolytope( "1 0 0 0 0 0; 0 .5 0 .5 0 0; 0 0 1 0 0 0; 0 .5 0 .5 0 0; 0 0 0 0 .5 .5; 0 0 0 0 .5 .5",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY TRANSFORM
                                  4, 19, "19", "", "boundary67" );

   // DATA boundary9A
   m_list[20] = BoundaryPolytope( ".5 0 0 0 .5 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 .5 0 .5; .5 0 0 0 .5 0; 0 0 0 .5 0 .5",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY TRANSFORM
                                  4, 20, "20", "", "boundary9A" );

   // DATA boundaryCD
   m_list[21] = BoundaryPolytope( ".5 0 0 0 0 .5; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 .5 .5 0; 0 0 0 .5 .5 0; .5 0 0 0 0 .5",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY TRANSFORM
                                  4, 21, "21", "", "boundaryCD" );

   // Case boundary12 g1=g2=g3
   m_list[22] = BoundaryPolytope( "0.3333333333 0.3333333333 0.3333333333 0 0 0; 0.3333333333 0.3333333333 0.3333333333 0 0 0; 0.3333333333 0.3333333333 0.3333333333 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1.",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY TRANSFORM
                                  3, 22, "22", "", "boundary12 g1=g2=g3" );

//C      case boundary 8F g4=-g2, g1+g2+g4+g5+g6 = 0
   m_list[23] = BoundaryPolytope( "0.6666666667 0 0 0 -0.3333333333 -0.3333333333; 0 .5 0 -.5 0 0; 0 0 1. 0 0 0; 0 -.5 0 .5 0 0; -0.3333333333 0 0 0 0.6666666667 -0.3333333333; -0.3333333333 0 0 0 -0.3333333333 0.6666666667",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY TRANSFORM
                                  1, 23, "23", "", "boundary 8F g4=-g2, g1+g2+g4+g5+g6 = 0" );

//C      case boundary BF g5=-g1, g1+g2+g4+g5+g6 = 0
   m_list[24] = BoundaryPolytope( ".5 0 0 0 -.5 0; 0 0.6666666667 0 -0.3333333333 0 -0.3333333333; 0 0 1. 0 0 0; 0 -0.3333333333 0 0.6666666667 0 -0.3333333333; -.5 0 0 0 .5 0; 0 -0.3333333333 0 -0.3333333333 0 0.6666666667",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY TRANSFORM
                                  1, 24, "24", "", "boundary BF g4=-g2, g1+g2+g4+g5+g6 = 0" );

//C      case boundary EF g6=-g1, g1+g2+g4+g5+g6 = 0
   m_list[25] = BoundaryPolytope( "0.5 0 0 0 0 -0.5; 0 0.6666666667 0 -0.3333333333 -0.3333333333 0; 0 0 1. 0 0 0; 0 -0.3333333333 0 0.6666666667 -0.3333333333 0; 0 -0.3333333333 0 -0.3333333333 0.6666666667 0; -0.5 0 0 0 0 0.5",
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY PRJHAT
                                  "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0", // NULL TEMPORARY TRANSFORM
                                  1, 25, "25", "", "boundary EF g6=-g1, g1+g2+g4+g5+g6 = 0" );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: ~BoundaryPolytopeList()
// Description: destructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
BoundaryPolytopeList::~BoundaryPolytopeList( void )
{
   m_list.clear();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: operator[]()
// Description: accessor for the n-th polytope
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
const BoundaryPolytope& BoundaryPolytopeList::operator[] ( const int n ) const
{
   int nn = n;
   if ( n < 0 || n > size()-1 ) nn = 0;
   return( m_list[nn] ); 
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: DegreesofFreedom()
// Description: I BELIEVE THIS IS NOT CORRECT !!!!!!!!!!!!!!!!!!!!!!!!!
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int BoundaryPolytopeList::DegreesofFreedom( const int ITDESG )
{
   static const
      int idf[] = { 0,1,0,1,0,1,1,2,1,3,
      1,1,2,3,1,2,3,1,2,3,
      1,1,2,1,3,2,3,3,3,3,
      5,2,3,3,3,2,3,2,3,2,
      3,2,3,5 };
   return( idf[std::max(1,std::min(abs(ITDESG),44))] );                            
}

