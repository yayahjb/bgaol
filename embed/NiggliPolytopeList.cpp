#include "NiggliPolytopeList.h"

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
   size_t size( void ) const                           == returns the number of items in the List
   const NiggliPolytope& operator[] ( const int n ) == returns the n-th item in the list
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
NiggliPolytopeList::NiggliPolytopeList(void)
   : m_list(45)
{
   m_list[00] = NiggliPolytope( "invalid", 0, "invalid", "", "invalid", 1, "0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0" ); // null invalid polytope

   m_list[01] = NiggliPolytope( "cF",  1, "44C", "(r, r, r, r, r, r)",              "face-centered cubic",          6, "1 1 1 1 1 1;      1 1 1 1 1 1;       1 1 1 1 1 1;       1 1 1 1 1 1;           1 1 1 1 1 1;        1 1 1 1 1 1" );
   m_list[02] = NiggliPolytope( "hR",  2, "49C", "(r, r, r, s, s, s)",              "rhomobohedral",                3, "1 1 1 0 0 0;      1 1 1 0 0 0;       1 1 1 0 0 0;       0 0 0 1 1 1;           0 0 0 1 1 1;        0 0 0 1 1 1" );
   m_list[03] = NiggliPolytope( "cP",  3, "44A", "(r, r, r, 0, 0, 0)",              "primitive cubic",              3, "1 1 1 0 0 0;      1 1 1 0 0 0;       1 1 1 0 0 0;       0 0 0 0 0 0;           0 0 0 0 0 0;        0 0 0 0 0 0" );
   m_list[04] = NiggliPolytope( "hR",  4, "49D", "(r, r, r, -s, -s, -s)",           "rhomobohedral",                3, "1 1 1 0 0 0;      1 1 1 0 0 0;       1 1 1 0 0 0;       0 0 0 1 1 1;           0 0 0 1 1 1;        0 0 0 1 1 1" );
   m_list[05] = NiggliPolytope( "cI",  5, "44B", "(r, r, r, -2r/3, -2r/3, -2r/3)",  "body-centered cubic",         34, "9 9 9 -6 -6 -6;   9 9 9 -6 -6 -6;    9 9 9 -6 -6       -6; -6 -6 -6 4 4 4;    -6 -6 -6 4 4 4;     -6 -6 -6 4 4 4" );
   m_list[06] = NiggliPolytope( "tI",  6, "45D", "(r, r, r, -r + s, -r + s, -2s)",  "body-centered tetragonal",    26, "6 6 6 -4 -4 -4;   6 6 6 -4 -4 -4;    6 6 6 -4 -4       -4; -4 -4 -4 7 7 -6;   -4 -4 -4 7 7 -6;    -4 -4 -4 -6 -6 20" );
   m_list[07] = NiggliPolytope( "tI",  7, "45D", "(r, r, r, -2s, -r + s, -r + s)",  "body-centered tetragonal",    26, "6 6 6 -4 -4 -4;   6 6 6 -4 -4 -4;    6 6 6 -4 -4       -4; -4 -4 -4 20 -6 -6; -4 -4 -4 -6 7 7;    -4 -4 -4 -6 7 7" );
   m_list[ 8] = NiggliPolytope( "oI",  8, "52A", "(r, r, r, -s, -t, -2r + s + t)",  "body-centered orthorhombic",  13, "3 3 3 -2 -2 -2;   3 3 3 -2 -2 -2;    3 3 3 -2 -2       -2; -2 -2 -2 10 -3 -3; -2 -2 -2 -3 10 -3;  -2 -2 -2 -3 -3 10" );
   m_list[ 9] = NiggliPolytope( "hR",  9, "49B", "(r, r, s, r, r, r)",              "rhomobohedral",                5, "1 1 0 1 1 1;      1 1 0 1 1 1;       0 0 5 0 0 0;       1 1 0 1 1 1;           1 1 0 1 1 1;        1 1 0 1 1 1" );
   m_list[10] = NiggliPolytope( "mC", 10, "55A", "(r, r, s, t, t, u)",              "C-centered monoclinic",        2, "1 1 0 0 0 0;      1 1 0 0 0 0;       0 0 2 0 0 0;       0 0 0 1 1 0;           0 0 0 1 1 0;        0 0 0 0 0 2" );
   m_list[11] = NiggliPolytope( "tP", 11, "45A", "(r, r, s, 0, 0, 0)",              "primitive tetragonal",         2, "1 1 0 0 0 0;      1 1 0 0 0 0;       0 0 2 0 0 0;       0 0 0 0 0 0;           0 0 0 0 0 0;        0 0 0 0 0 0" );
   m_list[12] = NiggliPolytope( "hP", 12, "48A", "(r, r, s, 0, 0, -r)",             "primitive hexagonal",          3, "1 1 0 0 0 -1;     1 1 0 0 0 -1;      0 0 3 0 0 0;       0 0 0 0 0 0;           0 0 0 0 0 0;       -1 -1 0 0 0 1" );
   m_list[13] = NiggliPolytope( "oC", 13, "50D", "(r, r, s, 0, 0, -t)",             "body-centered orthorhombic",   2, "1 1 0 0 0 0;      1 1 0 0 0 0;       0 0 2 0 0 0;       0 0 0 0 0 0;           0 0 0 0 0 0;        0 0 0 0 0 2" );
   m_list[14] = NiggliPolytope( "mC", 14, "55A", "(r, r, s, t, t, u)",              "C-centered monoclinic",        2, "1 1 0 0 0 0;      1 1 0 0 0 0;       0 0 2 0 0 0;       0 0 0 1 1 0;           0 0 0 1 1 0;        0 0 0 0 0 2" );
   m_list[15] = NiggliPolytope( "tI", 15, "45C", "(r, r, s, -r, -r, 0)",            "body-centered tetragonal",     4, "1 1 0 -1 -1 0;    1 1 0 -1 -1 0;     0 0 4 0 0 0;      -1 -1 0 1 1 0;         -1 -1 0 1 1 0;       0 0 0 0 0 0" );
   m_list[16] = NiggliPolytope( "oF", 16, "51A", "(r, r, s, -t, -t, -2r + 2t)",     "face-centered orthorhombic",  10, "3 3 0 -2 -2 -2;   3 3 0 -2 -2 -2;    0 0 10 0 0 0;     -2 -2 0 3 3 -2;        -2 -2 0 3 3 -2;     -2 -2 0 -2 -2 8" );
   m_list[17] = NiggliPolytope( "mC", 17, "57B", "(r, r, s, -t, -u, -2r + t + u)",  "C-centered monoclinic",       10, "3 3 0 -2 -2 -2;   3 3 0 -2 -2 -2;    0 0 10 0 0 0;     -2 -2 0 8 -2 -2;       -2 -2 0 -2 8 -2;    -2 -2 0 -2 -2 8" );
   m_list[18] = NiggliPolytope( "tI", 18, "45E", "(r, s, s, r/2, r, r)",            "body-centered tetragonal",    26, "8 0 0 4 8 8;      0 13 13 0 0 0;     0 13 13 0 0 0;     4 0 0 2 4 4;           8 0 0 4 8 8;        8 0 0 4 8 8" );
   m_list[19] = NiggliPolytope( "oI", 19, "52B", "(r, s, s, t, r, r)",              "body-centered orthorhombic",   6, "2 0 0 0 2 2;      0 3 3 0 0 0;       0 3 3 0 0 0;       0 0 0 6 0 0;           2 0 0 0 2 2;        2 0 0 0 2 2" );
   m_list[20] = NiggliPolytope( "mC", 20, "55B", "(r, s, s, t, u, u)",              "C-centered monoclinic",        2, "2 0 0 0 0 0;      0 1 1 0 0 0;       0 1 1 0 0 0;       0 0 0 2 0 0;           0 0 0 0 1 1;        0 0 0 0 1 1" );
   m_list[21] = NiggliPolytope( "tP", 21, "45B", "(r, s, s, 0, 0, 0)",              "primitive tetragonal",         2, "2 0 0 0 0 0;      0 1 1 0 0 0;       0 1 1 0 0 0;       0 0 0 0 0 0;           0 0 0 0 0 0;        0 0 0 0 0 0" );
   m_list[22] = NiggliPolytope( "hP", 22, "48B", "(r, s, s, -s, 0, 0)",             "primitive hexagonal",          3, "3 0 0 0 0 0;      0 1 1 -1 0 0;      0 1 1 -1 0 0;      0 -1 -1 1 0 0;         0 0 0 0 0 0;        0 0 0 0 0 0" );
   m_list[23] = NiggliPolytope( "oC", 23, "50E", "(r, s, s, -t, 0, 0)",             "C-centered orthorhombic",      2, "2 0 0 0 0 0;      0 1 1 0 0 0;       0 1 1 0 0 0;       0 0 0 2 0 0;           0 0 0 0 0 0;        0 0 0 0 0 0" );
   m_list[24] = NiggliPolytope( "hR", 24, "49E", "(r, s, s, -s+r/3, -2r/3, -2r/3)", "rhomobohedral",               53, "27 3 3 6 -18 -18; 3 18 18 -17 -2 -2; 3 18 18 -17 -2 -2; 6 -17 -17 19 -4 -4;  -18 -2 -2 -4 12 12; -18 -2 -2 -4 12 12" );
   m_list[25] = NiggliPolytope( "mC", 25, "55B", "(r, s, s, t, u, u)",              "C-centered monoclinic",        2, "2 0 0 0 0 0;      0 1 1 0 0 0;       0 1 1 0 0 0;       0 0 0 2 0 0;           0 0 0 0 1 1;        0 0 0 0 1 1" );
   m_list[26] = NiggliPolytope( "oF", 26, "51B", "(r, s, t, r/2, r, r)",            "face-centered corthorhombic", 13, "4 0 0 2 4 4;      0 13 0 0 0 0;      0 0 13 0 0 0;      2 0 0 1 2 2;           4 0 0 2 4 4;        4 0 0 2 4 4" );
   m_list[27] = NiggliPolytope( "mC", 27, "57C", "(r, s, t, u, r, r)",              "C-centered monoclinic",        3, "1 0 0 0 1 1;      0 3 0 0 0 0;       0 0 3 0 0 0;       0 0 0 3 0 0;           1 0 0 0 1 1;        1 0 0 0 1 1" );
   m_list[28] = NiggliPolytope( "mC", 28, "56A", "(r, s, t, u, r, 2u)",             "C-centered monoclinic",       10, "5 0 0 0 5 0;      0 10 0 0 0 0;      0 0 10 0 0 0;      0 0 0 2 0 4;           5 0 0 0 5 0;        0 0 0 4 0 8" );
   m_list[29] = NiggliPolytope( "mC", 29, "56C", "(r, s, t, u, 2u, r)",             "C-centered monoclinic",       10, "5 0 0 0 0 5;      0 10 0 0 0 0;      0 0 10 0 0 0;      0 0 0 2 4 0;           0 0 0 4 8 0;        5 0 0 0 0 5" );
   m_list[30] = NiggliPolytope( "mC", 30, "56B", "(r, s, t, s, u, 2u)",             "C-centered monoclinic",       10, "10 0 0 0 0 0;     0 5 0 5 0 0;       0 0 10 0 0 0;      0 5 0 5 0 0;           0 0 0 0 2 4;        0 0 0 0 4 8" );
   m_list[31] = NiggliPolytope( "aP", 31, "58S", "(r, s, t, u, v, w)",              "anorthic",                     1, "1 0 0 0 0 0;      0 1 0 0 0 0;       0 0 1 0 0 0;       0 0 0 1 0 0;           0 0 0 0 1 0;        0 0 0 0 0 1" );
   m_list[32] = NiggliPolytope( "oP", 32, "50C", "(r, s, t, 0, 0, 0)",              "body-centered orthorhombic",   1, "1 0 0 0 0 0;      0 1 0 0 0 0;       0 0 1 0 0 0;       0 0 0 0 0 0;           0 0 0 0 0 0;        0 0 0 0 0 0" );
   m_list[33] = NiggliPolytope( "mP", 33, "53A", "(r, s, t, 0, -u, 0)",             "primitive monoclinic",         1, "1 0 0 0 0 0;      0 1 0 0 0 0;       0 0 1 0 0 0;       0 0 0 0 0 0;           0 0 0 0 1 0;        0 0 0 0 0 0" );
   m_list[34] = NiggliPolytope( "mP", 34, "53C", "(r, s, t, 0, 0, -u)",             "primitive monoclinic",         1, "1 0 0 0 0 0;      0 1 0 0 0 0;       0 0 1 0 0 0;       0 0 0 0 0 0;           0 0 0 0 0 0;        0 0 0 0 0 1" );
   m_list[35] = NiggliPolytope( "mP", 35, "53B", "(r, s, t, -u, 0, 0)",             "primitive monoclinic",         1, "1 0 0 0 0 0;      0 1 0 0 0 0;       0 0 1 0 0 0;       0 0 0 1 0 0;           0 0 0 0 0 0;        0 0 0 0 0 0" );
   m_list[36] = NiggliPolytope( "oC", 36, "50A", "(r, s, t, 0, -r, 0)",             "body-centered orthorhombic",   2, "1 0 0 0 -1 0;     0 2 0 0 0 0;       0 0 2 0 0 0;       0 0 0 0 0 0;          -1 0 0 0 1 0;        0 0 0 0 0 0" );
   m_list[37] = NiggliPolytope( "mC", 37, "54C", "(r, s, t, -u, -r, 0)",            "C-centered monoclinic",        2, "1 0 0 0 -1 0;     0 2 0 0 0 0;       0 0 2 0 0 0;       0 0 0 2 0 0;          -1 0 0 0 1 0;        0 0 0 0 0 0" );
   m_list[38] = NiggliPolytope( "oC", 38, "50B", "(r, s, t, 0, 0, -r)",             "C-centered orthorhombic",      2, "1 0 0 0 0 -1;     0 2 0 0 0 0;       0 0 2 0 0 0;       0 0 0 0 0 0;           0 0 0 0 0 0;       -1 0 0 0 0 1" );
   m_list[39] = NiggliPolytope( "mC", 39, "54A", "(r, s, t, -u, 0, -r)",            "C-centered monoclinic",        2, "1 0 0 0 0 -1;     0 2 0 0 0 0;       0 0 2 0 0 0;       0 0 0 2 0 0;           0 0 0 0 0 0;       -1 0 0 0 0 1" );
   m_list[40] = NiggliPolytope( "oC", 40, "50F", "(r, s, t, -s, 0, 0)",             "C-centered orthorhombic" ,     2, "2 0 0 0 0 0;      0 1 0 -1 0 0;      0 0 2 0 0 0;       0 -1 0 1 0 0;          0 0 0 0 0 0;        0 0 0 0 0 0" );
   m_list[41] = NiggliPolytope( "mC", 41, "54B", "(r, s, t, -s, -u, 0)",            "C-centered monoclinic",        2, "2 0 0 0 0 0;      0 1 0 -1 0 0;      0 0 2 0 0 0;       0 -1 0 1 0 0;          0 0 0 0 2 0;        0 0 0 0 0 0" );
   m_list[42] = NiggliPolytope( "oI", 42, "52C", "(r, s, t, -s, -r, 0)",            "body-centered orthorhombic",   2, "1 0 0 0 -1 0;     0 1 0 -1 0 0;      0 0 2 0 0 0;       0 -1 0 1 0 0;         -1 0 0 0 1 0;        0 0 0 0 0 0" );
   m_list[43] = NiggliPolytope( "mC", 43, "57A", "(r, s, t, -s + u, -r + u, -2u)",  "C-centered monoclinic",       20, "11 1 0 1 -9 -4;   1 11 0 -9 1 -4;    0 0 20 0 0 0;      1 -9 0 11 1 -4;       -9 1 0 1 11 -4;     -4 -4 0 -4 -4 16" );
   m_list[44] = NiggliPolytope( "aP", 44, "58B", "(r, s, t, u, v, w)",              "anorthic",                     1, "1 0 0 0 0 0;      0 1 0 0 0 0;       0 0 1 0 0 0;       0 0 0 1 0 0;           0 0 0 0 1 0;        0 0 0 0 0 1" );
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: ~NiggliPolytopeList()
// Description: destructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
NiggliPolytopeList::~NiggliPolytopeList(void)
{
   m_list.clear();
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// Name: operator[]()
// Description: accessor for the n-th polytope
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
const NiggliPolytope& NiggliPolytopeList::operator[] ( const int n ) const
{
   int nn = n;
   if ( n < 0 || n > size()-1 ) nn = 0;
   return( m_list[nn] ); 
}
