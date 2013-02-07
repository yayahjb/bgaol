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
   int size( void ) const                           == returns the number of items in the List
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

/*

      data ngtype(1) /-3/
      data lattyp(1) /'cP'/
      data pjnorm(1) /3D0/
      data( projct(i,1),i=1,36) /
     1.0/3.0 1.0/3.0 1.0/3.0 0 0 0 ; 1.0/3.0 1.0/3.0 1.0/3.0 0 0 0 ; 1.0/3.0 1.0/3.0 1.0/3.0 0 0 0 ; 0 0 0 0 0 0 ; 0 0 0 0 0 0 ; 0 0 0 0 0 0


      data ngtype(2) /-5/
      data lattyp(2) /'cI'/
data pjnorm(2) /39D0/
data (projct(i,2),i=1,36) /

"9/39.0  9/39.0  9/39.0 -6/39.0 -6/39.0 -6/39.0; -6/39.0 -6/39.0 -6/39.0  4/39.0  4/39.0  4/39.0; 9/39.0  9/39.0  9/39.0 -6/39.0 -6/39.0 -6/39.0; 9/39.0  9/39.0  9/39.0 -6/39.0 -6/39.0 -6/39.0; -6/39.0 -6/39.0 -6/39.0  4/39.0  4/39.0  4/39.0; -6/39.0 -6/39.0 -6/39.0  4/39.0  4/39.0  4/39.0"


      data ngtype(3) /1/
      data lattyp(3) /'cF'/
      data pjnorm(3) /6D0/
      data (projct(i,3),i=1,36) /
       1,1,1,1,1,1,
       1,1,1,1,1,1,
       1,1,1,1,1,1,
       1,1,1,1,1,1,
       1,1,1,1,1,1,
       1,1,1,1,1,1 /



      data ngtype(4) /-11/
      data lattyp(4) /'tP'/
      data pjnorm(4) /2D0/
      data (projct(i,4),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(5) /-21/
      data lattyp(5) /'tP'/
      data pjnorm(5) /2D0/
      data (projct(i,5),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(6) /-15/
      data lattyp(6) /'tI'/
      data pjnorm(6) /4D0/
      data (projct(i,6),i=1,36) /
     1  1, 1, 0,-1,-1, 0,
     2  1, 1, 0,-1,-1, 0,
     3  0, 0, 4, 0, 0, 0,
     4 -1,-1, 0, 1, 1, 0,
     5 -1,-1, 0, 1, 1, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(7) /-6/
      data lattyp(7) /'tI'/
      data pjnorm(7) /26D0/
      data (projct(i,7),i=1,36) /
     1  6, 6, 6,-4,-4,-4,
     2  6, 6, 6,-4,-4,-4,
     3  6, 6, 6,-4,-4,-4,
     4 -4,-4,-4, 7, 7,-6,
     5 -4,-4,-4, 7, 7,-6,
     6 -4,-4,-4,-6,-6,20 /


      data ngtype(8) /-7/
      data lattyp(8) /'tI'/
      data pjnorm(8) /26D0/
      data (projct(i,8),i=1,36) /
     1  6, 6, 6,-4,-4,-4,
     2  6, 6, 6,-4,-4,-4,
     3  6, 6, 6,-4,-4,-4,
     4 -4,-4,-4,20,-6,-6,
     5 -4,-4,-4,-6, 7, 7,
     6 -4,-4,-4,-6, 7, 7 /


      data ngtype(9) /18/
      data lattyp(9) /'tI'/
      data pjnorm(9) /26D0/
      data (projct(i,9),i=1,36) /
     1  8, 0, 0, 4, 8, 8,
     2  0,13,13, 0, 0, 0,
     3  0,13,13, 0, 0, 0,
     4  4, 0, 0, 2, 4, 4,
     5  8, 0, 0, 4, 8, 8,
     6  8, 0, 0, 4, 8, 8 /


      data ngtype(10) /-12/
      data lattyp(10) /'hP'/
      data pjnorm(10) /3D0/
      data (projct(i,10),i=1,36) /
     1  1, 1, 0, 0, 0,-1,
     2  1, 1, 0, 0, 0,-1,
     3  0, 0, 3, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1,-1, 0, 0, 0, 1 /


      data ngtype(11) /-22/
      data lattyp(11) /'hP'/
      data pjnorm(11) /3D0/
      data (projct(i,11),i=1,36) /
     1  3, 0, 0, 0, 0, 0,
     2  0, 1, 1,-1, 0, 0,
     3  0, 1, 1,-1, 0, 0,
     4  0,-1,-1, 1, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(12) /9/
      data lattyp(12) /'hR'/
      data pjnorm(12) /5D0/
      data (projct(i,12),i=1,36) /
     1  1,1,0,1,1,1,
     2  1,1,0,1,1,1,
     3  0,0,5,0,0,0,
     4  1,1,0,1,1,1,
     5  1,1,0,1,1,1,
     6  1,1,0,1,1,1  /


      data ngtype(13) /2/
      data lattyp(13) /'hR'/
      data pjnorm(13) /3D0/
      data (projct(i,13),i=1,36) /
     1  1,1,1,0,0,0,
     2  1,1,1,0,0,0,
     3  1,1,1,0,0,0,
     4  0,0,0,1,1,1,
     5  0,0,0,1,1,1,
     6  0,0,0,1,1,1 /


      data ngtype(14) /-4/
      data lattyp(14) /'hR'/
      data pjnorm(14) /3D0/
      data (projct(i,14),i=1,36) /
     1  1,1,1,0,0,0,
     2  1,1,1,0,0,0,
     3  1,1,1,0,0,0,
     4  0,0,0,1,1,1,
     5  0,0,0,1,1,1,
     6  0,0,0,1,1,1 /


      data ngtype(15) /-24/
      data lattyp(15) /'hR'/
      data pjnorm(15) /53D0/
      data (projct(i,15),i=1,36) /
     1  27,  3,  3,  6,-18,-18,
     2   3, 18, 18,-17, -2, -2,
     3   3, 18, 18,-17, -2, -2,
     4   6,-17,-17, 19, -4, -4,
     5 -18, -2, -2, -4, 12, 12,
     6 -18, -2, -2, -4, 12, 12 /


      data ngtype(16) /-32/
      data lattyp(16) /'oP'/
      data pjnorm(16) /1D0/
      data (projct(i,16),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(17) /-36/
      data lattyp(17) /'oC'/
      data pjnorm(17) /2D0/
      data (projct(i,17),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0 /


      data ngtype(18) /-38/
      data lattyp(18) /'oC'/
      data pjnorm(18) /2D0/
      data (projct(i,18),i=1,36) /
     1  1, 0, 0, 0, 0,-1,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 0, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1, 0, 0, 0, 0, 1 /


      data ngtype(19) /-13/
      data lattyp(19) /'oC'/
      data pjnorm(19) /2D0/
      data (projct(i,19),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,2 /


      data ngtype(20) /-23/
      data lattyp(20) /'oC'/
      data pjnorm(20) /2D0/
      data (projct(i,20),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0 /


      data ngtype(21) /-40/
      data lattyp(21) /'oC'/
      data pjnorm(21) /2D0/
      data (projct(i,21),i=1,36) /
     1  2, 0, 0, 0, 0, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(22) /-16/
      data lattyp(22) /'oF'/
      data pjnorm(22) /10D0/
      data (projct(i,22),i=1,36) /
     1  3, 3, 0,-2,-2,-2,
     2  3, 3, 0,-2,-2,-2,
     3  0, 0,10, 0, 0, 0,
     4 -2,-2, 0, 3, 3,-2,
     5 -2,-2, 0, 3, 3,-2,
     6 -2,-2, 0,-2,-2, 8  /


      data ngtype(23) /26/
      data lattyp(23) /'oF'/
      data pjnorm(23) /13D0/
      data (projct(i,23),i=1,36) /
     1  4, 0, 0, 2, 4, 4,
     2  0,13, 0, 0, 0, 0,
     3  0, 0,13, 0, 0, 0,
     4  2, 0, 0, 1, 2, 2,
     5  4, 0, 0, 2, 4, 4,
     6  4, 0, 0, 2, 4, 4  /


      data ngtype(24) /-8/
      data lattyp(24) /'oI'/
      data pjnorm(24) /13D0/
      data (projct(i,24),i=1,36) /
     1  3, 3, 3,-2,-2,-2,
     2  3, 3, 3,-2,-2,-2,
     3  3, 3, 3,-2,-2,-2,
     4 -2,-2,-2,10,-3,-3,
     5 -2,-2,-2,-3,10,-3,
     6 -2,-2,-2,-3,-3,10  /


      data ngtype(25) /19/
      data lattyp(25) /'oI'/
      data pjnorm(25) /6D0/
      data (projct(i,25),i=1,36) /
     1  2,0,0,0,2,2,
     2  0,3,3,0,0,0,
     3  0,3,3,0,0,0,
     4  0,0,0,6,0,0,
     5  2,0,0,0,2,2,
     6  2,0,0,0,2,2  /


      data ngtype(26) /-42/
      data lattyp(26) /'oI'/
      data pjnorm(26) /2D0/
      data (projct(i,26),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(27) /-33/
      data lattyp(27) /'mP'/
      data pjnorm(27) /1D0/
      data (projct(i,27),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,1,0,
     6  0,0,0,0,0,0  /


      data ngtype(28) /-35/
      data lattyp(28) /'mP'/
      data pjnorm(28) /1D0/
      data (projct(i,28),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,1,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,0  /


      data ngtype(29) /-34/
      data lattyp(29) /'mP'/
      data pjnorm(29) /1D0/
      data (projct(i,29),i=1,36) /
     1  1,0,0,0,0,0,
     2  0,1,0,0,0,0,
     3  0,0,1,0,0,0,
     4  0,0,0,0,0,0,
     5  0,0,0,0,0,0,
     6  0,0,0,0,0,1  /


      data ngtype(30) /-39/
      data lattyp(30) /'mC'/
      data pjnorm(30) /2D0/
      data (projct(i,30),i=1,36) /
     1  1, 0, 0, 0, 0,-1,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 2, 0, 0,
     5  0, 0, 0, 0, 0, 0,
     6 -1, 0, 0, 0, 0, 1  /


      data ngtype(31) /-41/
      data lattyp(31) /'mC'/
      data pjnorm(31) /2D0/
      data (projct(i,31),i=1,36) /
     1  2, 0, 0, 0, 0, 0,
     2  0, 1, 0,-1, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0,-1, 0, 1, 0, 0,
     5  0, 0, 0, 0, 2, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(32) /-37/
      data lattyp(32) /'mC'/
      data pjnorm(32) /2D0/
      data (projct(i,32),i=1,36) /
     1  1, 0, 0, 0,-1, 0,
     2  0, 2, 0, 0, 0, 0,
     3  0, 0, 2, 0, 0, 0,
     4  0, 0, 0, 2, 0, 0,
     5 -1, 0, 0, 0, 1, 0,
     6  0, 0, 0, 0, 0, 0  /


      data ngtype(33) /10/
      data lattyp(33) /'mC'/
      data pjnorm(33) /2D0/
      data (projct(i,33),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,1,1,0,
     5  0,0,0,1,1,0,
     6  0,0,0,0,0,2  /


      data ngtype(34) /-14/
      data lattyp(34) /'mC'/
      data pjnorm(34) /2D0/
      data (projct(i,34),i=1,36) /
     1  1,1,0,0,0,0,
     2  1,1,0,0,0,0,
     3  0,0,2,0,0,0,
     4  0,0,0,1,1,0,
     5  0,0,0,1,1,0,
     6  0,0,0,0,0,2  /


      data ngtype(35) /20/
      data lattyp(35) /'mC'/
      data pjnorm(35) /2D0/
      data (projct(i,35),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,1,1,
     6  0,0,0,0,1,1  /


      data ngtype(36) /-25/
      data lattyp(36) /'mC'/
      data pjnorm(36) /2D0/
      data (projct(i,36),i=1,36) /
     1  2,0,0,0,0,0,
     2  0,1,1,0,0,0,
     3  0,1,1,0,0,0,
     4  0,0,0,2,0,0,
     5  0,0,0,0,1,1,
     6  0,0,0,0,1,1  /


      data ngtype(37) /28/
      data lattyp(37) /'mC'/
      data pjnorm(37) /10D0/
      data (projct(i,37),i=1,36) /
     1  5, 0, 0, 0, 5, 0,
     2  0,10, 0, 0, 0, 0,
     3  0, 0,10, 0, 0, 0,
     4  0, 0, 0, 2, 0, 4,
     5  5, 0, 0, 0, 5, 0,
     6  0, 0, 0, 4, 0, 8  /


      data ngtype(38) /30/
      data lattyp(38) /'mC'/
      data pjnorm(38) /10D0/
      data (projct(i,38),i=1,36) /
     1  10, 0, 0, 0, 0, 0,
     2   0, 5, 0, 5, 0, 0,
     3   0, 0,10, 0, 0, 0,
     4   0, 5, 0, 5, 0, 0,
     5   0, 0, 0, 0, 2, 4,
     6   0, 0, 0, 0, 4, 8  /


      data ngtype(39) /29/
      data lattyp(39) /'mC'/
      data pjnorm(39) /10D0/
      data (projct(i,39),i=1,36) /
     1  5, 0, 0, 0, 0, 5,
     2  0,10, 0, 0, 0, 0,
     3  0, 0,10, 0, 0, 0,
     4  0, 0, 0, 2, 4, 0,
     5  0, 0, 0, 4, 8, 0,
     6  5, 0, 0, 0, 0, 5  /


      data ngtype(40) /-43/
      data lattyp(40) /'mI'/
      data pjnorm(40) /20D0/
      data (projct(i,40),i=1,36) /
     1  11, 1, 0, 1,-9,-4,
     2   1,11, 0,-9, 1,-4,
     3   0, 0,20, 0, 0, 0,
     4   1,-9, 0,11, 1,-4,
     5  -9, 1, 0, 1,11,-4,
     6  -4,-4, 0,-4,-4,16  /


      data ngtype(41) /-17/
      data lattyp(41) /'mI'/
      data pjnorm(41) /10D0/
      data (projct(i,41),i=1,36) /
     1  3, 3, 0,-2,-2,-2,
     2  3, 3, 0,-2,-2,-2,
     3  0, 0,10, 0, 0, 0,
     4 -2,-2, 0, 8,-2,-2,
     5 -2,-2, 0,-2, 8,-2,
     6 -2,-2, 0,-2,-2, 8  /


      data ngtype(42) /27/
      data lattyp(42) /'mI'/
      data pjnorm(42) /3D0/
      data (projct(i,42),i=1,36) /
     1  1,0,0,0,1,1,
     2  0,3,0,0,0,0,
     3  0,0,3,0,0,0,
     4  0,0,0,3,0,0,
     5  1,0,0,0,1,1,
     6  1,0,0,0,1,1  /
     */
