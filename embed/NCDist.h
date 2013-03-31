//
//  CNCDist.c
//
//
//  Created by Herbert J. Bernstein on 3/26/13.
//
//

/* The projectors for the 15 base types (5-D boundaries
 in G6), plus a few extra for internal boundaries
 Note that the array indices are swapped from the
 Fortan versions */


#include <cmath>

#define NCD_min(a,b) (a<b?a:b)

static double prj[25][36]= {
    /* 1 */
    {0.5,0.5,0.0,0.0,0.0,0.0,
        0.5,0.5,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 2 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.5,0.0,0.0,0.0,
        0.0,0.5,0.5,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 3 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 4 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 5 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 6 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 7 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 8 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,-0.5,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,-0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 9 */
    {0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* A */
    {0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,1},
    /* B */
    {0.5,0.0,0.0,0.0,-0.5,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        -0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* C */
    {0.5,0.0,0.0,0.0,0.0,0.5,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.5,0.0,0.0,0.0,0.0,0.5},
    /* D */
    {0.5,0.0,0.0,0.0,0.0,0.5,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.5,0.0,0.0,0.0,0.0,0.5},
    /* E */
    {0.5,0.0,0.0,0.0,0.0,-0.5,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        -0.5,0.0,0.0,0.0,0.0,0.5},
    /* F */
    {0.8,-0.2,0.0,-0.2,-0.2,-0.2,
        -0.2,0.8,0.0,-0.2,-0.2,-0.2,
        0.0,0.0,1.0,0.0,0.0,0.0,
        -0.2,-0.2,0.0,0.8,-0.2,-0.2,
        -0.2,-0.2,0.0,-0.2,0.8,-0.2,
        -0.2,-0.2,0.0,-0.2,-0.2,0.8},
    /* case678X g2=g4=g5=g6=0  */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* case9ABCDEX g1=g4=g5=g6=0 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* caseFX g1=g2=g4=g5=g6=0 */
    {0.5,-0.5,0.0,0.0,0.0,0.0,
        -0.5,0.5,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 67 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.5,0.5,
        0.0,0.0,0.0,0.0,0.5,0.5},
    /* 9A */
    {0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.5,0.0,0.5,
        0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.5,0.0,0.5},
    /* CD */
    {0.5,0.0,0.0,0.0,0.0,0.5,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.5,0.5,0.0,
        0.0,0.0,0.0,0.5,0.5,0.0,
        0.5,0.0,0.0,0.0,0.0,0.5},
    /* 12 g1=g2=g3 */
    {.3333333333333333,.3333333333333333,.3333333333333333,0.0,0.0,0.0,
        .3333333333333333,.3333333333333333,.3333333333333333,0.0,0.0,0.0,
        .3333333333333333,.3333333333333333,.3333333333333333,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 8F g4=-g2, g1+g2+g4+g5+g6 = 0 */
    {.6666666666666667,0.0,0.0,0.0,-.3333333333333333,-.3333333333333333,
        0.0,0.5,0.0,-0.5,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,-0.5,0.0,0.5,0.0,0.0,
        -.3333333333333333,0.0,0.0,0.0,.6666666666666667,-.3333333333333333,
        -.3333333333333333,0.0,0.0,0.0,-.3333333333333333,.6666666666666667},
    /* BF g5=-g1.0, g1+g2+g4+g5+g6 = 0 */
    {0.5,0.0,0.0,0.0,-0.5,0.0,
        0.0,.6666666666666667,0.0,-.3333333333333333,0.0,-.3333333333333333,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,-.3333333333333333,0.0,.6666666666666667,0.0,-.3333333333333333,
        -0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,-.3333333333333333,0.0,-.3333333333333333,0.0,.6666666666666667},
    /* EF g6=-g1.0, g1+g2+g4+g5+g6 = 0 */
    {0.5,0.0,0.0,0.0,0.0,-0.5,
        0.0,.6666666666666667,0.0,-.3333333333333333,-.3333333333333333,0.0,
        0.0,0.0,1.0,0.0,0.0,0.0,
        0.0,-.3333333333333333,0.0,.6666666666666667,-.3333333333333333,0.0,
        0.0,-.3333333333333333,0.0,-.3333333333333333,.6666666666666666,0.0,
        -0.5,0.0,0.0,0.0,0.0,0.5}
};
static double prjperp[25][36] = {
    /* 1 */
    {0.5,-0.5,0.0,0.0,0.0,0.0,
        -0.5,0.5,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 2 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,-0.5,0.0,0.0,0.0,
        0.0,-0.5,0.5,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 3 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 4 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 5 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 6 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,-0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,-0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 7 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,-0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,-0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 8 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 9 */
    {0.5,0.0,0.0,0.0,-0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        -0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* A */
    {0.5,0.0,0.0,0.0,-0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        -0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* B */
    {0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* C */
    {0.5,0.0,0.0,0.0,0.0,-0.5,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        -0.5,0.0,0.0,0.0,0.0,0.5},
    /* D */
    {0.5,0.0,0.0,0.0,0.0,-0.5,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        -0.5,0.0,0.0,0.0,0.0,0.5},
    /* E */
    {0.5,0.0,0.0,0.0,0.0,0.5,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.5,0.0,0.0,0.0,0.0,0.5},
    /* F */
    {0.2,0.2,0.0,0.2,0.2,0.2,
        0.2,0.2,0.0,0.2,0.2,0.2,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.2,0.2,0.0,0.2,0.2,0.2,
        0.2,0.2,0.0,0.2,0.2,0.2,
        0.2,0.2,0.0,0.2,0.2,0.2},
    /* case678X g2=g4=g5=g6=0 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,1.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* case9ABCDEX g1=g4=g5=g6=0 */
    {1.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* FX g1=g2=g4=g5=g6=0 */
    {0.5,0.5,0.0,0.0,0.0,0.0,
        0.5,0.5,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,0.0,0.0,1.0},
    /* 67 */
    {0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,-0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,-0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.5,-0.5,
        0.0,0.0,0.0,0.0,-0.5,0.5},
    /* 9A */
    {0.5,0.0,0.0,0.0,-0.5,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.5,0.0,-0.5,
        -0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,0.0,0.0,-0.5,0.0,0.5},
    /* CD */
    {0.5,0.0,0.0,0.0,0.0,-0.5,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.5,-0.5,0.0,
        0.0,0.0,0.0,-0.5,0.5,0.0,
        -0.5,0.0,0.0,0.0,0.0,0.5},
    /* 12 g1=g2=g3 */
    {.6666666666666667,-.3333333333333333,-.3333333333333333,0.0, 0.0,0.0,
        -.3333333333333333,.6666666666666667,-.3333333333333333,0.0, 0.0,0.0,
        -.3333333333333333,-.3333333333333333,.6666666666666667,0.0, 0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0},
    /* 8F g4=-g2, g1+g2+g4+g5+g6 = 0 */
    {.3333333333333333,0.0,0.0,0.0,.3333333333333333,.3333333333333333,
        0.0,0.5,0.0,0.5,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.5,0.0,0.5,0.0,0.0,
        .3333333333333333,0.0,0.0,0.0,.3333333333333333,.3333333333333333,
        .3333333333333333,0.0,0.0,0.0,.3333333333333333,.3333333333333333},
    /* BF BF g5=-g1, g1+g2+g4+g5+g6 = 0 */
    {0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,.3333333333333333,0.0,.3333333333333333,0.0,.3333333333333333,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,.3333333333333333,0.0,.3333333333333333,0.0,.3333333333333333,
        0.5,0.0,0.0,0.0,0.5,0.0,
        0.0,.3333333333333333,0.0,.3333333333333333,0.0,.3333333333333333},
    /* EF g6=-g1.0, g1+g2+g4+g5+g6 = 0 */
    {0.5,0.0,0.0,0.0,0.0,0.5,
        0.0,.3333333333333334,0.0,.3333333333333333,.3333333333333333,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,.3333333333333333,0.0,.3333333333333334,.3333333333333333,0.0,
        0.0,.3333333333333333,0.0,.3333333333333333,.3333333333333333,0.0,
        0.5,0.0,0.0,0.0,0.0,0.5}
    
};

/* The following matrices are the transformation
 matrices that may be applied at the associated
 boundaries  */

static int MS[18][36] = {
    
    /* M_1 (g1 = g2, a -> b, b -> a) */
    {0,1,0,0,0,0,
        1,0,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,0,1,0,
        0,0,0,1,0,0,
        0,0,0,0,0,1 },
    
    /* M_2 (g2 = g3, b -> c, c -> b) */
    {1,0,0,0,0,0,
        0,0,1,0,0,0,
        0,1,0,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,0,1,
        0,0,0,0,1,0 },
    
    /* M_3 (g4 = 0, a -> -a) */
    {1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,1,0,0,
        0,0,0,0,-1,0,
        0,0,0,0,0,-1 },
    
    /* M_4 (g5 = 0, b -> -b) */
    {1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,-1,0,0,
        0,0,0,0,1,0,
        0,0,0,0,0,-1 },
    
    /* M_5 (g6 = 0, c -> -c) */
    {1,0,0,0,0,0,
        0,1,0,0,0,0,
        0,0,1,0,0,0,
        0,0,0,-1,0,0,
        0,0,0,0,-1,0,
        0,0,0,0,0,1 },
    
    /* M_6 (g2 = g4, g5 >= g6, b -> -b, c -> b - c) */
    {1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 1, 1,-1, 0, 0,
        0,-2, 0, 1, 0, 0,
        0, 0, 0, 0,-1, 1,
        0, 0, 0, 0, 0,-1 },
    
    /* M_7 (g2 = g4, g5 < g6, c -> b - c) */
    {1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 1, 1,-1, 0, 0,
        0, 2, 0,-1, 0, 0,
        0, 0, 0, 0,-1, 1,
        0, 0, 0, 0, 0, 1 },
    
    /* M_8 (g2 = -g4, a -> -a, c -> b + c) */
    {1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 1, 1, 1, 0, 0,
        0, 2, 0, 1, 0, 0,
        0, 0, 0, 0,-1,-1,
        0, 0, 0, 0, 0,-1 },
    
    /* M_9 (g1 = g5, g4 >= g6, b -> -b, c -> c - a) */
    {1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        1, 0, 1, 0,-1, 0,
        0, 0, 0,-1, 0, 1,
        -2, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0,-1 },
    
    /* M_A (g1 = g5, g4 < g6, c -> a - c) */
    {1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        1, 0, 1, 0,-1, 0,
        0, 0, 0,-1, 0, 1,
        2, 0, 0, 0,-1, 0,
        0, 0, 0, 0, 0, 1 },
    
    /* M_B (g1 = -g5, b -> -b, c -> a + c) */
    {1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        1, 0, 1, 0, 1, 0,
        0, 0, 0,-1, 0,-1,
        2, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0,-1 },
    
    /* M_C (g1 = g6, g4 >= g5, b -> -b, b -> b - a) */
    {1, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0,-1,
        0, 0, 1, 0, 0, 0,
        0, 0, 0,-1, 1, 0,
        0, 0, 0, 0,-1, 0,
        -2, 0, 0, 0, 0, 1 },
    
    /* M_D (g1 = g6, g4 < g5, b -> a - b) */
    {1, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0,-1,
        0, 0, 1, 0, 0, 0,
        0, 0, 0,-1, 1, 0,
        0, 0, 0, 0, 1, 0,
        2, 0, 0, 0, 0,-1 },
    
    /* M_E (g1 = -g6, b -> a + b, c -> -c ) */
    {1, 0, 0, 0, 0, 0,
        1, 1, 0, 0, 0, 1,
        0, 0, 1, 0, 0, 0,
        0, 0, 0,-1,-1, 0,
        0, 0, 0, 0,-1, 0,
        2, 0, 0, 0, 0, 1 },
    
    /* M_F (g1+g2+g3+g4+g5+g6 = g3, c -> -(a+b+c)) */
    {1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        1, 1, 1, 1, 1, 1,
        0,-2, 0,-1, 0,-1,
        -2, 0, 0, 0,-1,-1,
        0, 0, 0, 0, 0, 1 },
    
    /* M_1.M_2 (a -> b, b -> c, c -> a) */
    {0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 1, 0 },
    
    /* M_2.M_1 (a -> c, b -> a, c -> b) */
    {0, 1, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0,
        1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 1,
        0, 0, 0, 1, 0, 0 },
    
    /* M_2.M_1.M_2 (a -> c, c -> a) */
    {0, 0, 1, 0, 0, 0,
        0, 1, 0, 0, 0, 0,
        1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1,
        0, 0, 0, 0, 1, 0,
        0, 0, 0, 1, 0, 0}

};

static void cpyvn(int n, double src[], double dst[] ) {
    int i;
    for (i = 0; i < n; i++) {
        dst[i] = src[i];
    }
}


static void imv6 (double v1[6], int m[36], double v2[6]) {
    int i, j;
    double sum;
    for (i = 0; i < 6; i++) {
        sum = 0.0;
        for(j=0; j < 6; j++) {
            sum = sum + ((double)m[6*i+j])*v1[j];
        }
        v2[i] = sum;
    }
}

static void rmv6 (double v1[6], double m[36], double v2[6]) {
    int i, j;
    double sum;
    for (i = 0; i < 6; i++) {
        sum = 0.0;
        for(j=0; j < 6; j++) {
            sum = sum + m[6*i+j]*v1[j];
        }
        v2[i] = sum;
    }
}


/*     Map a G6 vector onto the boundaries after
       applying the 24-way unfolding */

static void bdmaps(double gvec[6],
            double vecs[24][6],
            double dists[15][24],
            double pgs[15][24][6],
            double mpgs[15][24][6]) {
    
    int ii, jj;
    double xtemp;
    
    
    /*
     0 --  5 +++
     6 -- 11 +--
     12 -- 17 -+-
     18 -- 23 --+
     */
    
    for (ii = 0; ii < 24; ii++) {
        cpyvn(6,gvec,vecs[ii]);
        if (ii >= 6 && ii <= 17) {
            vecs[ii][5] = -vecs[ii][5];
        }
        if ((ii >= 6 && ii <= 11)|| ii >= 18) {
            vecs[ii][4] = -vecs[ii][4];
        }
        if (ii>=12) {
            vecs[ii][3] = -vecs[ii][3];
        }
        jj = ii%6;
        if (jj==2||jj==3) {
            xtemp = vecs[ii][0];
            vecs[ii][0] = vecs[ii][1];
            vecs[ii][1] = xtemp;
            xtemp = vecs[ii][3];
            vecs[ii][3] = vecs[ii][4];
            vecs[ii][4] = xtemp;
        }
        if (jj==1||jj>=4) {
            xtemp = vecs[ii][1];
            vecs[ii][1] = vecs[ii][2];
            vecs[ii][2] = xtemp;
            xtemp = vecs[ii][4];
            vecs[ii][4] = vecs[ii][5];
            vecs[ii][5] = xtemp;
        }
        if (jj==4) {
            xtemp = vecs[ii][0];
            vecs[ii][0] = vecs[ii][1];
            vecs[ii][1] = xtemp;
            xtemp = vecs[ii][3];
            vecs[ii][3] = vecs[ii][4];
            vecs[ii][4] = xtemp;
        }
        if (jj==3||jj==5) {
            xtemp = vecs[ii][1];
            vecs[ii][1] = vecs[ii][2];
            vecs[ii][2] = xtemp;
            xtemp = vecs[ii][4];
            vecs[ii][4] = vecs[ii][5];
            vecs[ii][5] = xtemp;
        }
        dists[0][ii] = fabs(vecs[ii][1]-vecs[ii][0])/sqrt(2.);
        dists[1][ii] = fabs(vecs[ii][2]-vecs[ii][1])/sqrt(2.);
        dists[2][ii] = fabs(vecs[ii][3]);
        dists[3][ii] = fabs(vecs[ii][4]);
        dists[4][ii] = fabs(vecs[ii][5]);
        dists[5][ii] = fabs(vecs[ii][1]-vecs[ii][3])/sqrt(2.);
        dists[6][ii] = fabs(vecs[ii][1]-vecs[ii][3])/sqrt(2.);
        dists[7][ii] = fabs(vecs[ii][1]+vecs[ii][3])/sqrt(2.);
        dists[8][ii] = fabs(vecs[ii][0]-vecs[ii][4])/sqrt(2.);
        dists[9][ii] = fabs(vecs[ii][0]-vecs[ii][4])/sqrt(2.);
        dists[10][ii] = fabs(vecs[ii][0]+vecs[ii][4])/sqrt(2.);
        dists[11][ii] = fabs(vecs[ii][0]-vecs[ii][5])/sqrt(2.);
        dists[12][ii] = fabs(vecs[ii][0]-vecs[ii][5])/sqrt(2.);
        dists[13][ii] = fabs(vecs[ii][0]+vecs[ii][5])/sqrt(2.);
        dists[14][ii] = fabs(vecs[ii][0]+vecs[ii][1]+vecs[ii][3]+vecs[ii][4]+vecs[ii][5])/sqrt(5.);
        for (jj = 0; jj < 15; jj++ ) {
            rmv6(vecs[ii], prj[jj], pgs[jj][ii]);
            imv6(pgs[jj][ii], MS[jj], mpgs[jj][ii]);
        }
    }
}

/*    
    Map a G6 vector onto the intersection of the
    face diagonal and body diagonal boundaries.

    For a given g6 point, we need the distance
    to  8F directly or via 69 and M_6
        BF directly of via 69 and M_9
        EF firectly or via 6C and M_C
    in addition we need M_F applied to the
    F targets

    So the distances and pgs in order are
    ||P_8F_perp.g||  P_8F.g  M_8.P_8F.g
    ||P_BF_perp.g||  P_BF.g  M_B.P_BF.g
    ||P_EF_perp.g||  P_EF.g  M_E.P_EF.g
    ||P_6C_perp.g||  P_6C.g  M_C.P_6C.g
    ||P_69_perp.g||  P_69.g  M_6.P_69.g
                             M_9.P_69.g
    
     bfmaps must be called first to generate the vecs
     array.

     The viable pairs for a given g1 and g2 for
     consideration are {

     boundary  distances
     69:  ||P_69_perp.g1|| ||P_69_perp.g2|| P_69.g1 -- P_69.g2
          ||P_69_perp.g1|| ||P_8F_perp.g2|| P_69.g1 -- M_8.P_8F.g2
          ||P_69_perp.g1|| ||P_BF_perp.g2|| P_69.g1 -- M_B.P_BF.g2
          and with g1 and g2 exchanged
     6C:  ||P_6C_perp.g1|| ||P_6C_perp.g2|| P_6C.g1 -- P_6C.g2
          ||P_6C_perp.g1|| ||P_EF_perp.g2|| P_69.g1 -- M_C.P_EF.g2
          and with g1 and g2 exchanged
     8F:  ||P_8F_perp.g1|| ||P_8F_perp.g2|| P_8F.g1 -- P_8F.g2
          ||P_8F_perp.g1|| ||P_69_perp.g2|| P_8F.g1 -- M_6.P_69.g
     BF:  ||P_BF_perp.g1|| ||P_BF_perp.g2|| P_BF.g1 -- P_BF.g2
          ||P_BF_perp.g1|| ||P_69_perp.g2|| P_BF.g1 -- M_9.P_69.g
     EF:  ||P_EF_perp.g1|| ||P_EF_perp.g2|| P_EF.g1 -- P_EF.g2
          ||P_EF_perp.g1|| ||P_6C_perp.g2|| P_EF.g1 -- M_C.P_6C.g
 */

static void bdfmaps(double vecs[24][6],
             double dists[5][24],
             double pgs[5][24][6],
             double mpgs[6][24][6],
             int nmpgs[5]) {
      
      int ii, jj;
      
      double pgtemp[6];
      
      for (ii=0; ii < 24; ii++) {
          dists[3][ii] = sqrt((vecs[ii][0]-vecs[ii][5])*(vecs[ii][0]-vecs[ii][5])
                                   + (vecs[ii][1]-vecs[ii][3])*(vecs[ii][1]-vecs[ii][3]))/sqrt(2.);
          dists[4][ii] = sqrt((vecs[ii][0]-vecs[ii][4])*(vecs[ii][0]-vecs[ii][4])
                                   + (vecs[ii][1]-vecs[ii][3])*(vecs[ii][1]-vecs[ii][3]))/sqrt(2.);
          dists[0][ii] = sqrt( 2.*(vecs[ii][5]+vecs[ii][4]+vecs[ii][0])
                                   *(vecs[ii][5]+vecs[ii][4]+vecs[ii][0])
                                   +3.*(vecs[ii][3]+vecs[ii][1])*(vecs[ii][3]+vecs[ii][1]))/sqrt(6.);
          dists[1][ii] = sqrt( 2.*(vecs[ii][5]+vecs[ii][3]+vecs[ii][1])*(vecs[ii][5]+vecs[ii][3]+vecs[ii][1])
                                   +3.*(vecs[ii][4]+vecs[ii][0])*(vecs[ii][4]+vecs[ii][0]))/sqrt(6.);
          dists[2][ii] = sqrt( 2.*(vecs[ii][4]+vecs[ii][3]+vecs[ii][1])*(vecs[ii][4]+vecs[ii][3]+vecs[ii][1])
                                   +3.*(vecs[ii][5]+vecs[ii][0])*(vecs[ii][5]+vecs[ii][0]))/sqrt(6.);
          /*
           prj[22] is P_8F
           prj[23] is P_BF
           prj[24] is P_EF
           */
          for (jj = 0; jj < 3; jj++) {
              rmv6(vecs[ii], prj[22+jj], pgs[jj][ii]);
              imv6(pgs[jj][ii],MS[jj*3+7],mpgs[jj][ii]);
          }
          rmv6(vecs[ii],prj[5],pgtemp);
          rmv6(pgtemp,prj[11],pgs[3][ii]);
          rmv6(pgtemp,prj[8],pgs[4][ii]);
          
          imv6(pgs[3][ii],MS[11],mpgs[3][ii]);
          imv6(pgs[4][ii],MS[5],mpgs[4][ii]);
          imv6(pgs[4][ii],MS[8],mpgs[5][ii]);
      }
      nmpgs[0] = 1;
      nmpgs[1] = 1;
      nmpgs[2] = 1;
      nmpgs[3] = 1;
      nmpgs[4] = 2;
  }



/*     Compute the best distance between 2 G6 vectors
     allowing for cell-preserving sign changes in
     g4,5,6
*/
static double g456distsq(double v1[6], double v2[6]){
    
    double vtemp;
    double xdot;
    int ii;
    double dist;
    
    xdot = 0.;
    
    for (ii = 0; ii < 6; ii++ ) {
        vtemp = v1[ii]-v2[ii];
        xdot = xdot+vtemp*vtemp;
    }
    dist = (xdot+
                4.*NCD_min(NCD_min(NCD_min(0.,
                                  v1[3]*v2[3]+v1[4]*v2[4]),
                             v1[3]*v2[3]+v1[5]*v2[5]),
                        v1[4]*v2[4]+v1[5]*v2[5]));
    return dist;
}
static double g456dist(double v1[6], double v2[6]){

    return sqrt(g456distsq(v1,v2));
    
}

/*   Macro version of g456dist
 Compute the best distance between 2 G6 vectors
 allowing for cell-preserving sign changes in
 g4,5,6
 */

#define CNCM_g456distsq(v1,v2) \
    ( \
      (v1[0]-v2[0])*(v1[0]-v2[0])+\
      (v1[1]-v2[1])*(v1[1]-v2[1])+\
      (v1[2]-v2[2])*(v1[2]-v2[2])+\
      (v1[3]-v2[3])*(v1[3]-v2[3])+\
      (v1[4]-v2[4])*(v1[4]-v2[4])+\
      (v1[5]-v2[5])*(v1[5]-v2[5])+\
      4.*NCD_min(NCD_min(NCD_min(0.,       \
             v1[3]*v2[3]+v1[4]*v2[4]), \
             v1[3]*v2[3]+v1[5]*v2[5]), \
             v1[4]*v2[4]+v1[5]*v2[5]))

#define CNCM_g456dist(v1,v2) \
    sqrt(CNCM_g456distsq(v1,v2))


/*     Compute the best distance between 2 G6 vectors
     allowing for permulations of g1, g2, g3 as
     well as sign changes
 */


static double g123distsq(double v1[6],double v2[6]) {
    double vtemp[6];
    double distsq;
    int i;
    for (i = 0; i < 6; i++ ) {
        vtemp[i] = v2[i];
    }
    /*     123 */
    distsq = CNCM_g456distsq(v1,vtemp);
    /*     213 */
    vtemp[0]=v2[1];
    vtemp[1]=v2[0];
    vtemp[3]=v2[4];
    vtemp[4]=v2[3];
    distsq = NCD_min(distsq,CNCM_g456distsq(v1,vtemp));
    /*     231 */
    vtemp[1]=v2[2];
    vtemp[2]=v2[0];
    vtemp[4]=v2[5];
    vtemp[5]=v2[3];
    distsq = NCD_min(distsq,CNCM_g456distsq(v1,vtemp));
    /*     321 */
    vtemp[0]=v2[2];
    vtemp[1]=v2[1];
    vtemp[3]=v2[5];
    vtemp[4]=v2[4];
    distsq = NCD_min(distsq,CNCM_g456distsq(v1,vtemp));
    /*     312 */
    vtemp[2]=v2[1];
    vtemp[1]=v2[0];
    vtemp[5]=v2[4];
    vtemp[4]=v2[3];
    distsq = NCD_min(distsq,CNCM_g456distsq(v1,vtemp));
    /*     132 */
    vtemp[0]=v2[0];
    vtemp[1]=v2[2];
    vtemp[3]=v2[3];
    vtemp[4]=v2[5];
    distsq = NCD_min(distsq,CNCM_g456distsq(v1,vtemp));
    return distsq;
}


/*
     Compute the NCD_minimal distance between two Niggli-reduced
     vectors in the Niggli Cone following the embedding paths
     to the 15 boundaries
 */

double NCDist(double gvec1[6],double gvec2[6]) {
    double vecs1[24][6], dists1[15][24];
    double pgs1[15][24][6], mpgs1[15][24][6];
    double fdists1[5][24],fdists2[5][24];
    double fpgs1[5][24][6],fpgs2[5][24][6];
    double fmpgs1[6][24][6], fmpgs2[6][24][6];
    int nmpgs[5];
    double vecs2[24][6], dists2[15][24];
    double pgs2[15][24][6], mpgs2[15][24][6];
    double dpg1pg2;
    double distsq;
    int jx1, jx2, ix2;
    int jord[15] = {0,1,8,9,10,5,6,7,11,12,13,2,3,4,14},
    jord2[15] = {0,1,10,9,8,7,6,5,13,12,11,2,3,4,14};
    int i1,i2,j1,j2;
    
    bdmaps(gvec1,vecs1,dists1,pgs1,mpgs1);
    bdfmaps(vecs1,fdists1,fpgs1,fmpgs1,nmpgs);
    bdmaps(gvec2,vecs2,dists2,pgs2,mpgs2);
    bdfmaps(vecs2,fdists2,fpgs2,fmpgs2,nmpgs);
    
    distsq = g123distsq(gvec1,gvec2);
    
    for (i1=0; i1<24; i1++) {
        for (ix2 = 1; ix2 < 24; ix2++) {
            i2 = (i1+ix2)%24;
            for (jx1 = 0; jx1 < 15; jx1++) {
              j1 = jord[jx1];
              if (dists1[j1][i1]*dists1[j1][i1] < distsq) {
                jx2 = jx1;
                j2 = jord2[jx2];
                if (j1==j2) {
                    if((dists1[j1][i1]+dists2[j2][i2])*(dists1[j1][i1]+dists2[j2][i2]) < distsq) {
                        dpg1pg2 = NCD_min(NCD_min(NCD_min(CNCM_g456distsq(pgs1[j1][i1],pgs2[j2][i2]),
                                              CNCM_g456distsq(pgs1[j1][i1],mpgs2[j2][i2])),
                                          CNCM_g456distsq(mpgs1[j1][i1],pgs2[j2][i2])),
                                      CNCM_g456distsq(mpgs1[j1][i1],mpgs2[j2][i2]));
                        distsq = NCD_min(distsq,
                                   ((dists1[j1][i1]+dists2[j2][i2])
                                             *(dists1[j1][i1]+dists2[j2][i2])
                                             + dpg1pg2));
                    }
                    else
                        if ((dists1[j1][i1]+dists2[j1][i2])*(dists1[j1][i1]+dists2[j1][i2]) < distsq) {
                            dpg1pg2 =  CNCM_g456distsq(pgs1[j1][i1],pgs2[j1][i2]);
                            distsq = NCD_min(distsq,
                                       ((dists1[j1][i1]+dists2[j1][i2])
                                                 *(dists1[j1][i1]+dists2[j1][i2])
                                                 + dpg1pg2));
                        }
                    if ((dists1[j1][i1]+dists2[j2][i2])*(dists1[j1][i1]+dists2[j2][i2]) < distsq) {
                        dpg1pg2 =  CNCM_g456distsq(pgs1[j1][i1],mpgs2[j2][i2]);
                        distsq = NCD_min(distsq, ((dists1[j1][i1]+dists2[j2][i2])
                                               *(dists1[j1][i1]+dists2[j2][i2])
                                               + dpg1pg2));
                    }
                }
              }
            }
        }
        /*     69 */
        if ((fdists1[4][i1]+fdists2[4][i2])*(fdists1[4][i1]+fdists2[4][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[4][i1]+fdists2[4][i2])*(fdists1[4][i1]+fdists2[4][i2])
                                 + CNCM_g456distsq(fpgs1[4][i1],fpgs2[4][i2])));
        }
        if ((fdists1[4][i1]+fdists2[0][i2])*(fdists1[4][i1]+fdists2[0][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[4][i1]+fdists2[0][i2])
                                 *(fdists1[4][i1]+fdists2[0][i2])
                                 + CNCM_g456distsq(fpgs1[4][i1],fmpgs2[0][i2])));
        }
        if ((fdists1[4][i1]+fdists2[1][i2])*(fdists1[4][i1]+fdists2[1][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[4][i1]+fdists2[1][i2])
                                 *(fdists1[4][i1]+fdists2[1][i2])
                                 + CNCM_g456distsq(fpgs1[4][i1],fmpgs2[1][i2])));
        }
        if ((fdists1[0][i1]+fdists2[4][i2])*(fdists1[0][i1]+fdists2[4][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[0][i1]+fdists2[4][i2])
                                 *(fdists1[0][i1]+fdists2[4][i2])
                                 + CNCM_g456distsq(fmpgs1[0][i1],fpgs2[4][i2])));
        }
        if ((fdists1[1][i1]+fdists2[4][i2])*(fdists1[1][i1]+fdists2[4][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[1][i1]+fdists2[4][i2])
                                 *(fdists1[1][i1]+fdists2[4][i2])
                                 + CNCM_g456distsq(fmpgs1[1][i1],fpgs2[4][i2])));
        }
        /*     6C */
        if ((fdists1[3][i1]+fdists2[3][i2])*(fdists1[3][i1]+fdists2[3][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[3][i1]+fdists2[3][i2])
                                 *(fdists1[3][i1]+fdists2[3][i2])
                                 + CNCM_g456distsq(fpgs1[3][i1],fpgs2[3][i2])));
        }
        if ((fdists1[3][i1]+fdists2[2][i2])*(fdists1[3][i1]+fdists2[2][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[3][i1]+fdists2[2][i2])
                                 *(fdists1[3][i1]+fdists2[2][i2])
                                 + CNCM_g456distsq(fpgs1[3][i1],fmpgs2[2][i1])));
        }
        if ((fdists1[2][i1]+fdists2[3][i2])*(fdists1[2][i1]+fdists2[3][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[2][i1]+fdists2[3][i2])
                                 *(fdists1[2][i1]+fdists2[3][i2])
                                 + CNCM_g456distsq(fmpgs1[2][i1],fpgs2[3][i2])));
        }
        /*     8F */
        if ((fdists1[0][i1]+fdists2[0][i2])*(fdists1[0][i1]+fdists2[0][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[0][i1]+fdists2[0][i2])
                                 *(fdists1[0][i1]+fdists2[0][i2])
                                 + CNCM_g456distsq(fpgs1[0][i1],fpgs2[0][i2])));
        }
        if ((fdists1[0][i1]+fdists2[4][i2])*(fdists1[0][i1]+fdists2[4][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[0][i1]+fdists2[4][i2])
                                 *(fdists1[0][i1]+fdists2[4][i2])
                                 + CNCM_g456distsq(fpgs1[0][i1],fmpgs2[4][i2])));
        }
        if ((fdists1[4][i1]+fdists2[0][i2])*(fdists1[4][i1]+fdists2[0][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[4][i1]+fdists2[0][i2])
                                 *(fdists1[4][i1]+fdists2[0][i2])
                                 + CNCM_g456distsq(fmpgs1[4][i1],fpgs2[0][i2])));
        }
        /*     BF */
        if ((fdists1[1][i1]+fdists2[1][i2])*(fdists1[1][i1]+fdists2[1][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[0][i1]+fdists2[0][i2])
                                 *(fdists1[0][i1]+fdists2[0][i2])
                                 + CNCM_g456distsq(fpgs1[0][i1],fpgs2[0][i2])));
        }
        if ((fdists1[1][i1]+fdists2[4][i2])*(fdists1[1][i1]+fdists2[4][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[1][i1]+fdists2[4][i2])
                                 *(fdists1[1][i1]+fdists2[4][i2])
                                 + CNCM_g456distsq(fpgs1[1][i1],fmpgs2[5][i2])));
        }
        if ((fdists1[4][i1]+fdists2[1][i2])*(fdists1[4][i1]+fdists2[1][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[4][i1]+fdists2[1][i2])
                                 *(fdists1[4][i1]+fdists2[1][i2])
                                 + CNCM_g456distsq(fmpgs1[5][i1],fpgs2[1][i2])));
        }
        /*     EF */
        if ((fdists1[2][i1]+fdists2[2][i2])*(fdists1[2][i1]+fdists2[2][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[2][i1]+fdists2[2][i2])
                                 *(fdists1[2][i1]+fdists2[2][i2])
                                 + CNCM_g456distsq(fpgs1[2][i1],fpgs2[2][i1])));
        }
        if ((fdists1[2][i1]+fdists2[3][i2])*(fdists1[2][i1]+fdists2[3][i2])< distsq) {
            distsq = NCD_min( distsq,
                       ((fdists1[2][i1]+fdists2[3][i2])
                                 *(fdists1[2][i1]+fdists2[3][i2])
                                 + CNCM_g456distsq(fpgs1[2][i1],fmpgs2[3][i2])));
        }
        if ((fdists1[3][i1]+fdists2[2][i2])*(fdists1[3][i1]+fdists2[2][i2])< distsq) {
            distsq = NCD_min( distsq,((fdists1[3][i1]+fdists2[2][i2])
                                       *(fdists1[3][i1]+fdists2[2][i2])
                                       + CNCM_g456distsq(fmpgs1[3][i1],fpgs2[2][i1])));
        }
    }
    return sqrt(distsq);
}
