
#define ARMA_DONT_USE_BLAS
#define ARMA_DONT_USE_LAPACK
#include <armadillo>

#include "BasicDistance.h"
#include "Cell.h"
#include "MakeGaol.h"
#include "Reducer.h"

#include <cfloat>
#include <string>


/*  class Reducer
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
A class to implement Niggli reduction in the G6 space. Reduce returns a reduced cell
   and the matrix that converts the input cell to the reduced cell.

   Reduced(void)
                 == constructor -- nothing to do here
   static void Reduce( const arma::vec6& vi, arma::mat66& m, arma::vec6& vout, const double delta )
                 == returns vout as the reduced cell of vi and m, the conversion matrix
   static bool NearRed( const arma::vec6& gvec, const double delta );
                 == determines whether gvec is reduced within delta
   static bool Near2Red( const arma::vec6& gvec, const double delta, arma::vec6& vout, double& dist )
                 == determines whether gvec is reduced within delta, and returns the reduced cell
                    and how far from reduced
   static void Reporter( const std::string& text, const arma::vec6& vin, const arma::vec6& vout, const arma::mat66& m )
                 == prints information about each step in reduction (including standard presentation)

private:
   static void MKnorm( const arma::vec6& vi, arma::mat66& m, arma::vec6& vout, const double delta  )
                 == internal function to convert vi to standard presentation - the matrix
                    to implement that change is returned
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*/


double previousVolume = DBL_MAX;
const bool DEBUG_REDUCER(false);

//-----------------------------------------------------------------------------
// Name: g6sign()
// Description: returns the value of the first variable with the sign of the second
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
double g6sign( const double d, const double x )
{
   return( x>0.0 ? d : -d );
}

//-----------------------------------------------------------------------------
// Name: Reporter()
// Description: prints out the intermediate step values and looks for
//              changes in the volume (indication of error)
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void Reducer::Reporter( const std::string& text, const arma::vec6& vin, const arma::vec6& vout, const arma::mat66& m )
{
   if ( ! DEBUG_REDUCER ) return;

   const double volume = Cell(vout).Volume( );
   const  double previousVolume = Cell(vin).Volume( );
   const bool volumeChanged = ( fabs(volume-previousVolume) / std::max(volume, previousVolume) ) > 1.0E-12;
   if ( volumeChanged )
   {
      printf( "**************************************************************************************\n" );
   }

   if ( text.empty( ) )
   {
      const int i19191 = 19191;
   }

   MakeGaol::PrintG6( std::string("vin ")+text.c_str(), vin );
   printf("\n");
   MakeGaol::PrintG6( "vout ", vout );

   if ( volumeChanged )
   {
      printf( "\nVolume(vout) %f", volume );
      printf( "    PREVIOUS VOLUME %f\n", previousVolume );
   }
   else
   {
      printf( "\nVolume %f\n", volume );
   }

   printf("\n");
   MakeGaol::PrintM66_INT( "m ", m );
   printf("\n");
} // end Reporter

//-----------------------------------------------------------------------------
// Name: Reducer()
// Description: constructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
Reducer::Reducer(void)
{
}

//-----------------------------------------------------------------------------
// Name: Reducer()
// Description: destructor
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
Reducer::~Reducer(void)
{
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
// These are the matrices that can be used to convert a vector to standard
// presentation. They are used in MKnorm. There are NO OTHER POSSIBILITIES.
// The numbering is from Andrews and Bernstein, 1988
   const static arma::mat66 spnull( "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1" );

   const static arma::mat66 sp1( "0 1 0 0 0 0; 1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 1 0 0; 0 0 0 0 0 1" );
   const static arma::mat66 sp2( "1 0 0 0 0 0; 0 0 1 0 0 0; 0 1 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 0 1; 0 0 0 0 1 0" );


   const static arma::mat66 sp34a( "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0 -1 0; 0 0 0 0 0  1" );
   const static arma::mat66 sp34b( "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 -1 0 0; 0 0 0 0  1 0; 0 0 0 0 0 -1" );
   const static arma::mat66 sp34c( "1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0  1 0 0; 0 0 0 0 -1 0; 0 0 0 0 0 -1" );
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


//-----------------------------------------------------------------------------
// Name: MKnorm()
// Description: changes a G6 vector to standard presentation (often called
//              normalization in the literature) and returns the standard
//              vector and the matrix that changes the input vector to the
//              standard one
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void Reducer::MKnorm( const arma::vec6& vi, arma::mat66& m, arma::vec6& vout, const double delta  )
{
   bool again = true;
   int ncycle = 0;
   arma::mat66 mat( spnull );
   arma::vec6 vin;

   vin = vi;

   // assure that g1<=g2<=g3
   while ( again && (ncycle < 5) )
   {
      ++ncycle;
      again = false;

      std::string sptext;
      if ( (fabs(vin[0]) > fabs(vin[1])+delta) ||
           (vin[0] ==vin[1] && delta<1.0E-6 && fabs(vin[3])>fabs(vin[4])+delta ) )
      { // SP1
         mat    = sp1;
         again  = true;
         sptext = "SP1";
      }
      else if ( (fabs(vin[1]) > fabs(vin[2])+delta) ||
           (vin[1] ==vin[2] && delta<1.0E-6 && fabs(vin[4])>fabs(vin[5])+delta ) )
      { // SP2
         mat    = sp2;
         again  = true;
         sptext = "SP2";
      }

      if ( again )
      {
         // Accumulate the total transformation from the input vector into vout and
         // the total transformation itself into matrix m.
         const arma::mat66 mtemp = mat*m;
         m = mtemp;
         vout = mat*vin;
         Reporter( sptext, vin, vout, mat );
         vin = vout;
      }
   }

   // now we assure (within delta) that the vector is +++ or ---

   // HJB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // DELTA IS NO LONGER USED HERE -- what to do?
   int bMinusPattern = 0;
   if( vin[3] < 1.0E-10 ) bMinusPattern += 4;
   if( vin[4] < 1.0E-10 ) bMinusPattern += 2;
   if( vin[5] < 1.0E-10 ) bMinusPattern += 1;
   std::string sptext2( "ERROR" );;

   switch( bMinusPattern )
   {
   case 0:
      {
         mat = spnull;
         sptext2 = "no mknorm action sp1,sp2-0";
         break;
      }
   case 1:
      {
         mat = sp34a;
         sptext2 = "SP34a-1";
         break;
      }
   case 2:
      {
         mat = sp34b;
         sptext2 = "SP34b-2";
         break;
      }
   case 3:
      {
         mat = sp34c;
         sptext2 = "SP34c-3";
         break;
      }
   case 4:
      {
         mat = sp34c;
         sptext2 = "SP34c-4";
         break;
      }
   case 5:
      {
         mat = sp34b;
         sptext2 = "SP34b-5";
         break;
      }
   case 6:
      {
         mat = sp34a;
         sptext2 = "SP34a-6";
         break;;
      }
   case 7:
      {
         mat = spnull;
         sptext2 = "no mknorm action sp1,sp2-7";
         break;
      }
   }

   // Accumulate the total transformation from the input vector into vout and
   // the total transformation itself into matrix m.
   const arma::mat66 mtemp = mat*m;
   m = mtemp;
   vout = mat*vin;
   Reporter( sptext2, vin, vout, mat );
}  // end MKnorm

//-----------------------------------------------------------------------------
// Name: Reduce()
// Description: performs Niggli reduction and returns the standard reduced
//              vector and the matrix that changes the input vector to the
//              reduced one
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool Reducer::Reduce( const arma::vec6& vi, arma::mat66& m, arma::vec6& vout, const double delta )
{
   arma::vec6 vin;
   arma::mat66 m1;
   int ncycle = 0;
   bool again = true;
   const bool debug = true;

   m1.eye( );
   vin = vi;

   /* Mapping from Fortran indices:

   1  2  3  4  5  6
   0  6 12 18 24 30

   7  8  9 10 11 12
   1  7 13 19 25 31

   13 14 15 16 17 18
   2  8 14 20 26 32

   19 20 21 22 23 24
   3  9 15 21 27 33

   25 26 27 28 29 30
   4 10 16 22 28 34

   31 32 33 34 35 36
   5 11 17 23 29 35

   */

   // Try some number of times to reduce the input vector and accumulate
   // the changing vector and total transformation matrix
   // The limit on the number of cycles is because (whether because of
   // floating point rounding or algorithm issues) there might be a
   // case of an infinite loop. The designations R5-R12 are from
   // Andrews and Bernstein, 1988
   while ( again && ncycle < 26 )
   {
      m1.eye( );
      MKnorm( vin, m1, vout, delta );
      vin = vout;
      const arma::mat66 m2temp = m1*m;
      m = m2temp;
      m1.eye( );

      if ( fabs(vin[3]) > fabs(vin[1])+delta )
      { // R5
         m1[8]  = 1;
         m1[20] = -g6sign( 1.0, vin[3] );
         m1[9]  = 2.0*m1[20];
         m1[34] = m1[20];
         again  = true;
         const arma::mat66 m2 = m1*m;
         m = m2;
         vout = m1*vin;
         Reporter( "R5", vin, vout, m1 );
      }
      else if ( fabs(vin[4]) > fabs(vin[0])+delta )
      { // R6
         m1[2]  = 1.0;
         m1[26] = -g6sign(1.0, vin[4] );
         m1[33] = m1[26];
         m1[4]  = 2.0*m1[26];
         again  = true;
         const arma::mat66 m2 = m1*m;
         m = m2;
         vout = m1*vin;
         Reporter( "R6", vin, vout, m1 );
      }
      else if ( fabs(vin[5]) > fabs(vin[0])+delta )
      { // R7
         m1[1]  = 1.0;
         m1[31] = -g6sign( 1.0, vin[5] );
         m1[27] = m1[31];
         m1[5]  = 2.0*m1[31];
         again  = true;
         const arma::mat66 m2 = m1*m;
         m = m2;
         vout = m1*vin;
         Reporter( "R7", vin, vout, m1 );
      }
      else if ( vin[3]+vin[4]+vin[5]+fabs(vin[0])+fabs(vin[1])+delta < 0.0 )
      { //R8
         m1[ 2] = 1.0;
         m1[ 8] = 1.0;
         m1[14] = 1.0;
         m1[20] = 1.0;
         m1[26] = 1.0;
         m1[32] = 1.0;
         m1[ 9] = -2.0;
         m1[21] = -1.0;
         m1[33] = -1.0;
         m1[ 4] = -2.0;
         m1[28] = -1.0;
         m1[34] = -1.0;
         again = true;
         const arma::mat66 m2 = m1*m;
         m = m2;
         vout = m1*vin;
         Reporter( "R8", vin, vout, m1 );
      }
      else if ( (fabs(vin[3]-vin[1])<=delta && 2.0*vin[4]-delta<vin[5] ) ||
         (fabs(vin[3]+vin[1])<=delta &&
         vin[5]-delta<=0.0 ) )
      { // R9  There is an error in the paper says "2g5<g5" should be "2g5<g6"
         m1[ 8] = 1.0;
         m1[ 9] = 2.0;
         m1[28] = -1.0;

         if ( vin[3] <= delta )
         {
            m1[20] = 1.0;
            m1[34] = -1.0;
            m1[35] = -1.0;
         }
         else
         {
            m1[20] = -1.0;
            m1[21] = -1.0;
            m1[34] = 1.0;
         }
         again = true;
         const arma::mat66 m2 = m1*m;
         m = m2;
         vout = m1*vin;
         Reporter( "R9", vin, vout, m1 );
      }
      else if ( (fabs(vin[4]-vin[0])<=delta && 2.0*vin[3]-delta<vin[5] ) ||
         (fabs(vin[4]+vin[0])<=delta &&
         vin[5]-delta<=0.0 ) )
      { //R10
         m1[ 2] = 1.0;
         m1[21] = -1.0;
         m1[ 4] = 2.0;
         if( vin[4] <= delta )
         {
            m1[26] = 1.0;
            m1[33] = -1.0;
            m1[35] = -1.0;
         }
         else
         {
            m1[26] = -1.0;
            m1[28] = -1.0;
            m1[33] = 1.0;
         }
         again = true;
         const arma::mat66 m2 = m1*m;
         m = m2;
         vout = m1*vin;
         Reporter( "R10", vin, vout, m1 );
      }
      else if ( (fabs(vin[5]-vin[0])<=delta && 2.0*vin[3]-delta<vin[4] ) ||
         (fabs(vin[5]+vin[0])<=delta &&
         vin[4]-delta<=0.0 ) ) // paper says g6<0, but it seems wrong
      { // R11
         m1[ 1] = 1.0;
         m1[ 5] = 2.0;
         m1[21] = -1.0;
         if( vin[5] <= delta )
         {
            m1[31] = 1.0;
            m1[27] = -1.0;
            m1[28] = -1.0;
         }
         else
         {
            m1[31] = -1.0;
            m1[27] = 1.0;
            m1[35] = -1.0;
         }
         again = true;
         const arma::mat66 m2 = m1*m;
         m = m2;
         vout = m1*vin;
         Reporter( "R11", vin, vout, m1 );
      }
      else if ( fabs(vin[3]+vin[4]+vin[5]+fabs(vin[0])+fabs(vin[1])) <= delta && (2.0*(fabs(vin[0])+vin[4])+vin[5]>delta) )
      { //R12 (same as R8)
         m1[ 2] = 1.0;
         m1[ 8] = 1.0;
         m1[14] = 1.0;
         m1[20] = 1.0;
         m1[26] = 1.0;
         m1[32] = 1.0;
         m1[ 9] = 2.0;
         m1[33] = 1.0;
         m1[ 4] = 2.0;
         m1[34] = 1.0;
         again = true;
         //m1 = arma::mat66( "1 0 0 0 0 0; 0 1 0 0 0 0; 1 1 1 1 1 1; 0 2 0 1 0 1; 2 0 0 0 1 1; 0 0 0 0 0 1");
         const arma::mat66 m2 = m1*m;
         m = m2;
         vout = m1*vin;
         Reporter( "R12", vin, vout, m1 );
      }
      else
      {
         again = false;
         vout = vin;
      }

      // probably don't need to do this group of code when again==false !!!!!!!!!!!!!!
      m1.eye( );
      MKnorm( vout, m1, vin, delta );
      const arma::mat66 mtemp = m1*m;
      m = mtemp;
      Reporter( "vout after MKnorm at end of reduced cycle", vout, vin, m1 );
      vout = vin;

      if ( vin[0] < 0.0 || vin[1] < 0.0 || vin[2] < 0.0 )
      {
         // ERROR ERROR ERROR
         fprintf(stderr," Negative sq, axis %d \n", ncycle);
         fprintf(stderr," vin: [%g,%g,%g,%g,%g,%g]\n",
            vin[0],vin[1],vin[2],vin[3],vin[4],vin[5]);
          fprintf(stderr," vi: [%g,%g,%g,%g,%g,%g]\n",
                  vi[0],vi[1],vi[2],vi[3],vi[4],vi[5]);
         return( false );
      }

      if ( debug )
      {
         //printf( "%d    %f %f %f %f %f %f\n", ncycle, vout[0], vout[1], vout[2], vout[3], vout[4], vout[5] );
      }

      ++ncycle;
   }

   return( true );

} // end of Reduce


//-----------------------------------------------------------------------------
// Name: NearRed()
// Description: test whether a vector is nearly Niggli reduced
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool Reducer::NearRed( const arma::vec6& gvec, const double delta )
{
//C
//C     RETURNS .true. IF GVEC IS NEARLY NIGGLI REDUCED
//C     ALLOWING A NON_REDUCTION ERROR OF DELTA.  NO
//C     MATRICES OR VECTORS ARE KEPT.
//C
//C     IF DELTA .EQ. 0.D0, THE TESTS ARE ON REDUCTION
//C     RATHER THAN NEAR REDUCTION
//C
//C     ALL CASES OF BEING ON THE WRONG SIDES OF THE
//C     FOLDING BOUNDARIES ARE ACCEPTED AS NEAR
//C     REDUCED
//C----------------------------------------------------------------------C
//C     TEST FOR G1, G2, G3 OUT OF BOUNDS OR WRONG ORDER
//C
   if ( gvec[0] < -delta        ||
        gvec[1] < -delta        ||
        gvec[2] < -delta        ||
        gvec[0] > gvec[1]+delta ||
        gvec[1] > gvec[2]+delta )
   {
      return( false );
   }
//C
//C     TEST FOR NON-REDUCED SIGN COMBINATIONS IN
//C     G4, G5 AND G6
//C
   if ( (gvec[3] <= delta ||
         gvec[4] <= delta ||
         gvec[5] <= delta ) &&
        (gvec[4] >  delta ||
         gvec[4] >  delta ||
         gvec[5] >  delta ) )
   {
      return( false );
   }
//C
//C     TEST ABS(G{4,5,6}) AGAINST G{1,2,3}
//C
   if ( fabs(gvec[3]) > fabs(gvec[2])+delta ||
        fabs(gvec[4]) > fabs(gvec[0])+delta ||
        fabs(gvec[5]) > fabs(gvec[0])+delta )
   {
      return( false );
   }
//C
//C     TEST THE BODY DIAGONAL
//C
   if ( gvec[3]+gvec[4]+gvec[5] + fabs(gvec[0])+fabs(gvec[1])+delta < 0.0 ) return( false );
   if ( delta > 0.0 ) return( true );
//C
//C     TEST THE 678, 9AB, CDE BOUNDARY FOLDS
//C
   if ( (gvec[3]         == gvec[1] && 2.0*gvec[4] < gvec[5] ) ||
        (gvec[4]         == gvec[0] && 2.0*gvec[3] < gvec[5] ) ||
        (gvec[5]         == gvec[0] && 2.0*gvec[3] < gvec[4] ) ||
        (gvec[3]+gvec[1] == 0.0     && gvec[5]     <= 0.0 )    ||
        (gvec[4]+gvec[0] == 0.0     && gvec[5]     <= 0.0 )    ||
        (gvec[5]+gvec[0] == 0.0     && gvec[4]     <= 0.0 ) ) return( false );
//C
//C     TEST THE F BOUDARY FOLD
//C
   if ( fabs(gvec[3]+gvec[4]+gvec[5]+gvec[0]+gvec[1]) <= delta ) return( false );

   return( true );
} // end NearRed


//-----------------------------------------------------------------------------
// Name: NearRed()
// Description: test whether a vector is nearly Niggli reduced
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
bool Reducer::Near2Red( const arma::vec6& gvec, const double delta, arma::vec6& vout, double& dist )
{
//C
//C     RETURNS .true. IF GVEC IS NEARLY NIGGLI REDUCED
//C     ALLOWING A NON_REDUCTION ERROR OF DELTA.
//C     A VECTOR VOUT IS PRODUCED WHICH IS REDUCED, IF
//C     POSSIBLE.  NO MATRIX IS KEPT.  DIST IS THE DISTANCE
//C     FROM GVEC TO VOUT.
//C
//C     IF DELTA .EQ. 0.D0, THE TESTS ARE ON REDUCTION
//C     RATHER THAN NEAR REDUCTION
//C
//C     ALL CASES OF BEING ON THE WRONG SIDES OF THE
//C     FOLDING BOUNDARIES ARE ACCEPTED AS NEAR
//C     REDUCED
//C
//C----------------------------------------------------------------------C

   const double xdelta = 1.0E-6 * delta;
//C
//C     TEST FOR G1, G2, G3 OUT OF BOUNDS OR WRONG ORDER
//C
   if ( gvec[0] < -xdelta ||
        gvec[1] < -xdelta ||
        gvec[2] < -xdelta ||
        gvec[0] > gvec[1]+xdelta ||
        gvec[1] > gvec[2]+xdelta )
   {

      if ( gvec[0] < 0.0 ) vout[0] = 0.0;
      if ( gvec[1] < 0.0 ) vout[1] = 0.0;
      if ( gvec[2] < 0.0 ) vout[2] = 0.0;
      if ( vout[2]> std::max(vout[0],vout[1] ) )
      {
         if ( vout[1] < vout[0] )
         {
            vout[1] = (vout[0]+vout[1]) / 2.0;
            vout[0] = vout[1];
         }
      }
      else if ( vout[0] < std::min( vout[1],vout[2] ) )
      {
         if ( vout[2] < vout[1] )
         {
            vout[2] = (vout[1]+vout[2]) / 2.0;
            vout[1] = vout[2];
         }
         else
         {
            vout[0] = (vout[0]+vout[1]+vout[2]) / 3.0;
            vout[1] = vout[0];
            vout[2] = vout[0];
         }
      }
   }
//C
   //C     TEST FOR NON-REDUCED SIGN COMBINATIONS IN
   //C     G4, G5 AND G6
   //C
   double s456 = 1.0;
   if ( gvec[3] <= xdelta ||
        gvec[4] <= xdelta ||
        gvec[5] <= xdelta )
   {
      s456 = -1.0;
      if ( gvec[3] > 0.0 ) vout[3] = 0.0;
      if ( gvec[4] > 0.0 ) vout[4] = 0.0;
      if ( gvec[5] > 0.0 ) vout[5] = 0.0;
   }

//C     TEST ABS(G{4,5,6}) AGAINST G{1,2,3}
//C
   if ( fabs(gvec[3]) > fabs(gvec[1])+xdelta ||
        fabs(gvec[4]) > fabs(gvec[0])+xdelta ||
        fabs(gvec[5]) > fabs(gvec[0])+xdelta )
   {
      vout[3] = s456 * std::min( fabs(gvec[1]) , fabs(vout[3]) );
      vout[4] = s456 * std::min( fabs(gvec[0]) , fabs(vout[4]) );
      vout[5] = s456 * std::min( fabs(gvec[0]) , fabs(vout[5]) );
   }
//C
//C     TEST THE BODY DIAGONAL
//C
   if ( gvec[3]+gvec[4]+gvec[5]+fabs(gvec[0])+fabs(gvec[1])+xdelta < 0.0 )
   {
      vout[5] = -vout[0]-vout[1]-vout[3]-vout[4];
   }

//C     IF DELTA IS NON-ZERO, WE STOP HERE

   if ( delta <= 0.0 )
   {
//C
//C     TEST THE 678, 9AB, CDE BOUNDARY FOLDS
//C
   if ( (gvec[3]         ==gvec[1] && 2.0*gvec[4] < gvec[5]) ||
        (gvec[4]         ==gvec[0] && 2.0*gvec[3] < gvec[5]) ||
        (gvec[5]         ==gvec[0] && 2.0*gvec[3] < gvec[4]) ||
        (gvec[5]         == 0.0    &&     gvec[4] <= 0.0 ) ||
        (gvec[3]+gvec[1] == 0.0    &&     gvec[5] <= 0.0 ) ||
        (gvec[4]+gvec[0] == 0.0    &&     gvec[5] <= 0.0) )
        return( false );
//C
//C     TEST THE F BOUNDARY FOLD
//C
   if ( fabs(gvec[3]+gvec[4]+gvec[5]+gvec[0]+gvec[1]) <= xdelta &&
      2.0*(gvec[0]+gvec[4])+gvec[5] > 0.0 )
      return( false );
   }
   dist = BasicDistance::g456dist( gvec, vout );

//C      if (DIST.GT.DELTA) near2red = .false.
//C      write(*,*)"NEAR2RED",delta,near2red,dist
//C      call printg6("GVEC",GVEC)
//C      call printg6("Vout",vout)
   return( dist <= delta );
}
