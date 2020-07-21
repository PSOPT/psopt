/*********************************************************************************************

This file is part of the DMatrix library, a C++ tool for numerical linear algebra

Copyright (C) 2009-2020 Victor M. Becerra

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA,
or visit http://www.gnu.org/licenses/

Author:    Professor Victor M. Becerra
Address:   University of Portsmouth
           School of Energy and Electronic Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/

#ifdef UNIX


extern "C" {


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <cstdio>

}

#else

#ifdef WIN32
#include <stdio.h>
#include <math.h>
#include <string.h>
#else
#include <stdio.h>
#include <math.h>
#include <cstdio>
#endif

#endif


#ifdef SPARSE_MATRIX
extern "C" {

#include "commonlib.h"
#include "myblas.h"
#include "lusol.h"
#include "lusolio.h"
#include "lusolmain.h"

}

#endif


extern "C" {

#include <stdlib.h>

}


#ifdef UNIX

#include "dmatrixv.h"

#ifdef MATLAB_MEX_FILE

#include "mex.h"

#endif  /*MEX*/

#else

#include "dmatrixv.h"

#ifdef MATLAB_MEX_FILE

extern "C" {

#include "mex.h"

}

#endif

#endif



#ifdef MATLAB_MEX_FILE

#define ERROR_MESSAGE error_message

#define PRINTF mexPrintf


#else

#define ERROR_MESSAGE error_message

#define PRINTF printf

#endif


#ifndef MAX
#define MAX(a, b) ( (a)>(b)?  (a):(b) )
#endif
#ifndef MIN
#define MIN(a, b) ( (a)<(b)?  (a):(b) )
#endif


DMatrix*  DMatrix::GetTempPr(int i) {


   auxPr[i].mtype = 0;
   auxPr[i].mt    = NULL;
   auxPr[i].rowIndx = NULL;
   auxPr[i].colIndx = NULL;
   return &auxPr[i];

}

inline long ChkTmpIndx( long taindx )
{


   if ( taindx >= N_TEMP_OBJECTS )

     ERROR_MESSAGE(" Temporary arrays error: Increase N_TEMP_OBJECTS");

#ifdef DEBUG_TEMPS

   printf("\n Temp Created --> indx: %d", taindx );

#endif

   return  taindx;


}



/* Definition of static member variables of class DMatrix */


DEC_THREAD DMatrix* DMatrix::auxPr = NULL;        // Pointer to auxiliary arrays

DEC_THREAD int   DMatrix::initFlag = 0;           // To be allocated in AllocateAuxArr() */

DEC_THREAD int   DMatrix::noAuxArr=N_TEMP_OBJECTS;// Number of auxiliary arrays

DEC_THREAD long  DMatrix::dimAux=D_TEMP_OBJECTS;  // Dimension of auxiliary arrays

DEC_THREAD int   DMatrix::auxIndx=-1;             // Index of used auxiliary arrays

DEC_THREAD int   DMatrix::memberFlag = 0;         // Member function flag

DEC_THREAD int   DMatrix::errorFlag  = 0;	       // Error flag

DEC_THREAD int   DMatrix::print_level = 1;        // Print level

DEC_THREAD clock_t DMatrix::start_clock = 0;

DEC_THREAD double DMatrix::MACH_EPS = 0;

DEC_THREAD time_t DMatrix::start_time = 0;

DEC_THREAD long* DMatrix::seed = NULL;

// long DMatrix::seed[RAND_STREAMS] = {RAND_DEFAULT};

long generate_seed(time_t *timer);


DEC_THREAD int  DMatrix::stream        = 0;                    // stream index


#ifdef OTISS_CASE

DEC_THREAD double DMatrix::axp[N_TEMP_OBJECTS][D_TEMP_OBJECTS];

#endif


char* num2str(double num);

long generate_seed(time_t *timer)
{

   long seed;
   tm* st = localtime(timer);

   seed	 = st->tm_yday*365*24*3600;
   seed += st->tm_hour*3600;
   seed += st->tm_min*60;
   seed += st->tm_sec;

   seed *=37;

   return seed;
}


/* Declaration of Dummy object to allocate and de-allocate   */
/* auxiliary arrays                                           */

#ifndef MATLAB_MEX_FILE
// #ifndef DECLARED_TEMPS

InitializeDMatrixClass Dummy;

// #endif
#endif

#define NR_END 1

#define FREE_ARG char*

void DMatrix::initVars()
{
	a = NULL;
	n = 0;
	m = 0;
	asize = 0;
	atype = 0;
	mtype = 0;
	auxFlag = 0;
	allocated = 0;
	DMatrix* mt = NULL;
	rowIndx = NULL;
	colIndx = NULL;
}


DMatrix::DMatrix(void)
// Default constructor
{
    
    initVars();

    SetReferencedDMatrixPointer( NULL );
    SetRowIndexPointer(NULL);
    SetColIndexPointer(NULL);
    SetMType(0);
    allocated = false;

}



DMatrix::DMatrix(long Initn, long Initm)
//  Matrix constructor using memory allocation
{
    
  initVars();

  n=Initn;

  m=Initm;

  asize = n*m;

  if ( Initm<0 || Initn<0 ) {

      ERROR_MESSAGE("Attempt to create a matrix with \
                     negative dimension in \
                     DMatrix::DMatrix(long,long)");
  }

  if (Initm!=0 && Initn!=0) {

    DMatrix::Allocate(asize);

    memset(a, 0, asize*sizeof(double) );


  }


  SetReferencedDMatrixPointer( NULL );
  SetRowIndexPointer(NULL);
  SetColIndexPointer(NULL);
  SetMType(0);

}

DMatrix::DMatrix( long vDim, double* v, long Initn, long Initm )
// Matrix constructor using pre-allocated array or pointer
{

   initVars();
   
   n = Initn;

   m = Initm;

   a = v;

   if ( Initm<0 || Initn<0 ) {

      ERROR_MESSAGE("Attempt to create a matrix with  \
                     negative dimension in \
                     DMatrix::DMatrix(long, double*, long,long)");
   }
   if ( n*m>vDim ) {

      ERROR_MESSAGE("Error in DMatrix::DMatrix(long,double*,long,long):\
                     increase dimensions of pre-allocated array or pointer");
      fprintf(stderr, "Requested size %ld elements", n*m );
      fprintf(stderr, "Pre-allocated size %ld elements", vDim );
   }

   asize = vDim;

   atype = 1;
   auxFlag=0;

   SetReferencedDMatrixPointer( NULL );
   SetRowIndexPointer(NULL);
   SetColIndexPointer(NULL);
   SetMType(0);
   allocated = false;

}



DMatrix::DMatrix( long Initn)
{
// Matrix constructor with Initn rows and only one column

  if ( Initn <0 ) {

      ERROR_MESSAGE("Attempt to create a matrix with  \
                     negative dimension in \
                     DMatrix::DMatrix(long)");
  }

  initVars();
  
  n=Initn;

  m=1;

  asize = n;

  if ( Initn!=0) {

    DMatrix::Allocate(n);

    memset(a, 0, n*sizeof(double) );

  }

  atype = 0;
  auxFlag = 0;


  SetReferencedDMatrixPointer( NULL );
  SetRowIndexPointer(NULL);
  SetColIndexPointer(NULL);
  SetMType(0);

}


DMatrix::DMatrix(long rows,long columns,double a11,...)
{
// Matrix constructor with assigment

   long i,j;
   
   initVars();

   n = rows;

   m = columns;

   asize = n*m;

   atype = 0;


   Allocate(asize);

   va_list List;

   va_start(List,a11);

   elem(1,1) = a11;

   for (i=1; i<=n; i++) {

     for (j=1; j<=m; j++) {

        if (i*j != 1) {

          elem(i,j) = va_arg(List, double);

        }

     }

   }

   va_end(List);


   SetReferencedDMatrixPointer( NULL );
   SetRowIndexPointer(NULL);
   SetColIndexPointer(NULL);
   SetMType(0);
   auxFlag=0;

}



DMatrix::DMatrix( const DMatrix& A)
// copy constructor
{

  initVars();
  
  n=A.n;

  m=A.m;

  asize = n*m;

  DMatrix::Allocate(asize);

  memcpy( a, A.a , n*m*sizeof(double) );

  atype = 0;


  SetReferencedDMatrixPointer( NULL );
  SetRowIndexPointer(NULL);
  SetColIndexPointer(NULL);
  SetMType(0);
  auxFlag=0;

  if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

}

void DMatrix::PrintInfo(const char *str) const
{
    fprintf(stderr,"\nInformation about DMatrix object: %s", str);
    fprintf(stderr,"\nNumber of Rows    :\t\t%ld",GetNoRows() );
    fprintf(stderr,"\nNumber of Columns :\t\t%ld",GetNoCols() );
    fprintf(stderr,"\nAllocation type   :\t\t%d",atype);
    fprintf(stderr,"\nAllocation size   :\t\t%ld",asize);
}


void DMatrix::AllocateAuxArr( void )
{
// Routine to allocate the temporary objects at program start

   int i;


   // const DEC_THREAD double DMatrix::MACH_EPS = MachEps();

    DMatrix::MACH_EPS = MachEps();

    DMatrix::start_time = time(NULL);

    DMatrix::seed = (long*) my_calloc(RAND_STREAMS, sizeof(long));

    for (i=0;i<RAND_STREAMS;i++ ) {
       DMatrix::seed[i] = generate_seed(&DMatrix::start_time);
    }


   if (DMatrix::GetInitFlag()==0 ) {

    *(DMatrix::GetAuxPr()) = (DMatrix*) my_calloc( N_TEMP_OBJECTS , sizeof(DMatrix) );

    for (i =0; i< GetNoAuxArr(); i++ )
    {

#ifdef DECLARED_TEMPS

      DMatrix::auxPr[i].a = axp[i];
      DMatrix::auxPr[i].atype = 1;



#endif /* DECLARED_TEMPS */

      DMatrix::auxPr[i].Resize(GetDimAux(), 1 );
	  DMatrix::auxPr[i].auxFlag = 1;
	  DMatrix::auxPr[i].asize = GetDimAux();

    }

    DMatrix::SetInitFlag( 1 );

   } // End if


}

void DMatrix::DeAllocateAuxArr( void )
{
// Routine to deallocate temporary objects at program termination

    int i;
    
    if (DMatrix::seed) {
        mxFree(DMatrix::seed);
        DMatrix::seed = NULL;
    }

    for ( i=0; i< GetNoAuxArr(); i++ )
    {

       DMatrix::auxPr[i].~DMatrix();

    }



    mxFree( *(DMatrix::GetAuxPr()) );


}


DMatrix& reshape(DMatrix& A, long nnrow, long nncol) {

    DMatrix* Temp;

    Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

    Temp->Resize( A.GetNoRows(), A.GetNoCols() );

    *Temp = A;

    if (A.GetNoRows()*A.GetNoCols() != nnrow*nncol) {

      error_message("reshape(): The number of elements must be preserved");

    }

    Temp->Resize(nnrow,nncol);

    return (*Temp);

}


void DMatrix::Resize( long nnrow, long nncol )
{
// Resizes the row and column dimensions of a DMatrix object
// additional memory is allocated if necessary
   long i;

   if (nnrow==0 || nncol==0) {
	   n = nnrow;
	   m = nncol;
	   if (atype==0 && auxFlag!=1) {
	     DMatrix::DeAllocate();
	     asize = 0;
	   }
	   return;
   }

   if (n==nnrow && m==nncol ) return;

   if ( n*m == nnrow*nncol)
   {
          n = nnrow;
		  m = nncol;
		  return;
   }

   if (nnrow < 0 || nncol < 0 ) { error_message("Negative arguments in Resize()"); }

   if (asize >= nnrow*nncol && atype==1 )
   {
	   for (i=nnrow*nncol+1;i<=asize;i++) {
		   a[i-1]=0.0;
	   }
	   n=nnrow;
	   m=nncol;

	   return;
   }


   if ( a!=NULL && atype!=1 && auxFlag!=1 )
   {

	 double* atemp = (double *) my_calloc( n*m, sizeof(double) );

         if (!atemp) {
	      ERROR_MESSAGE("\nError allocating double array in DMatrix::Resize()");
         }

         memcpy( atemp, a, n*m*sizeof(double) );

         DMatrix::DeAllocate();

         asize = nnrow*nncol;

	 DMatrix::Allocate(asize);

         memcpy( a, atemp, std::min( n*m, asize )*sizeof(double) );

         n = nnrow;

         m = nncol;

         mxFree( atemp ) ;

	 return;

   }

   if ( a!=NULL && auxFlag==1 ) {

	   if (nnrow*nncol <= dimAux )
	   {
		  n = nnrow;
		  m = nncol;
	   }

	   else {
			ERROR_MESSAGE("\nError resizing Temporary array in DMatrix::Resize()");
	   }

   }


   if ( a==NULL && atype==0 )
   {
	   n = nnrow;
	   m = nncol;
	   asize = nnrow*nncol;
	   DMatrix::Allocate(asize);
	   return;
   }


   if ( atype == 1 && nnrow*nncol > asize  )
   {

       ERROR_MESSAGE("Size error in DMatrix::Resize");

   }


}


DMatrix& DMatrix::Column( long icol ) const
{
// Extracts one column from a DMatrix object
// Eg.  C = A.Column(1);

  DMatrix* Temp;

  int i;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n , 1 );

  for ( i= 1; i<= n ; i++ ) {

     (*Temp)(i) = elem( i, icol );

  }
  return *Temp;

}

DMatrix& DMatrix::Row( long irow ) const
{
// Extracts one row from a DMatrix object
// Eg. R = A.Row( 1 );

  DMatrix* Temp;

  int i;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( 1 , m );

  for ( i= 1; i<= m ; i++ ) {

     (*Temp)(1,i) = elem( irow, i );

  }

  return *Temp;

}

void DMatrix::SwapRows( int i , int j )
{
// Swaps two rows of a DMatrix object
// Eg.  A.SwapRows( 1,2 );


    if ( (i<1) || (i>n) || (j<1) || (j>n) )
    {

        error_message( "Range error in SwapRows() ");

    }

    for (long ii = 1; ii<= m; ii++ ) {

      Swap( &elem(i,ii) , &elem( j, ii ) );

    }

}

void DMatrix::SwapColumns( int i , int j )
{
// Swaps two columns of a DMatrix object
// Eg.  A.SwapColumns( 1,2 );

    if ( (i<1) || (i>m) || (j<1) || (j>m) )
    {

        error_message( "Range error in SwapColumns() ");

    }

    for (long  ii = 1; ii<= n; ii++ ) {

      Swap( &elem(ii, i) , &elem(ii, j) );

    }

}

void DMatrix::diag( const DMatrix& dd )
{
// Assigns values to the diagonal of a DMatrix object
// Eg. A.diag( V );

   int i;

   FillWithZeros();

   if ( n!=m || dd.GetNoRows()!=n || dd.GetNoCols()!=1 ) {

      error_message( "Dimension error in DMatrix::diag() ");

   }

   for ( i=1; i<= n; i++ ) {

      elem( i, i ) = dd(i);

   }

}

void DMatrix::SetColumn( const DMatrix& dd, int icol )
{
// Sets values to a column of a DMatrix object
// Eg. A.SetColumn( 1, C );

   int i;

   if (  dd.GetNoRows()!=n || dd.GetNoCols()!=1 ) {

      error_message( "Dimension error in DMatrix::SetColumn() ");

   }

   for ( i=1; i<= n; i++ ) {

      elem( i, icol ) = dd(i);

   }

}

void DMatrix::SetRow( const DMatrix& dd, int irow )
{
// Sets values to a row of a DMatrix object
// Eg. A.SetRow( 1, R );

   int i;

   if (  dd.GetNoRows()!=m || dd.GetNoCols()!=1 ) {

      error_message( "Dimension error in DMatrix::SetRow() ");

   }

   for ( i=1; i<= m; i++ ) {

      elem( irow, i ) = dd(i);

   }

   if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

}

void DMatrix::rowMult( long r, double x )
{
// Multiplies a row of a DMatrix object by a scalar
// Eg. A.rowMult( 1, 10.0 );

   int i;

   if (  r < 1 || r> n ) {

      error_message( "Range error in DMatrix::rowMult() ");

   }

   for ( i=1; i<= m; i++ ) {

      elem( r, i ) *= x;

   }

}

void DMatrix::colMult( long c, double x )
{
// Multiplies a column of a DMatrix object by a scalar
// Eg. A.colMult( 1, 10.0 );

   int i;

   if (  c < 1 || c> m ) {

      error_message( "Range error in DMatrix::colMult() ");

   }

   for ( i=1; i<= n; i++ ) {

      elem( i ,c ) *= x;

   }

}

void DMatrix::assign(long rows,long columns,double a11,...)
{
// Assign values to a DMatrix object.
// the object is resized if necessary
// Eg. A.assign( 2,2 , 1.0, 0.0,
//		      -1.0, 4.0  );

   long i,j;

   DMatrix::Resize(rows, columns);

   va_list List;

   va_start(List,a11);

   elem(1,1) = a11;

   for (i=1; i<=n; i++) {

     for (j=1; j<=m; j++) {

        if (i*j != 1) {

          elem(i,j) = va_arg(List, double);

        }

     }

   }

   va_end(List);

}




void DMatrix::MemCpyArray(double * aptr )
{
// Copies the elements of an array to the array used for
// matrix storage in a DMatrix object
// Eg.  double v[4] = { 1.0, 2.0, 4.0 , 9.0 }
//      DMatrix A(2,2);
//	...
//      A.MemCpyArray( v );

  memcpy( a, aptr, n*m*sizeof(double) );

}

void DMatrix::input_matrix()
{
// Reads a matrix from the standard input
  long i,j;

  for (i=1;i<=n;i++)

  {

    for (j=1;j<=m;j++)

    {

      PRINTF("Input DMatrix(%ld,%ld):",i,j);

      scanf("%lf", &elem( i, j) );

    }

  }

}


void DMatrix::SetPrintLevel( int plevel )
{
	 print_level=plevel;
	 return;
}


int DMatrix::PrintLevel() {
	return print_level;
}



void DMatrix::Print(const char* text) const
{
// Prints a matrix
// Eg. A.Print(" A = ");


  long i,j;
  if (print_level) {
#ifndef MATLAB_MEX_FILE
  fprintf(OUTPUT_STREAM,"\n%s\n",text);
  for (i=1;i<=n;i++)
  {
    for (j=1;j<=m;j++)
    {
      fprintf(OUTPUT_STREAM," %s ", num2str(elem(i,j)) );
      fprintf(OUTPUT_STREAM,"\t");
    }
    fprintf(OUTPUT_STREAM,"\n");
  }
#else
  mexPrintf("\n%s\n",text);
  for (i=1;i<=n;i++)
  {
    for (j=1;j<=m;j++)
    {
      mexPrintf(" %s ", num2str(elem(i,j)) );
      mexPrintf("\t");
    }
    mexPrintf("\n");
  }
#endif
  }
  if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

}


void error_message(const char *error_text)
{
// Error handling routine

  string m1;
  string m2;

  m1 =  "\n**** Run time error";
  m1 += "\n**** To trace this error, set up your debugger to break at";
  m1 += "\n**** function 'error_message', then do a backtrace.";
  m1 += "\n**** A diagnostic message is given below:";
  m2 = error_text;
  m2 =  "\n**** ====> " + m2 + " <====\n\n";

  if ( DMatrix::PrintLevel() ) {
#ifndef MATLAB_MEX_FILE
  fprintf(OUTPUT_STREAM,"%s", m1.c_str());
  fprintf(OUTPUT_STREAM,"%s", m2.c_str());
  FILE* err_file = fopen("error_message.txt","w");
  fprintf(err_file,"%s", m1.c_str() );
  fprintf(err_file,"%s", m2.c_str() );
  fclose(err_file);
#else
  mexPrintf("\n%s\n",m1.c_str());
  mexPrintf("%s\n",m2.c_str());
#endif
 }

  DMatrix::RiseErrorFlag();

  throw ErrorHandler(m1+m2);


}

ErrorHandler::ErrorHandler(const string m)
{
  error_message = m;
}

DMatrix& DMatrix::operator+ (const DMatrix& Other_matrix) const
{
// Matrix addition operator: Adds two matrices
// Example C = D + E; where C,D,E are DMatrix objects

  long i;

  const double* OtherPr = Other_matrix.GetConstPr();

  DMatrix* Temp;

  long nelem = n*m;

  if (m != Other_matrix.m || n != Other_matrix.n)

  { ERROR_MESSAGE("Incoherent matrix dimensions in addition"); }


  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );


  for (i=0; i< nelem; i++)
  {

     Temp->GetPr()[i] = a[i] + OtherPr[i];

  }

  return *Temp;

}


DMatrix& DMatrix::operator+ (double arg) const
{
// Matrix  + scalar operator: Adds scalar to each element of the matrix.
// Example C = D + f; where C,D are DMatrix objects and f is a double

  long i;


  DMatrix* Temp;

  long nelem = n*m;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );


  for (i=0; i< nelem; i++)
  {

     Temp->GetPr()[i] = a[i] + arg;

  }

  return *Temp;

}

DMatrix& DMatrix::operator- (double arg) const
{
// Matrix  - scalar operator: Subtracts scalar to each element of the matrix.
// Example C = D - f; where C,D are DMatrix objects and f is a double

  long i;



  DMatrix* Temp;

  long nelem = n*m;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );


  for (i=0; i< nelem; i++)
  {

     Temp->GetPr()[i] = a[i] - arg;

  }

  return *Temp;

}




DMatrix& DMatrix::operator += (const DMatrix &rval)
{
// Matrix addition with substitution
// Eg. A += B;  ***** Equivalent to A = A+B; ******

   long nelem = n*m;

   if(m != rval.m || n != rval.n )
   {

     ERROR_MESSAGE("Bad dimensions in DMatrix op. +=");

   }

   for (long i=0; i< nelem; i++ )
   {

      a[i] += rval.a[i];

   }

   if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

   return *this;

}




DMatrix& DMatrix::operator- (const DMatrix& Other_matrix) const
{
// Matrix substraction operator: substracts two matrices
// Example C = D - E; where C,D,E are DMatrix objects

  long i;

  long nelem = n*m;

  const double *OtherPr = Other_matrix.GetConstPr();

  DMatrix* Temp;

  if (m != Other_matrix.m || n != Other_matrix.n)

  { ERROR_MESSAGE("Incoherent matrix dimensions in substraction"); }

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  for (i=0; i< nelem; i++)
  {

     Temp->GetPr()[i] = a[i] - OtherPr[i];

  }

  return *Temp;
}

DMatrix& DMatrix::operator -= (const DMatrix &rval)
{
// Matrix substraction with substitution
// Eg. A -= B;  ***** Equivalent to A = A-B; ******

   long nelem = n*m;

   if(n != rval.GetNoRows() || m != rval.GetNoCols() )
   {

     ERROR_MESSAGE("Bad dimensions in DMatrix op. +=");

   }

   for (long i=0; i< nelem; i++ )
   {

      a[i] -= rval.GetConstPr()[i];

   }

   if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

  return *this;

}


DMatrix& operator- (const DMatrix& A)
{
// Matrix unary minus operator
// Eg. C =  -A;


  long nelem = A.GetNoRows()*A.GetNoCols();
  const double *Apr = A.GetConstPr();
  DMatrix* Temp;




  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( A.GetNoRows(), A.GetNoCols() );


  for(long i=0; i< nelem; i++)
  {

     Temp->GetPr()[i] = -Apr[i];

  }

  return *Temp;

}



DMatrix& DMatrix::operator* (const DMatrix& B) const
{
// Matrix product operator, multiplies two matrices
// Eg. C = A*B;

   return Product( *this, B );

}

DMatrix& Product(const DMatrix& A, const DMatrix& B)
{
// Matrix product optimised for vector columnwise storage

  long  i,j,k,na,nb,mb;

  long  rowA, colB;

  double sum;

  DMatrix *Temp;

  na=A.n;

  mb=B.m;

  nb=B.n;

  const double* apr = A.GetConstPr();

  const double* bpr = B.GetConstPr();

  if (A.m != nb )

  { ERROR_MESSAGE("Incoherent matrix dimensions in product"); }


  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( A.n, mb );

  for (i=0; i< na; i++)
  {


    for (j=0; j<mb; j++)
    {

      colB = j*B.n;

      sum=0.0;

      rowA = i;

      for (k=0; k<nb; k++)
      {

          sum += apr[ rowA ]*bpr[ colB + k ];
          rowA += A.n;

      }

      Temp->elem(i+1,j+1) = sum;

    }

  }

  return *Temp;

}

DMatrix& crossProduct(const DMatrix& x, const DMatrix& y)
{

  if ( !x.isVector() || !y.isVector() )
      error_message("Both arguments must be vectors in crossProduct()");

  if (length(x) != 3 ) {
     error_message("length of first argument is not 3 in crossProduct()");
  }

  if (length(y) != 3 ) {
     error_message("length of second argument is not 3 in crossProduct()");
  }

  DMatrix *Temp;

  DMatrix A(3,3);

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));


  Temp->Resize( x.n, x.m );

   A.FillWithZeros();

  A(1,2) = -x(3);  A(1,3)= x(2);
  A(2,1) =  x(3);  A(2,3)=-x(1);
  A(3,1) = -x(2);  A(3,2)= x(1);

  (*Temp) = A*y;

  return *Temp;

}


DMatrix& DMatrix::operator*= (const DMatrix& Other_matrix)
{
// Matrix product operator with substitution
// Eg. C *= A;     ***** Equivalent to C = C*A *******

  long  i,j,k,na,nb,mb;

  double sum;

  DMatrix* Temp;

  na=n;

  mb=Other_matrix.m;

  nb=Other_matrix.n;

  if (m != nb )

  { ERROR_MESSAGE("wrong matrix dimensions in op. *="); }


  double *pr;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( na, mb );

  pr = Temp->GetPr();

  for (i=1; i<=na; i++)
  {

    for (j=1; j<=mb; j++)
    {

      sum=0.0;

      for (k=1; k<=nb; k++)
      {

           sum += elem(i,k)*Other_matrix.elem(k,j);

      }


      Temp->elem(i,j) = sum;

    }

  }

  memcpy(a,pr, n*m*sizeof(double));

   if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

  return *this;

}





DMatrix& DMatrix::operator* (double Arg) const
{
// Matrix by scalar product.
// Eg. C = A*x, where C & A are DMatrix objects and x is a double
  long i;

  long nelem = n*m;

  DMatrix* Temp;

  double *pr;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  pr = Temp->GetPr();

  for( i=0; i< nelem; i++ )
  {

     pr[i] = a[i]*Arg;

  }

  return *Temp;

}


DMatrix& DMatrix::operator*= (double Arg)
{
// Matrix by scalar product with substitution
// Eg. C *= x, where C is a DMatrix object and x is a double
// Equivalent to: C = C*x;

  long i;

  long nelem = n*m;

  for( i=0; i< nelem; i++ )
  {

     a[i] *= Arg;

  }

  return *this;

}

DMatrix& DMatrix::operator/ (const DMatrix& rval) const
{
// Right division implementation, B/A == B*inv(A)

  int localFlag = 0;

  if (DMatrix::GetMemberFlag()) localFlag = 1;

  DMatrix::SetMemberFlag( 1 );

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, rval.m );

  *Temp = tra( tra(rval)%tra(*this) );

  if(!localFlag) DMatrix::SetMemberFlag( 0 );

  return *Temp;

}


DMatrix& DMatrix::operator/ (double Arg) const
{
// Matrix scalar division.
// Eg. C = A/x, where C & A are DMatrix objects and x is a double

  long i;

  long nelem = n*m;

  DMatrix* Temp;

  double *pr;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  pr = Temp->GetPr();

  for( i=0; i< nelem; i++ )
  {

     pr[i] = a[i]/Arg;

  }

  return *Temp;

}


DMatrix& DMatrix::operator/= (double Arg)
{
// Matrix by scalar division with substitution
// Eg. C /= x, where C is a DMatrix object and x is a double
// Equivalent to: C = C/x;

  long i;

  long nelem = n*m;

  for( i=0; i< nelem; i++ )
  {

     a[i] /= Arg;

  }

  return *this;

}


DMatrix& operator *(double Arg, const DMatrix& A)
{
// Scalar by Matrix product.
// Eg. C = x*A, where C & A are DMatrix objects and x is a double


  long i,n,m,nelem;

  const double * Apr = A.GetConstPr();

  n=A.GetNoRows();

  m=A.GetNoCols();

  nelem = n*m;

  DMatrix* Temp;

  double *pr;


  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  pr   = Temp->GetPr();


  for (i=0; i< nelem; i++ )
  {

     pr[i] = Apr[i]*Arg;

  }

  return *Temp;

}


DMatrix& DMatrix::operator <  ( double val ) const
{
// Relational < operator
// Elementwise comparision
// A DMatrix object is returned contaning zeros and ones
// Eg D = C < x ; where D and C are DMatrix objects and x is a double

  long nelem = n*m;

  long i;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  for ( i=0; i< nelem; i++ )
  {

     if ( a[i] < val ) Temp->a[i] = 1.0;

     else Temp->a[i] = 0.0;

  }

  return *Temp;

}

DMatrix& DMatrix::operator <=  ( double val ) const
{
// Relational <= operator
// Elementwise comparision
// A DMatrix object is returned contaning zeros and ones
// Eg D = C <= x ; where D and C are DMatrix objects and x is a double

  long nelem = n*m;

  long i;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  for ( i=0; i< nelem; i++ )
  {

     if ( a[i] <= val ) Temp->a[i] = 1.0;

     else Temp->a[i] = 0.0;

  }

  return *Temp;

}

DMatrix& DMatrix::operator >  ( double val ) const
{
// Relational > operator
// Elementwise comparision
// A DMatrix object is returned contaning zeros and ones
// Eg D = C > x ; where D and C are DMatrix objects and x is a double

  long nelem = n*m;

  long i;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  for ( i=0; i< nelem; i++ )
  {

     if ( a[i] > val ) Temp->a[i] = 1.0;

     else Temp->a[i] = 0.0;

  }

  return *Temp;

}

DMatrix& DMatrix::operator >=  ( double val ) const
{
// Relational >= operator
// Elementwise comparision
// A DMatrix object is returned contaning zeros and ones
// Eg D = C >= x ; where D and C are DMatrix objects and x is a double

  long nelem = n*m;

  long i;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  for ( i=0; i< nelem; i++ )
  {

     if ( a[i] >= val ) Temp->a[i] = 1.0;

     else Temp->a[i] = 0.0;

  }

  return *Temp;

}

DMatrix& DMatrix::operator ==  ( double val ) const
{
// Relational == operator
// Elementwise comparision
// A DMatrix object is returned contaning zeros and ones
// Eg D = ( C == x ); where D and C are DMatrix objects and x is a double

  long nelem = n*m;

  long i;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  for ( i=0; i< nelem; i++ )
  {

     if ( a[i] == val ) Temp->a[i] = 1.0;

     else Temp->a[i] = 0.0;

  }

  return *Temp;

}



DMatrix& DMatrix::operator = (const char* str)
{
    char cnumber[20];
    long rowcount = 1;
    long colcount = 0;
    long i = 0;
    double value;
    if ( str[0]!='[') {

//        fprintf(stderr,"\n1: str[0] = %c", str[0]);
        error_message("Input error in function DMatrix::operator = (const char* str)");
    }

    int pos = 1;
    int j = 0;
    int dotcount = 0;


    bool nflag = false;

    // First pass to find out the number of rows and columns.

    while(1) {
        if (str[pos]==' ') {
           pos++;
        }
    	while ( str[pos]=='1' || str[pos]=='2' ||str[pos]=='3' || str[pos]=='4' ||str[pos]=='5' || str[pos]=='6' ||str[pos]=='7' || str[pos]=='8' || str[pos]=='9' || str[pos] == '0' || str[pos] == '.' || str[pos]=='-' || str[pos]=='+' || str[pos]=='e' || str[pos]=='E' )
        {
            cnumber[j] = str[pos];
            j++;
            pos++;
            if (str[pos]=='.') {
                 dotcount++;
            }
            if (dotcount > 1) {
               error_message("Input error in function DMatrix::operator = (const char* str)");
            }
            nflag = true;
        }
        if( nflag && dotcount ) {
        	cnumber[j]='\0';
        }
        if( nflag && dotcount==0 ) {
                cnumber[j] = '.';
        	cnumber[j+1]='\0';
        }
        if ( nflag && rowcount==1  ) {
               colcount++;
        }
        j = 0;
        dotcount = 0;
        if (str[pos]==';') {
           rowcount++;
        }
        else if (str[pos]==']') {
             break;
        }
        else if (str[pos] != ' ' && str[pos]!=',') {
            error_message("Input error in function DMatrix::operator = (const char* str)");
        }
        if (nflag) {
            nflag = false;
        }
        pos++;
    }

    this->Resize(rowcount,colcount);

    pos = 1;
    j = 0;
    dotcount = 0;

    nflag = false;

    // Second pass to actually fill the values

    rowcount= 1;
    colcount= 0;


    while(1) {
        if (str[pos]==' ' || str[pos]==',') {
           pos++;
        }
    	while ( str[pos]=='1' || str[pos]=='2' ||str[pos]=='3' || str[pos]=='4' ||str[pos]=='5' || str[pos]=='6' ||str[pos]=='7' || str[pos]=='8' || str[pos]=='9' || str[pos] == '0' || str[pos] == '.' || str[pos]=='-' || str[pos]=='+' || str[pos]=='e' || str[pos]=='E')
        {
            cnumber[j] = str[pos];
            j++;
            pos++;
            if (str[pos]=='.') {
                 dotcount++;
            }
            if (dotcount > 1) {
               error_message("Input error in function DMatrix::operator = (const char* str)");
            }
            nflag = true;
        }
        if( nflag && dotcount ) {
        	cnumber[j]='\0';
        }
        if( nflag && dotcount==0 ) {
                cnumber[j] = '.';
        	cnumber[j+1]='\0';
        }
        if ( nflag  ) {
               colcount++;
        }
        j = 0;
        dotcount = 0;
        if (nflag) {
        	sscanf(cnumber, "%lf", &value);
		this->elem(rowcount,colcount) = value;
        	i++;
                nflag = false;
        }

        if (str[pos]==';') {
           rowcount++;
           colcount = 0;
        }
        else if (str[pos]==']') {
             break;
        }
        else if (str[pos] != ' ' && str[pos] !=',' ) {
            error_message("Input error in function DMatrix::operator = (const char* str)");
        }
        if (nflag) {
            nflag = false;
        }
        pos++;
    }

    return (*this);
}

DMatrix& DMatrix::operator !=  ( double val ) const
{
// Relational != operator
// Elementwise comparision
// A DMatrix object is returned contaning zeros and ones
// Eg D = ( C != x ); where D and C are DMatrix objects and x is a double

  long nelem = n*m;

  long i;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  for ( i=0; i< nelem; i++ )
  {

     if ( a[i] != val ) Temp->a[i] = 1.0;

     else Temp->a[i] = 0.0;

  }

  return *Temp;

}


DMatrix& DMatrix::compMat( const DMatrix& mtx, char op ) const
{
// Elementwise matrix comparision function

  long nelem = n*m;

  long i;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  if ( n!=mtx.n || m!=mtx.m ) error_message("\n Dimension error in operator function DMatrix::compMat()");

  for ( i=0; i< nelem; i++ )
  {
     switch( op ) {
     case 1:
       if ( a[i] > mtx.a[i] ) Temp->a[i] = 1.0;
       else Temp->a[i] = 0.0; break;
     case 2:
       if ( a[i] >= mtx.a[i] ) Temp->a[i] = 1.0;
       else Temp->a[i] = 0.0; break;
     case 3:
       if ( a[i] < mtx.a[i] ) Temp->a[i] = 1.0;
       else Temp->a[i] = 0.0; break;
     case 4:
       if ( a[i] <= mtx.a[i] ) Temp->a[i] = 1.0;
       else Temp->a[i] = 0.0; break;
     case 5:
       if ( a[i] == mtx.a[i] ) Temp->a[i] = 1.0;
       else Temp->a[i] = 0.0; break;
     case 6:
       if ( a[i] != mtx.a[i] ) Temp->a[i] = 1.0;
       else Temp->a[i] = 0.0; break;

     } /* end switch */

  } /* end for */

  return *Temp;

}

DMatrix& DMatrix::operator >  ( const DMatrix& mtx ) const
{
   return this->compMat( mtx, 1 );
}

DMatrix& DMatrix::operator >=  ( const DMatrix& mtx ) const
{
   return this->compMat( mtx, 2 );
}

DMatrix& DMatrix::operator <  ( const DMatrix& mtx ) const
{
   return this->compMat( mtx, 3 );
}

DMatrix& DMatrix::operator <=  ( const DMatrix& mtx ) const
{
   return this->compMat( mtx, 4 );
}

DMatrix& DMatrix::operator ==  ( const DMatrix& mtx ) const
{
   return this->compMat( mtx, 5 );
}

DMatrix& DMatrix::operator !=  ( const DMatrix& mtx ) const
{
   return this->compMat( mtx, 6 );
}

void DMatrix::SetMType( int arg )
{
     if (arg !=0 && arg !=1)  {
            ERROR_MESSAGE("Incorrect argument value in DMatrix::SetMType");
     }
     mtype = arg;
}


int any( const DMatrix& mtx )
{
   int retval;
   if (MaxAbs( mtx ) > 0.0) retval = 1;
   else retval = 0;
   if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );
   return retval;
}

DMatrix& tra(const DMatrix& A)
{
// Matrix transposition, friend implementation
// Eg. C = tra(A);
  long i,j,n,m;

  n=A.n;

  m=A.m;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( m, n );


  for (i=1;i<=n;i++)
  {

    for( j=1;j<=m;j++)
    {

        Temp->elem(j,i) = A.elem(i,j);

    }

  }

  return *Temp;

}

void DMatrix::Transpose(void)
{
// Transposes and modifies a matrix
// Eg. A.Transpose();

  int i,j,tmp;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( m, n );

  for (i=1;i<=n;i++)
  {

    for( j=1;j<=m;j++)
    {

        Temp->elem(j,i) = elem(i,j);

    }

  }

  tmp = n;
  n = m;
  m = tmp;

  memcpy( a, Temp->a, n*m*sizeof(double) );


}



DMatrix& TProduct(const DMatrix& A, const DMatrix& B)
{
// Returns the product tra(A)*B
// Eg. C = Tproduct(A,B);

  long  i,j,k,na,nb,mb;

  double sum;

  DMatrix *Temp;

  na=A.m;

  mb=B.m;

  nb=B.n;

  if (A.n != nb )

  { ERROR_MESSAGE("Incoherent matrix dimensions in TProduct()"); }


  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( A.m, mb );

  for (i=1; i<=na; i++)
  {

    for (j=1; j<=mb; j++)
    {

      sum=0.0;

      for (k=1; k<=nb; k++)
      {

          sum += A.elem(k,i)*B.elem(k,j);

      }

      Temp->elem(i,j) = sum;

    }

  }

  return *Temp;

}

DMatrix& ProductT(const DMatrix& A, const DMatrix& B)
{
// Returns the product A*tra(B)
// Eg. C = productT(A,B);


  long  i,j,k,na,nb,mb;

  double sum;

  DMatrix *Temp;

  na=A.n;

  mb = B.n;

  nb=  B.m;

  if (A.m != nb )

  { ERROR_MESSAGE("Incoherent matrix dimensions in ProductT"); }


  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( A.n, mb );

  for (i=1; i<=na; i++)
  {

    for (j=1; j<=mb; j++)
    {

      sum=0.0;

      for (k=1; k<=nb; k++)
      {

          sum += A.elem(i,k)*B.elem(j,k);

      }

      Temp->elem(i,j) = sum;

    }

  }

  return *Temp;

}

DMatrix& TProductT(const DMatrix& A, const DMatrix& B)
{
// Matrix TProductT optimised for vector storage
// Returns the product tra(A)*tra(B)
// Eg. C = TProductT(A,B);

  long  i,j,k;

  long  colA, rowB;

  double sum;

  DMatrix *Temp;

  const double* apr = A.GetConstPr();

  const double* bpr = B.GetConstPr();

  if (A.n != B.m )

  { ERROR_MESSAGE("Incoherent matrix dimensions in TProductT()"); }


  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( A.m, B.n );

  for (i=0; i< A.m; i++)
  {


    for (j=0; j<B.n; j++)

    {

      colA = i*A.n;

      sum=0.0;

      rowB = j;

      for (k=0; k< A.n; k++)
      {

          sum += apr[ colA+k ]*bpr[ rowB ];
          rowB += B.n;

      }

      Temp->elem(i+1,j+1) = sum;

    }

  }

  return *Temp;

}

DMatrix& inv(const DMatrix& A)
{
// Matrix inversion
// Uses Gauss-Jordan elimination
// Eg. C = inv(A);

  long i,p,k,n,m;

  double* a;

  double pivot;

  DMatrix* Temp;

  n=A.n;

  m=A.m;

  if (n!=m) ERROR_MESSAGE("matrix must be square in inv()");

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  memcpy( Temp->a, A.a, A.n*A.m*sizeof(double) );

  a = Temp->a;

  for (p=0;p<n;p++)
  {

    pivot = a[p*n+p];

    for ( i=0;i<n;i++)
    {

      if ( i!=p )
      {

        for ( k=0; k<n; k++)
        {

          if ( k!=p )
          {

            if (  pivot==0.0) ERROR_MESSAGE("singular matrix");
//            Temp->elem(i,k)=Temp->elem(i,k)-Temp->elem(i,p)*Temp->elem(p,k)/Temp->elem(p,p);
              a[k*n+i] +=  -a[p*n+i]*a[k*n+p]/pivot;

          }

        }

      }

    }

    for (i=0;i<n;i++)
    {

      if (i!=p)
      {

//        Temp->elem(i,p)=-Temp->elem(i,p)/Temp->elem(p,p);
          a[p*n+i] = -a[p*n+i]/pivot;

//        Temp->elem(p,i)=Temp->elem(p,i)/Temp->elem(p,p);
          a[i*n+p]=a[i*n+p]/pivot;


      }

    }


    if (a[p*n+p]==0.0) ERROR_MESSAGE("singular matrix");

//    Temp->elem(p,p)=1/Temp->elem(p,p);
    a[p*n+p]=1/pivot;


  }

  return *Temp;

}



DMatrix& DMatrix::operator= (const DMatrix& Other_matrix)
{
// Matrix Assignment operator
// Eg. A = B;


 if ( GetMType() == 1 ) {

      return (AssignmentToColonReference( Other_matrix ));

 }

 else {

  if (n > 0 && m>0) {

    if (n != Other_matrix.n || m != Other_matrix.m)

       this->Resize(Other_matrix.n,Other_matrix.m);

  }

  else {

         n = Other_matrix.n;

         m = Other_matrix.m;
         if (atype==0 && auxFlag!=1)
           DMatrix::Allocate(n*m);

  }

  memcpy( a, Other_matrix.a, n*m*sizeof(double) );

   if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

  return *this;

 } /* End if */

}



void DMatrix::DeAllocate()
/* free the arrays  allocated with Allocate() */
{

   if (n!=0 && m!=0 ) {

          if (allocated==true) {
		  mxFree((FREE_ARG) (a) );
                  allocated=false;
          }

   }

}


void DMatrix::Read(FILE *filex)
{

  long i,j;

  for (i=1;i<=n;i++) {

     for (j=1;j<=m;j++) {

       fscanf(filex,"%lf",&elem(i,j));

     }
  }

}



DMatrix& identity(long ndim)
{

  long i;

  double *pr;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( ndim, ndim );

  pr   = Temp->GetPr();

  memset(pr, 0, ndim*ndim*sizeof(double) );

  for( i=0; i<ndim; i++ ) {

    pr[ i * ndim + i ] = 1.0;

  }

  return *Temp;

}

DMatrix& identity(long ndim,long mdim)
{

  long i;

  double *pr;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( ndim, mdim );

  pr   = Temp->GetPr();

  memset(pr, 0, ndim*mdim*sizeof(double) );

  for( i=0; i<MIN(ndim,mdim); i++ ) {

    pr[ i * ndim + i ] = 1.0;

  }

  return *Temp;
}

DMatrix& diag(const DMatrix& A )
{
// if A is a matrix this function extracts a column vector
//   with the diagonal values of A
// if A is a vector this function returns a matrix having the
//   elements of A in the diagonal
  long i;

  long ndim;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  if ( A.GetNoRows() == 1 || A.GetNoCols() == 1 ) {

       DMatrix* v;

       v = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

       *v = vec(A);

       ndim = A.GetNoRows()*A.GetNoCols();

       Temp->Resize(ndim,ndim);

       Temp->FillWithZeros();

       for(i=1;i<=ndim;i++) {

          (*Temp)(i,i) = (*v)(i);

       }

  }

  else {

    ndim = MIN( A.GetNoRows(), A.GetNoCols() );

    Temp->Resize( ndim , 1 );

    for( i=0; i<ndim; i++ ) {

      Temp->a[ i ] = A( i+1, i+1 );

    }

  } /* End if */

  return *Temp;


}

DMatrix& vec(const DMatrix& A )
{

  long i;

  int  ndim = A.n*A.m;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( ndim , 1 );

  for( i=0; i<ndim; i++ ) {

    Temp->a[ i ] = A.a[i];

  }

  return *Temp;

}

void DMatrix::FillWithZeros( void )
{

  memset( a , 0, n*m*sizeof(double) );
  if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

}


DMatrix& zeros(long nrows, long ncols)
{


  long nelem = nrows*ncols;

  double *pr;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( nrows, ncols );

  pr = Temp->GetPr();

  memset(pr, 0, nelem*sizeof(double) );

  return *Temp;

}

DMatrix& ones(long nrows, long ncols)
{

  long i;

  long nelem = nrows*ncols;

  double *pr;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( nrows, ncols );

  pr   = Temp->GetPr();

  for ( i=0; i< nelem; i++ ) {

     pr[i] = 1.0;

  }

  return *Temp;

}



DMatrix& expm(const DMatrix& B)
#ifdef EXPM_RECURSIVE
{
   long  k;

   int e;

   long nrow = B.getn();

   long ncol = B.getm();

   int localFlag = 0;

      if (DMatrix::GetMemberFlag()) localFlag = 1;

      DMatrix::SetMemberFlag( 1 );


  DMatrix  *E, *A;

  E    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
  A    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  E->Resize(nrow, ncol);
  A-> Resize(nrow,ncol);

  *A = eye(nrow);

  *E = eye(nrow);


  for(k=1; k<=EXPM_RECURSIVE; k++) {

       *E *= B;

       *A += (*E)/factorial(k);

  }

  if( !localFlag ) DMatrix::SetMemberFlag( 0 );

  return (*A);

}

#else
{
   // Uses pade approximation method to find
   // the exponential matrix.


   double c;

   long q, p, s, k;

   int e;

   long nrow = B.getn();

   long ncol = B.getm();

   int localFlag = 0;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );


  DMatrix  *E, *X, *cX, *D, *A;

  E    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
  X    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
  cX   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
  D    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
  A    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  E->Resize(nrow, ncol);
  X-> Resize(nrow,ncol);
  cX->Resize(nrow,ncol);
  D-> Resize(nrow,ncol);
  A-> Resize(nrow,ncol);

     *A = B;



     if ( nrow!= ncol ) {

       ERROR_MESSAGE("expm: Matrix must be square");

     }


     frexp( InfNorm(*A), &e );

     s = MAX( 0 , e+1 );

     (*A) /= pow(2.0, (double) s);
//     DMatrix::SetAuxIndx( -1 );

     *X = *A;
//     DMatrix::SetAuxIndx( -1 );

     c = 0.5;

     *E = *A;
//     DMatrix::SetAuxIndx( -1 );
     *E *= c;
//     DMatrix::SetAuxIndx( -1 );
     *E += identity( nrow );

//     DMatrix::SetAuxIndx( -1 );

     (*D) = *A;
//     DMatrix::SetAuxIndx( -1 );
     (*D) *= -c;
//     DMatrix::SetAuxIndx( -1 );
     (*D) += identity( nrow );
//     DMatrix::SetAuxIndx( -1 );

     q = 6;

     p = 1;

     for (k =2; k<=q; k++ )
     {

       c = c*(q-k+1) / (k*(2*q-k+1));

       *X = (*A)*(*X);

//       DMatrix::SetAuxIndx( -1 );

       *cX = c*(*X);

//       DMatrix::SetAuxIndx( -1 );

       (*E) += (*cX);

//       DMatrix::SetAuxIndx( -1 );

       if (p)
       {

                (*D) += (*cX);
//                DMatrix::SetAuxIndx( -1 );

       }

       else
       {

                (*D) -= (*cX);
//                DMatrix::SetAuxIndx( -1 );

       }

       if ( p == 0 ) p =1;

       else { p = 0; }

     }

     (*E) = (*D)%(*E);

//     DMatrix::SetAuxIndx( -1 );

     for (k =1; k<= s; k++ )
     {

       (*E) *= (*E);
//       DMatrix::SetAuxIndx( -1 );

     }



   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *E;

}

#endif


DMatrix& sin(const DMatrix& A)
{
   // Element-wise sin of a matrix

   long i;

   long nrow = A.getn();

   long ncol = A.getm();

   int localFlag = 0;

   double* temp_pr;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix  *Temp;

   Temp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize(nrow, ncol);

   temp_pr = Temp->GetPr();

   for( i =0; i< nrow*ncol; i++)
   {

      temp_pr[i] = sin( A.a[i] );

   }


   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *Temp;

}

DMatrix& cos(const DMatrix& A)
{
   // Element-wise cosine of a matrix

   long i;

   long nrow = A.getn();

   long ncol = A.getm();

   int localFlag = 0;

   double* temp_pr;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix  *Temp;

   Temp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize(nrow, ncol);

   temp_pr = Temp->GetPr();

   for( i =0; i< nrow*ncol; i++)
   {

      temp_pr[i] = cos( A.a[i] );

   }


   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *Temp;

}

DMatrix& tan(const DMatrix& A)
{
   // Element-wise tangent of a matrix

   long i;

   long nrow = A.getn();

   long ncol = A.getm();

   int localFlag = 0;

   double* temp_pr;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix  *Temp;

   Temp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize(nrow, ncol);

   temp_pr = Temp->GetPr();

   for( i =0; i< nrow*ncol; i++)
   {

      temp_pr[i] = tan( A.a[i] );

   }


   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *Temp;

}

DMatrix& exp(const DMatrix& A)
{
   // Element-wise exponential of a matrix

   long i;

   long nrow = A.getn();

   long ncol = A.getm();

   int localFlag = 0;

   double* temp_pr;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix  *Temp;

   Temp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize(nrow, ncol);

   temp_pr = Temp->GetPr();

   for( i =0; i< nrow*ncol; i++)
   {

      temp_pr[i] = exp( A.a[i] );

   }


   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *Temp;

}


DMatrix& tanh(const DMatrix& A)
{
   // Element-wise hyperbolic tangent of a real matrix

   long i;

   long nrow = A.getn();

   long ncol = A.getm();

   int localFlag = 0;

   double* temp_pr;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix  *Temp;

   Temp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize(nrow, ncol);

   temp_pr = Temp->GetPr();

   for( i =0; i< nrow*ncol; i++)
   {

      temp_pr[i] = tanh( A.a[i] );

   }


   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *Temp;

}

DMatrix& sinh(const DMatrix& A)
{
   // Element-wise hyperbolic sine of a real matrix

   long i;

   long nrow = A.getn();

   long ncol = A.getm();

   int localFlag = 0;

   double* temp_pr;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix  *Temp;

   Temp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize(nrow, ncol);

   temp_pr = Temp->GetPr();

   for( i =0; i< nrow*ncol; i++)
   {

      temp_pr[i] = sinh( A.a[i] );

   }


   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *Temp;

}

DMatrix& cosh(const DMatrix& A)
{
   // Element-wise hyperbolic cosine of a real matrix

   long i;

   long nrow = A.getn();

   long ncol = A.getm();

   int localFlag = 0;

   double* temp_pr;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix  *Temp;

   Temp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize(nrow, ncol);

   temp_pr = Temp->GetPr();

   for( i =0; i< nrow*ncol; i++)
   {

      temp_pr[i] = cosh( A.a[i] );

   }


   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *Temp;

}


DMatrix& log(const DMatrix& A)
{
   // Element-wise natural logarithm of a matrix
   // Elements of A must be positive

   long i;

   long nrow = A.getn();

   long ncol = A.getm();

   int localFlag = 0;

   double* temp_pr;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix  *Temp;

   Temp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize(nrow, ncol);

   temp_pr = Temp->GetPr();

   for( i =0; i< nrow*ncol; i++)
   {

      if (A.a[i]<0.0) error_message("Attempt to compute the natural logarithm of a negative number in log(const DMatrix&)");

      temp_pr[i] = log( A.a[i] );

   }


   if( !localFlag ) DMatrix::SetMemberFlag( 0 );

   return *Temp;

}



DMatrix& DMatrix::operator ||(const DMatrix& B) const
{ // Creates a matrix by putting a matrix beside other matrix


  long i,j;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  if (n==0 || m==0) {

	Temp->Resize( B.n, B.m);

	*Temp = B;

	return *Temp;

  }

  if (n!=B.n) ERROR_MESSAGE("Wrong matrix dimensions in operator ||");



  Temp->Resize( n, m + B.m );

  for (i=1;i<=n;i++) {

    for (j=1;j<=m;j++) {

        Temp->elem(i,j) = elem( i, j );

    }

  }

  for (i=1;i<=n;i++) {

    for (j=(m+1);j<=(m+B.m);j++) {

       Temp->elem(i,j) = B.elem( i, j-m );

    }

  }

  return *Temp;

}




DMatrix& DMatrix::operator &&(const DMatrix& B) const
{ // Creates a matrix by putting a matrix below other matrix

  long i,j;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));


   if (n==0 || m==0) {

	Temp->Resize( B.n, B.m);

	*Temp = B;

	return *Temp;

  }


  if (m!=B.m) ERROR_MESSAGE("Wrong matrix dimensions in operator &&");


  Temp->Resize( n+B.n , m );

  for (i=1;i<=n;i++) {

    for (j=1;j<=m;j++) {


        Temp->elem( i, j ) = elem( i, j );

    }

  }

  for (i=(n+1);i<=(n+B.n);i++) {

    for (j=1;j<=m;j++) {

       Temp->elem( i, j ) = B.elem( i-n, j );

    }

  }

  return *Temp;

}







double enorm(const DMatrix& A)  // euclidean norm function
{

   return ( ArrayNorm( A.n*A.m, A.a ) );

}



#ifdef LAPACK
double norm(const DMatrix& A)  // matrix norm function
{
   if ( A.n > 1 || A.m > 1 )
     return  Max( SVD( A ) );
   else
     return  enorm( A );
}
#endif /* LAPACK */

double Fnorm(const DMatrix& A)  // F-norm function
{

   return sqrt( trace( TProduct(A,A) ) );

}

DMatrix& DMatrix::operator &(const DMatrix& B) const
{
    return elemProduct(*this,B);
}

DMatrix& elemProduct( const DMatrix& A, const DMatrix& B )
{

   // element by element matrix product
   DMatrix* Temp;

   Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize( A.n, A.m );

   if ( (A.n != B.n) || (A.m != B.m) ) {

      error_message( "Dimension error in elemProduct()");

   }

   for ( long i=0; i< A.n*A.m; i++ ) {


     Temp->a[i] = A.a[i] * B.a[i] ;

   }

   return *Temp;

}

DMatrix& DMatrix::operator |(const DMatrix& B) const
{
    return elemDivision(*this,B);
}


DMatrix& elemDivision( const DMatrix& A, const DMatrix& B )
{

   // element by element matrix division
   DMatrix* Temp;

   Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize( A.n, A.m );

   if ( (A.n != B.n) || (A.m != B.m) ) {

      error_message( "Dimension error in elemDivision()");

   }

   for ( long i=0; i< A.n*A.m; i++ ) {


     Temp->a[i] = A.a[i] / B.a[i] ;

   }

   return *Temp;

}




DMatrix& kronProduct( const DMatrix& A, const DMatrix& B )
{

   // kronecker product of two matrices

   DMatrix* Temp;

   Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize( A.n*B.n, A.m*B.m );

   for (long  i = 1; i<= A.n; i++ ) {

     for (long j = 1; j<= A.m; j++ ) {

            Temp->SetSubMatrix( (i-1)*B.n + 1, (j-1)*B.m + 1 ,

	    		A.elem(i,j)*B  );

     }

   }

   return (*Temp);

}

DMatrix& MatrixSign( const DMatrix& A )
{

   DMatrix* Temp;

   Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize( A.n, A.m );

   for (long i= 0; i< A.n*A.m; i++ ) {

      if ( A.a[i] > 0 ) {

         Temp->a[ i ] = 1.0;

      }
      else {

         if ( A.a[i] < 0 ) {

            Temp->a[ i ] = -1.0;

         }
         else  Temp->a[ i ] = 0.0;

      }

   }

   return *Temp;

}

double dotProduct(const DMatrix& A, const DMatrix& B)
{
// dot Product between vectors

  long i;

  double sum;

  long nelem = A.n*A.m;

  if ( !A.isVector() || !B.isVector() )
       error_message("A and B must be vectors in dotProduct()");

  if (  A.n*A.m != B.n*B.m ) error_message("dotProduct() dimension error");

  sum=0.0;


  for( i=1 ; i<= nelem; i++ )
  {

      sum += A(i)*B(i);

  }

  return sum;

}

double trace(const DMatrix& A)  // matrix trace function
{

  long i;

  double sum;

  long n = A.n;

  if (  A.n != A.m ) error_message("trace() dimension error");

  sum=0.0;


  for( i=1 ; i<= n; i++ )
  {

      sum += A(i,i);

  }

  return sum;

}


DMatrix& sum(const DMatrix& A)  // matrix sum function
{
// Returns a row vector containing the sum of values of each column

  long i, j;

  double sum;

  DMatrix* Temp;

  sum=0.0;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( 1, A.m );


  for ( j = 1; j<=A.m; j++ ) {


    for( i=1 ; i<= A.n; i++ )
    {

        sum += A.elem(i,j);

    }

    Temp->elem(1,j) = sum;

    sum=0;

  }

  return *Temp;

}

DMatrix& mean(const DMatrix& A)  // matrix mean function
{

   DMatrix* Temp;

   int localFlag = 0;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   *Temp =  ( sum( A )/( (double) A.n ) );

   if(!localFlag) DMatrix::SetMemberFlag( 0 );

   return (*Temp);

}

DMatrix& prod(const DMatrix& A)  // matrix prod function
{
// Returns a row vector containing the product of elements of each column

  long i, j;

  double sum;

  DMatrix* Temp;

  sum=1.0;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( 1, A.m );


  for ( j = 1; j<=A.m; j++ ) {


    for( i=1 ; i<= A.n; i++ )
    {

        sum *= A.elem(i,j);

    }

    Temp->elem(1,j) = sum;

    sum=1.0;

  }

  return *Temp;

}


DMatrix& Std(const DMatrix& A, int ntype)
{
// matrix standard deviation function
// Computes the standard deviation of each column of A
// The result is returned as a row vector

   int localFlag = 0;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

  long i , j ;

  long n = A.n;

  long m = A.m;

  DMatrix* Temp;

  DMatrix* ave;

  const double *Apr = A.GetConstPr();

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  ave  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( 1 ,m );

  ave->Resize( 1, m );

  Temp->FillWithZeros();

  *ave = mean( A );

  if ( m == 1 || n == 1) {

       m = MAX(m,n);

       (*Temp)(1,1) = enorm( A - ((sum(A))/m)(1,1) );

  }
  else {

      *Temp = zeros( ave->GetNoRows(), ave->GetNoCols() );

      for( i=1; i<= ave->GetNoCols(); i++ ) {

           (*Temp)(i) = 0.0;

           for(j=1; j<= A.GetNoRows(); j++) {

                (*Temp)(i) += pow( A(j,i) - (*ave)(i) , 2.0);

           }

           (*Temp)(i) = sqrt((*Temp)(i));

      }


  }

  if ( n == 1 ) Temp->FillWithZeros();

  else {

      if (ntype==0) {

         *Temp /= sqrt((double) (n-1));

      }

      else {

         *Temp /= sqrt( (double) (n) );

      }

  }

  if(!localFlag) DMatrix::SetMemberFlag( 0 );

  return *Temp;

}

DMatrix& Std(const DMatrix& A) {
      return Std(A, 0);
}

DMatrix& cov(const DMatrix& A, int ntype)
{
// Computes the covariance matrix of a data matrix where
// the N rows correspond to samples and the M columns are variables
// The result is returned as an M x M matrix. If ntype=0 (default)
// then the result is normalised with N-1. Otherwise, if ntype=1, the
// result is normalised with N.

   int localFlag = 0;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

  DMatrix::SetMemberFlag( 1 );

  long i , j , k;

  long n = A.n;

  long m = A.m;

  DMatrix* Temp;

  DMatrix* ave;

  DMatrix* As;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  As   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  ave  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( m ,m );

  ave->Resize( 1, m );

  As->Resize( n, m );

  Temp->FillWithZeros();

  *ave = mean( A );

  *As = A;

  for (j=1; j<=m; j++ ) {

    (*As)(colon(), j ) = (*As)(colon(),j) - (*ave)(j);

  }

  for (k=1; k<=n; k++) {

    for (i=1; i<= m; i++ ) {

      for (j=1; j<=m; j++ )  {

        Temp->elem(i,j) += As->elem(k,i)*As->elem(k,j);

      }

    }

  }

  if (n > 1) {

     if (ntype == 0) {

         *Temp = (*Temp)/( (double) (n-1) );

     }

     else {

         *Temp = (*Temp)/( (double) (n) );

     }

  }

  if(!localFlag) DMatrix::SetMemberFlag( 0 );

  return *Temp;

}

DMatrix& cov(const DMatrix& A) {
      return cov(A, 0);
}

DMatrix& cov( DMatrix& X, DMatrix& Y, int ntype)
{
// Covariance of two vectors
// This is equivalent to cov( [x(:) y(:)], ntype );
// See the description of cov(DMatrix& A, int ntype) above

   return cov( ( X(colon()) || Y(colon()) ), ntype );

}



DMatrix& var(DMatrix& x, int ntype)
{
// Variance of the colunms of a matrix
// Returns a row vector (a DMatrix object reference) with
// the values of the variance of each column of the input matrix A.

   return  tra(diag( cov( x, ntype ) ));

}


DMatrix& Abs(const DMatrix& A)  // matrix absolute value function
{
// Returns a matrix containing the abs values of each element of A

  long i;

  long nelem = A.n*A.m;

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( A.n , A.m );

  for ( i = 0; i< nelem; i++ ) Temp->a[i] = fabs( A.a[i] );

  return *Temp;

}

void DMatrix::SetSubMatrix( long row, long col, const DMatrix& A )
{


        if ( row+A.n-1> n || col+A.m-1> m )
        {

          error_message("Range error in SetSubMatrix()");

        }

        for( int i = 1; i<= A.n; i++) {

                for( int j=1; j<= A.m; j++) {

                    elem(row+i-1, col+j-1) = A(i,j);

                }

        }

        if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

}

DMatrix& DMatrix::sub_matrix(long r1, long r2, long c1, long c2) const
{

  long i,j;
  if(r1<1||r1>r2||r2>n||c1<1||c1>c2||c2>m) {
    printf("\n r1 = %ld", r1);
    printf("\n r2 = %ld", r2);
    printf("\n c1 = %ld", c1);
    printf("\n c2 = %ld", c2);
    printf("\n n = %ld", n);
    printf("\n m = %ld", m);
    ERROR_MESSAGE("index error in DMatrix::sub_matrix");

  }

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( r2-r1+1 , c2-c1+1 );

  for (i=1; i<=r2-r1+1; i++) {

    for (j=1; j<=c2-c1+1; j++) {

      Temp->elem(i,j)=elem(i+r1-1,j+c1-1);

    }

  }

  return *Temp;

}




double& DMatrix::operator() (long row, long col )
{


  if (row > n || col > m )

	  DMatrix::Resize(row,col);

  if (row < 0 || col < 0)

    ERROR_MESSAGE("negative index error in A(row,col)");

  return elem( row, col );

}

double& DMatrix::operator() (long row, const char* end )
{

  long col = this->getm();

  if (row > n || col > m )

	  DMatrix::Resize(row,col);

  if (row < 0 || col < 0)

    ERROR_MESSAGE("negative index error in A(row,col)");

    return elem( row, col );

}

double& DMatrix::operator() (const char* end, long col )
{

  long row = this->getn();

  if (row > n || col > m )

	  DMatrix::Resize(row,col);

  if (row < 0 || col < 0)

    ERROR_MESSAGE("negative index error in A(row,col)");

    return elem( row, col );

}



double DMatrix::operator() (long row, long col ) const
{

  if (row > n || row < 1 || col > m || col < 1)

    ERROR_MESSAGE("index error in A(row,col)");

    return elem( row, col );

}

double DMatrix::operator() (long row,const char* end ) const
{

  long col = this->getm();

  if (row > n || row < 1 || col > m || col < 1)

    ERROR_MESSAGE("index error in A(row,col)");

  return elem( row, col );

}

double DMatrix::operator() (const char* end, long col ) const
{

  long row = this->getn();

  if (row > n || row < 1 || col > m || col < 1)

    ERROR_MESSAGE("index error in A(row,col)");

    return elem( row, col );

}



double DMatrix::element(long row, long col)
{

  if (row > n || row < 1 || col > m || col < 1)

    ERROR_MESSAGE("index error in DMatrix::element(row,col)");

    return elem( row, col );

}


double& DMatrix::operator() (long index )
{


  if (index > asize && (n==1 || m==1) )
  {
	  if (n==1)
	  {
	     	DMatrix::Resize(1,index);
	  }
	  else if (m==1)
	  {
		DMatrix::Resize(index,1);
	  }
  }

  if ( (n==0 || m==0) )
  {
	  DMatrix::Resize(index,1);
  }

  if (index > asize && (n>1 && m>1) )
  {
     ERROR_MESSAGE("Non-vector matrix cannot be resized in call to A(I) in LHS");
  }

  if (index > asize || index < 1 ) ERROR_MESSAGE("Index error in operator()(long)");


  return a[ index - 1 ];


}


double& DMatrix::operator() (const char* end )
{

  long index = (this->getn())*(this->getm());

  if (index > asize && (n==1 || m==1) )
  {
	  if (n==1)
	  {
	     DMatrix::Resize(1,index);
	  }
	  else if (m==1)
	  {
				DMatrix::Resize(index,1);
	  }
  }

  if ( (n==0 || m==0) )
  {
	  DMatrix::Resize(index,1);
  }

  if (index > asize && (n>1 && m>1) )
  {
     ERROR_MESSAGE("Non-vector matrix cannot be resized in call to A(I) in LHS");
  }

  if (index > asize || index < 1 ) ERROR_MESSAGE("Index error in operator()(long)");

  return a[ index - 1 ];


}


double DMatrix::operator() (const char* end) const
{

  long index = (this->getn())*(this->getm());

  if (index > asize || index < 1 ) ERROR_MESSAGE("Index error in operator()(long)");

  return a[ index - 1 ];

}

double DMatrix::operator() (long index ) const
{

  if (index > asize || index < 1 ) ERROR_MESSAGE("Index error in operator()(long)");

  return a[ index - 1 ];

}




DMatrix::~DMatrix()
{

  if (atype==0 && a!=NULL )
  {

      DMatrix::DeAllocate();

  }


}


void DMatrix::Save(const char * FileName )
{

  long i,j;

  FILE *fp;

  if ( (fp = fopen(FileName,"w")) == NULL )

  {  ERROR_MESSAGE( "Error opening file in DMatrix::Save()"); }


    for (i=1;i<=n;i++) { for (j=1;j<=m;j++) {

    fprintf(fp,"%g\t",elem(i,j)); } fprintf(fp,"\n"); }

  fclose(fp);

}

void DMatrix::Load(const char * FileName )
{

  long i,j;

  FILE *fp;

  if ( (fp = fopen(FileName,"r")) == NULL )

  {  ERROR_MESSAGE( "Error opening file in DMatrix::Load()"); }


    for (i=1;i<=n;i++) { for (j=1;j<=m;j++) {

    fscanf(fp,"%lf", &elem(i,j)); } }


  fclose(fp);


}


void DMatrix::Fprint( FILE *fp )
{

  long i,j;

    for (i=1;i<=n;i++) { for (j=1;j<=m;j++) {

    fprintf(fp,"%f\t",elem(i,j)); } fprintf(fp,"\n"); }


}

double Max(const DMatrix & A) {
      return Max(A, nullptr, nullptr);
}
                    
double Max(const DMatrix & A, int* rindx, int* cindx )
{

   double mx;
   const double *Apr = A.GetConstPr();

   int  ri, ci;

   long inx;

   mx = Apr[0];

   inx = 0;

   for (long i=1; i<A.n*A.m; i++)
   {

     if (Apr[i] > mx) { mx = Apr[i]; inx = i; }

   }

   ci = (inx)/A.n + 1;

   ri = (inx+1) - (ci-1)*A.n;


   if ( rindx != NULL ) *rindx = ri;
   if ( cindx != NULL ) *cindx = ci;


   return mx;

}

double InfNorm(const DMatrix & A )
{

  // Infinity norm of matrix
  return Max( sum( Abs(tra(A)) ) );

}



double Min(const DMatrix & A, int* rindx, int* cindx )
{

   double mx;

   int ri, ci;

   long inx;


   const double *Apr = A.GetConstPr();

   mx = Apr[0];

   inx = 0;

   for (long i=1; i<A.n*A.m; i++)
   {

     if (Apr[i] < mx) { mx = Apr[i]; inx = i; }

   }

   ci = (inx)/A.n + 1;

   ri = (inx+1) - (ci-1)*A.n;

   if ( rindx != NULL ) *rindx = ri;

   if ( cindx != NULL ) *cindx = ci;

   return mx;

}

double MaxAbs(const DMatrix& A) {
      return MaxAbs(A, nullptr, nullptr);
}
                    
double MaxAbs(const DMatrix & A, int* rindx, int* cindx )
{

   double mx;

   int ri, ci;

   long inx;

   const double *Apr = A.GetConstPr();

   mx = fabs(Apr[0]);

   inx = 0;

   for (long i=1; i<A.n*A.m; i++)
   {

     if (fabs(Apr[i]) > mx) { mx = fabs(Apr[i]); inx = i; }

   }

   ci = (inx)/A.n + 1;

   ri = (inx+1) - (ci-1)*A.n;

   if ( rindx != NULL ) *rindx = ri;

   if ( cindx != NULL ) *cindx = ci;

   return mx;

}

double MinAbs(const DMatrix & A, int* rindx, int* cindx )
{

   double mx;

   int ri, ci;

   long inx;

   const double *Apr = A.GetConstPr();

   mx = fabs(Apr[0]);

   inx = 0;

   for (long i=1; i<A.n*A.m; i++)
   {

     if (fabs(Apr[i]) < mx) { mx = fabs(Apr[i]), inx=i; }

   }

   ci = (inx)/A.n + 1;

   ri = (inx+1) - (ci-1)*A.n;

   if ( rindx != NULL ) *rindx = ri;

   if ( cindx != NULL ) *cindx = ci;

   return mx;

}


void  sort( DMatrix & A, int indx[]  )
{
   /* Sorts vector A in ascending order              */
   /* Returns the array of sorted indices (optional) */

   double *a;
   long   n = MAX(A.n, A.m);

   if ( A.n!=1 && A.m != 1 ) {

      error_message("Argument must be a column or row vector in sort() ");

   }


   unsigned long i,j,inc;
   int vi;
   double v;
   inc=1;

   a = A.a - 1;

   if (indx!= NULL) { for( i=1; i<= n; i++ )  indx[i] = i; }


        do {
                inc *= 3;
                inc++;
        } while (inc <= n);
        do {
                inc /= 3;
                for (i=inc+1;i<=n;i++) {
                        v=a[i];
                        if (indx!=NULL) vi=indx[i];
                        j=i;
                        while (a[j-inc] > v) {
                                a[j]=a[j-inc];
                                if ( indx!=NULL ) indx[j]=indx[j-inc];
                                j -= inc;
                                if (j <= inc) break;
                        }
                        a[j]=v;
                        if ( indx!=NULL ) indx[j]=vi;

                }
        } while (inc > 1);


}



void  sort( DMatrix & A, DMatrix& indx  )
{
   /* Sorts vector A in ascending order              */
   /* Returns the array of sorted indices  */

   double *a;
   long   n = MAX(A.n, A.m);

   if ( A.n!=1 && A.m != 1 ) {

      error_message("Argument must be a column or row vector in sort() ");

   }


   unsigned long i,j,inc;
   int vi;
   double v;
   inc=1;

   a = A.a - 1;

   indx.Resize(A.n,A.m);

   if (1) { for( i=1; i<= n; i++ )  indx(i) = i; }


        do {
                inc *= 3;
                inc++;
        } while (inc <= n);
        do {
                inc /= 3;
                for (i=inc+1;i<=n;i++) {
                        v=a[i];
                        vi=(int) indx(i);
                        j=i;
                        while (a[j-inc] > v) {
                                a[j]=a[j-inc];
                                indx(j)=indx(j-inc);
                                j -= inc;
                                if (j <= inc) break;
                        }
                        a[j]=v;
                        indx(j)=vi;

                }
        } while (inc > 1);


}


DMatrix& DMatrix::operator^(double x)
{
//  Elementwise matrix power operator A^p
//  Example:  C = A^2;


  int i;

  int localFlag = 0;

  if (DMatrix::GetMemberFlag()) localFlag = 1;

  DMatrix::SetMemberFlag( 1 );

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  Temp->FillWithZeros();

  for (i=0; i<n*m; i++) {

    Temp->a[i] = pow(this->a[i], x);

  }

  if(!localFlag) DMatrix::SetMemberFlag( 0 );

  return *Temp;

}

DMatrix& DMatrix::mpow(int p)
{
//  Integer matrix power function A^p
//  Example:  C = A.mpow(2);

  int i;

  int localFlag = 0;

  if (DMatrix::GetMemberFlag()) localFlag = 1;

  DMatrix::SetMemberFlag( 1 );

  DMatrix inverse;

  if (n!=m) ERROR_MESSAGE("Matrix must be square in mpow()");

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  Temp->FillWithZeros();

  for (i=1; i<=n; i++) {

    Temp->elem(i,i) = 1.0;

  }

  if (p<0) {
      inverse = inv(*this);
  }

  for (i=1; i<=abs(p); i++) {

     if ( p>0 ) {
     	*Temp = (*Temp)*(*this);
     }
     else if (p<0) {
	*Temp = (*Temp)*inverse;
     }

  }

  if(!localFlag) DMatrix::SetMemberFlag( 0 );

  return *Temp;

}

DMatrix& mpow(DMatrix& A, int p)
{
//  Integer power function A^p, friend implementation
//  Example:  C = mpow(A,2);
   return A.mpow( p );

}

void DMatrix::Allocate(long size)
{
        /* Memory allocation function */



        a = (double *) my_calloc( size, sizeof(double) );

        if (!a) ERROR_MESSAGE("allocation failure in DMatrix::Allocate()");

        // a = a-nl+NR_END;

        asize = size;

        allocated = true;

}



// extern "C"   char *gcvt( double, int, char*);


char* num2str(double num)
{


  static char str[20];
#ifndef MEX_OR_XPC
  int sig = 10;
  gcvt( num, sig, str);
#endif

  return str;

}

/*  ***************************************** */
/*                                            */
/*  DATA INTERFACE BETWEEN MATLAB AND DMATRIX */
/*                                            */
/*  ***************************************** */

#ifdef MATLAB_MEX_FILE


void mxArray2DMatrix( DMatrix& A, mxArray* mx )
{

   long nrow, ncol;

   nrow = mxGetM( mx );

   ncol = mxGetN( mx );

   A.Resize( nrow, ncol );

   memcpy( A.geta(), mxGetPr(mx), nrow*ncol*sizeof(double) );

}

void DMatrix2mxArray(mxArray* mx, DMatrix& A)
{


   mx = mxCreateDoubleMatrix( A.getn(), A.getm(), mxREAL );

   memcpy( mxGetPr(mx), A.geta(), A.getn()*A.getm()*sizeof(double) );

}

#endif


#define TINY 1.0e-20;

void LU_Decomposition(DMatrix& a, int n, DMatrix& indx, double *d, DMatrix& vv)
{
        int i,imax,j,k,l;
        double mc,dum,sum;

        *d=1.0;

        for (i=1;i<=n;i++) {

                mc = 0.0;

                for(l=1;l<=n;l++) {

                   if( fabs( a(l,i)) > mc ) mc = fabs(a(l,i));

                }

                if (mc == 0.0) error_message("Singular matrix in routine LU_Decomposition");

                vv(i)=1.0/mc;

        }

        for (j=1;j<=n;j++) {

                for (i=1;i<j;i++) {

                        sum=a(i,j);

                        for (k=1;k<i;k++) sum -= a(i,k)*a(k,j);

                        a(i,j)=sum;

                }

                mc=0.0;

                for (i=j;i<=n;i++) {

                        sum=a(i,j);

                        for (k=1;k<j;k++)

                                sum -= a(i,k)*a(k,j);

                        a(i,j)=sum;

                        if ( (dum=vv(i)*fabs(sum)) >= mc ) {

                                mc=dum;

                                imax=i;

                        }

                }

                if (j != imax) {

                        a.SwapRows(j,imax);

                        *d = -(*d);

                        vv(imax)=vv(j);

                }

                indx(j)=imax;

                if (a(j,j) == 0.0) a(j,j)=TINY;

                if (j != n) {

                        dum=1.0/(a(j,j));

                        for (i=j+1;i<=n;i++) a(i,j) *= dum;

                }

        }

}
#undef TINY




void LU_Back_Substitution(const DMatrix& a, int n,
                          const DMatrix& indx, DMatrix& b)
{

        int i,ii=0,ip,j;
        double sum;

        for (i=1;i<=n;i++) {

                ip= (int) indx(i);

                sum=b(ip);

                b(ip)=b(i);

                if (ii)

                        for (j=ii;j<=i-1;j++) sum -= a(i,j)*b(j);

                else if (sum) ii=i;

                b(i)=sum;

        }

        for (i=n;i>=1;i--) {

                sum=b(i);

                for (j=i+1;j<=n;j++) sum -= a(i,j)*b(j);

                b(i)=sum/a(i,i);

        }
}


DMatrix& LUSolve( const DMatrix& A, const DMatrix &b )
{

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix *ALU;
     DMatrix *xx;
     DMatrix *vv;
     DMatrix *indx;
     DMatrix *rr;

     int nn = A.GetNoRows();
     double dd;
     int i,j;

     if ( nn !=A.GetNoCols() ) {

        error_message("\n Matrix must be square in LUSolve ");

     }

     ALU = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     xx  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     vv  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     indx= DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     rr  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     ALU->Resize( nn , A.GetNoCols() );
     xx-> Resize( nn , 1 );
     vv-> Resize( nn , 1 );
     indx->Resize( nn, 1 );
     rr-> Resize( nn , b.GetNoCols() );


     *ALU = A;

     LU_Decomposition(*ALU, nn, *indx, &dd, *vv);

     for ( i=1; i<= b.GetNoCols(); i++ ) {

       for( j=1; j<= nn; j++ )  (*xx)(j) = b(j,i);

       LU_Back_Substitution(*ALU, nn, *indx, *xx );

       for( j=1; j<= nn; j++ )  (*rr)(j,i) = (*xx)(j);

     }

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *rr;

}

DMatrix& LUFSolve( const DMatrix& ALU, const DMatrix &b )
{

     /* Solution to A x = b by LU factorization */
     /* Alu obtained from ALU = LU( A ) */
     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix *xx;
     DMatrix *indx;
     DMatrix *rr;

     int nn = ALU.GetNoRows();

     int i;
     int j;

     if ( nn != ALU.GetNoCols()-1 ) {

        error_message("\n Dimension Error in LUFSolve ");

     }

     xx  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     indx= DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     rr  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));


     xx->  Resize( nn , 1 );
     rr->  Resize( nn, b.GetNoCols() );

     for ( i = 1; i<= nn; i++ ) {

            (*indx)( i ) =  ALU( i, ALU.GetNoCols() );

     }


     for ( i=1; i<= b.GetNoCols(); i++ ) {

        for( j=1; j<=nn; j++ ) (*xx)(j) = b(j,i);

        LU_Back_Substitution( ALU, nn, *indx, *xx );

        for( j=1; j<=nn; j++ ) (*rr)(j,i) = (*xx)(j);

     }

     if( !localFlag ) DMatrix::SetMemberFlag( 0 );

     return *rr;

}

double det( const DMatrix& A )
{
     /* Determinant of a matrix */
     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );


     DMatrix *ALU;
     DMatrix *vv;
     DMatrix *indx;

     int nn = A.GetNoRows();
     double dd;
     int i;

     ALU = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     vv  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     indx= DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     ALU->Resize( nn , A.GetNoCols() );
     vv-> Resize( nn , 1 );
     indx->Resize( nn, 1 );

     if ( nn !=A.GetNoCols() ) {

        error_message("\n Matrix must be square in det() ");

     }


     *ALU = A;


     LU_Decomposition(*ALU, nn, *indx, &dd, *vv);

     for ( i = 1; i<=  nn; i++ ) {

       dd*= (*ALU)(i,i);

     }

     if( !localFlag ) DMatrix::SetMemberFlag( 0 );

     return dd;

}


DMatrix& LU( const DMatrix& A  )
{

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );


     DMatrix *ALU;
     DMatrix *vv;
     DMatrix* indx;

     int nn = A.GetNoRows();
     int i,j;
     double dd;

     ALU = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     vv  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     indx= DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     ALU->Resize( nn , A.GetNoCols() + 1 );
     vv-> Resize( nn , 1 );
     indx->Resize( nn, 1 );

     if ( nn !=A.GetNoCols() ) {

        error_message("\n Matrix must be square in LU() ");

     }


     for ( i=1; i<= nn; i++ ) {

       for ( j= 1; j<=nn; j++ ) {

         (*ALU)(i,j) = A(i,j);

       }

     }

     LU_Decomposition(*ALU, nn, *indx, &dd, *vv);

     for ( i = 1; i<= nn; i++ ) {
     /* Stores indx vector in last column of ALU */
         (*ALU)( i, nn + 1 )=  (*indx)(i);
     }

     if( !localFlag ) DMatrix::SetMemberFlag( 0 );

     return *ALU;

}


DMatrix& DMatrix::operator % (const DMatrix& rval) const
{

  // Left division implementation, B%A == inv(B)*A

  int localFlag = 0;

  if (DMatrix::GetMemberFlag()) localFlag = 1;

  DMatrix::SetMemberFlag( 1 );

  DMatrix *Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( m, rval.m );

if ( n == m ) {

  int i, j;

  double dd;

  DMatrix* Alu;

  DMatrix* indx;

  DMatrix* vv;

  DMatrix* xx;


  Alu  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  indx = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  vv   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  xx   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Alu ->Resize( n,  m );

  indx->Resize( n, 1);

  vv->Resize( n, 1 );

  xx->Resize( n, 1 );

  for( i=0; i< n*n; i++ ) Alu->GetPr()[i] = this->GetConstPr()[i];

  LU_Decomposition(*Alu, n, *indx, &dd, *vv);


  for( i = 1; i<= rval.m; i++ ) {

     for ( j=1; j<=n ; j++ ) (*xx)(j) = rval(j,i);

     LU_Back_Substitution(*Alu, n, *indx, *xx );

     for ( j=1; j<=n; j++ ) (*Temp)(j,i) = (*xx)(j);

  }

}

else {

#ifdef LAPACK

 *Temp = LSMNSolve( *this, rval );

#else

  error_message(" Dimension error in DMatrix:: operator %");

#endif

}

if(!localFlag) DMatrix::SetMemberFlag( 0 );

return *Temp;

}


void CholeskyDecomp(DMatrix& a, int n, DMatrix& pM)
{
        int i,j,k;

        double sum;

        double *p = pM.GetPr()-1;

        for (i=1;i<=n;i++) {

                for (j=i;j<=n;j++) {

                        sum  = a.elem(i,j);

                        for (k=i-1;k>=1;k--) {

                           sum -= a.elem(i,k)*a.elem(j,k);

                        }

                        if (i == j) {

                                if (sum <= 0.0)

                                        error_message("CholeskyDecomp() failed");

                                p[i]=sqrt(sum);

                        } else a.elem(j,i)=sum/p[i];

                }

        }

}


void CholeskySolution(const DMatrix& a, int n, const DMatrix& pM, const DMatrix& bM, DMatrix& xM)
{
        int i,k;

        double sum;

        const double *p = pM.GetConstPr()-1;

        const double *b = bM.GetConstPr()-1;

        double *x = xM.GetPr()-1;

        for (i=1;i<=n;i++) {

                sum = b[i];

                for (k=i-1;k>=1;k--) {

                   sum -= a.elem(i,k)*x[k];

                }

                x[i]=sum/p[i];

        }

        for (i=n;i>=1;i--) {

                  sum = x[i];

                for (k=i+1;k<=n;k++) {


                  sum -= a.elem(k,i)*x[k];

                }

                x[i]=sum/p[i];

        }

}


DMatrix& CholSolve( const DMatrix& A, const DMatrix &b )
{

     /* Ax = b solution using Cholesky factorization  */
     /* A must be a positive definite symmetric matrix */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix *Achol;
     DMatrix *xx;
     DMatrix *pp;
     DMatrix *rr;
     DMatrix *bb;

     int nn = A.GetNoRows();

     int i,j;

     Achol = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     xx    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     pp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     rr    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     bb    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     Achol->Resize( nn , A.GetNoCols() );
     xx   ->Resize( nn , 1             );
     pp   ->Resize( nn , 1             );
     rr   ->Resize( nn , b.GetNoCols() );
     bb   ->Resize( nn , 1 );

     if ( nn !=A.GetNoCols() ) {

        error_message("\n Matrix must be square in CholSolve() ");

     }

     *Achol = A;

     CholeskyDecomp( *Achol , nn, *pp );

     for ( i=1; i<= b.GetNoCols(); i++) {

        for( j=1; j<= nn; j++ )  (*bb)(j) = b(j,i);

        CholeskySolution( *Achol,  nn, *pp, *bb, *xx);

        for( j=1; j<= nn; j++ )  (*rr)(j,i) = (*xx)(j);

     }

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *rr;

}

DMatrix& Chol( const DMatrix& A  )
{

     /* Cholesky factorization  */
     /* A must be a positive definite symmetric matrix */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix *Achol;
     DMatrix *pp;

     int nn = A.GetNoRows();

     int i, j;

     Achol = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     pp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     Achol->Resize( nn , A.GetNoCols()+1 );
     pp   ->Resize( nn , 1             );

     if ( nn !=A.GetNoCols() ) {

        error_message("\n Matrix must be square in Chol() ");

     }

     for ( i = 1; i<= nn; i++ ) {

       for ( j =1 ; j<= nn; j++ ) {

          (*Achol)(i,j) = A(i,j);

       }

     }

     CholeskyDecomp( *Achol , nn, *pp );

     for ( i = 1; i<= nn; i++ ) {

        (*Achol)( i, nn + 1 ) = (*pp)(i);

     }

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *Achol;

}

DMatrix& CholeskyRoot( const DMatrix& A  )
{

     /* Cholesky Root  */
     /* A must be a positive definite symmetric matrix */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix *Achol;
     DMatrix *Acholrtn;
     DMatrix *pp;

     int nn = A.GetNoRows();

     int i, j;

     Achol    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     Acholrtn = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     pp       = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     Achol->Resize( nn , A.GetNoCols()+1 );
     Acholrtn->Resize( nn, nn );
     pp   ->Resize( nn , 1             );

     if ( nn !=A.GetNoCols() ) {

        error_message("\n Matrix must be square in Chol() ");

     }

     for ( i = 1; i<= nn; i++ ) {

       for ( j =1 ; j<= nn; j++ ) {

          (*Achol)(i,j) = A(i,j);

       }

     }

     CholeskyDecomp( *Achol , nn, *pp );


     Acholrtn->FillWithZeros();

     for(i = 1; i<= nn; i++ ) {

       for( j = 1; j< i; j++ ) {

          (*Acholrtn)(i,j) = (*Achol)(i,j);

       }


     }

     for(i=1; i<= nn; i++) (*Acholrtn)(i,i) = (*pp)(i);

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *Acholrtn;

}





DMatrix& CholFSolve( const DMatrix& Achol, const DMatrix &b )
{

     /* Ax = b solution using Cholesky factorization  */
     /* A must be a positive definite symmetric matrix */
     /* Achol obtained from Achol = CHOL( A ); */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix *xx;
     DMatrix *pp;
     DMatrix *rr;
     DMatrix *bb;

     int nn = Achol.GetNoRows();

     int i,j;

     xx    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     pp    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     rr    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     bb    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     xx   ->Resize( nn , 1             );
     pp   ->Resize( nn , 1             );
     rr   ->Resize( nn , b.GetNoCols() );
     bb   ->Resize( nn , 1             );

     if ( nn !=Achol.GetNoCols()-1 ) {

        error_message("\n Dimension error in CholFSolve() ");

     }

     /* Extract pp from last column of Achol */

     for (i=1; i<=nn; i++ ) {

      (*pp)(i) = Achol( i, nn+1 );

     }

     for (i=1; i<= b.GetNoCols(); i++ ) {

        for ( j=1; j<=nn; j++ )  (*bb)(j) = b(j,i);

        CholeskySolution( Achol,  nn, *pp, *bb, *xx);

        for ( j=1; j<=nn; j++ )  (*rr)(j,i) = (*xx)(j);

     }

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *rr;

}



void Householder(int j, int k,double *aux,double *d)
{
   double s;

   if(j > k) error_message("Dim error in Householder");

   if(j == k) s = fabs(aux[j]);

   else  s = ArrayNorm(k-j+1,aux+j);

   if(s==0.)
   {
      *d=0.;
      for(int i=j;i<=k;i++) aux[i]=0.;
      return;
   }
   double h,w;

   w = sqrt((1. + fabs(aux[j])/s)*.5);

   if(aux[j] >= 0.)
   { h = 2.*s*w; *d = -s; }
   else
   { h = -2.*s*w; *d = s;}

   aux[j] = (w);

   for(int i = j+1;i <= k;i++)

      aux[i] = (aux[i]/h);

}

void HouseholderApplyLeft(int ri, int rf,int ci,int cf,
           DMatrix& a,double *aux)
{
   int k;
   for(int l = ci;l <= cf;l++)
      {
      double s = 0.;
      for(k = ri;k <= rf;k++)s += aux[k]*a.elem(k,l);
      s += s;
      for(k = ri;k <= rf;k++)
        {
        double *temp = &a.elem(k,l);
        *temp = (*temp - s*aux[k]);
        }
      }
}

char QRFactorization(int m,int n,DMatrix& a,double *d,
                     double *aux, double* nrm)
{
   if(m < n)
      error_message("nrows < ncols in QRFactorization");
      const double SIN_TINY = 10.*DMatrix::GetEPS();
   int i,j;
   for(j = 1;j <= n;j++)
      {
      for(i = 1;i <= m;i++)
        aux[i] = a.elem(i,j);
      nrm[j] = ArrayNorm(n,aux+1);
      }
   char sing = 0;
   for(j = 1;j <= n;j++)
      {
      for(i = j;i <= m;i++)
        aux[i] = a.elem(i,j);
      Householder(j,m,aux,&d[j]);
      if(fabs(d[j]) <= SIN_TINY*nrm[j] || nrm[j] == 0.)
        {
        sing = 1;
        d[j] = 0.;
        error_message("\nSingular matrix in QRDec");
        for(i = j;i <= m;i++)
           { aux[i] = a.elem(i,j) = 0.;}
        aux[j] = a.elem(i,j) = 1.;
        }
      else
        for(i = j;i <= m;i++)
           a.elem(i,j)=aux[i];
      HouseholderApplyLeft(j,m,j+1,n,a,aux);
      }
   return sing;
}

void QRSolution(int m,int n,const DMatrix& a,
              double *d,double *b,double *x)
{
   int j,k,i;
   double s;
   if(m < n) error_message("nrows < ncols in QRSolution");
   for(j=1;j<=n;j++)
      {
      s=0.;
      for(k = j;k <= m;k++)
        s -= a.elem(k,j)*b[k];
      s += s;
      for(k = j;k <= m;k++)
        b[k] = (b[k] + s*a.elem(k,j));
      }
   for(i = n;i >= 1;i--)
      {
      if(d[i]==0.)
        {
        error_message("\nSingular matrix in QRSolution!");
        x[i]=0.;
        }
      else
        {
        s = b[i];
        for(k = i+1;k <= n;k++)
           s -= a.elem(i,k)*x[k];
        x[i] = (s/d[i]);
        }
      }
}


DMatrix& QRSolve( const DMatrix& A, const DMatrix &b )
{

     /* Ax = b solution using QR decomposition  */
     /* A must have nRows >= nCols */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     int i;

     DMatrix *Aqr;
     DMatrix *xx;
     DMatrix *aux;
     DMatrix *nrm;
     DMatrix *dd;
     DMatrix *bb;
     DMatrix *xout;

     int nn = A.GetNoRows();
     int mm = A.GetNoCols();


     Aqr   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     xx    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     aux   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     dd    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     nrm   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     bb    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     xout  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     Aqr  ->Resize( nn   , A.GetNoCols() );
     xx   ->Resize( mm   , 1             );
     aux  ->Resize( nn+1 , 1             );
     dd   ->Resize( mm , 1  );
     nrm  ->Resize( mm + 1 , 1  );
     bb   ->Resize( nn  , 1  );
     xout ->Resize( A.m, b.m );


     memcpy(Aqr->a, A.a, A.n*A.m*sizeof(double));

     QRFactorization(A.GetNoRows(),A.GetNoCols(),*Aqr,dd->GetPr()-1,
                     aux->GetPr()-1, nrm->GetPr()-1 );

     for( i=1; i<=b.m; i++ ) {

       memcpy( bb->a, b.a+(i-1)*b.n, b.n*sizeof(double) );

       QRSolution(A.GetNoRows(),A.GetNoCols(),*Aqr,
               dd->GetPr()-1,bb->GetPr()-1,xx->GetPr()-1 );

       memcpy( xout->a + (i-1)*xout->n, xx->a, xout->n*sizeof(double) );

     }


     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *xout;

}

DMatrix& QRFSolve( const DMatrix& Aqr, const DMatrix &b )
{

     /* Ax = b solution using QR decomposition  */
     /* Aqr is obtained from Aqr = QR( A )      */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     int i;

     DMatrix *xx;
     DMatrix *dd;
     DMatrix *bb;
     DMatrix *xout;

     int nn = Aqr.GetNoRows();
     int mm = Aqr.GetNoCols()  - 1;


     xx    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     dd    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     bb    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     xout  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     xx   ->Resize( mm   , 1           );
     dd   ->Resize( mm+1 , 1           );
     bb   ->Resize( nn , 1             );
     xout ->Resize( mm , b.m           );


     /* Extract dd from Aqr */

    for ( i=1; i<=mm; i++ ) {

     (*dd)(i) = Aqr( i, mm + 1 );

    }


    for( i=1; i<=b.m; i++ ) {

       memcpy( bb->a, b.a+(i-1)*b.n, b.n*sizeof(double) );

       QRSolution( nn , mm,  Aqr,
               dd->GetPr()-1,bb->GetPr()-1,xx->GetPr()-1 );

       memcpy( xout->a + (i-1)*xout->n, xx->a, xout->n*sizeof(double) );

    }


    if(!localFlag) DMatrix::SetMemberFlag( 0 );

    return *xout;

}


DMatrix& QR( const DMatrix& A  )
{

     /*  QR decomposition  */
     /* A must have nRows>nCols  */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix *Aqr;
     DMatrix *aux;
     DMatrix *nrm;
     DMatrix *dd;

     int i,j;

     int nn = A.GetNoRows();
     int mm = A.GetNoCols();


     Aqr   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     aux   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     nrm   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     dd    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     Aqr   ->Resize( nn , mm + 1 );
     aux   ->Resize( nn +1 , 1             );
     nrm   ->Resize( mm +1 , 1             );
     dd    ->Resize( mm , 1             );

     for ( i = 1; i<= nn; i ++ ) {

        for ( j = 1; j<= mm; j++ ) {

         (*Aqr)( i, j )  = A( i, j );

        }

     }

     QRFactorization(A.GetNoRows(),A.GetNoCols(),*Aqr,dd->GetPr()-1,
                     aux->GetPr()-1, nrm->GetPr()-1 );

     for ( i = 1; i<= mm; i ++ ) {

        (*Aqr)( i, mm + 1 ) = (*dd)(i);

     }

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *Aqr;

}


#ifdef LAPACK

DMatrix& LQ( const DMatrix& A, DMatrix* Q  )
{

     /*  LQ decomposition  */


     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix *Alq;
     DMatrix *wk;
     DMatrix *tau;
     DMatrix *H;
     DMatrix *v;


     int i,j;

     integer nn = A.GetNoRows();
     integer mm = A.GetNoCols();

     integer nt = min( nn, mm );

     integer nw = 10 * nn;

     integer info;


     Alq   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     wk    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     tau   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     H     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     v     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));


     Alq   ->Resize( nn , mm );
     wk    ->Resize( 10*nn , 1             );
     tau   ->Resize( nt, 1 );
     Q     ->Resize( mm, mm );
     H     ->Resize( mm, mm );
     v     ->Resize( mm, 1 );

     for ( i = 1; i<= nn; i ++ ) {

        for ( j = 1; j<= mm; j++ ) {

           (*Alq)( i, j )  = A( i, j );


        }

     }


     /* int dgelqf_(integer *m, integer *n, doublereal *a, integer *
        lda, doublereal *tau, doublereal *work, integer *lwork, integer *info)
      */


     dgelqf_(&nn, &mm, Alq->GetPr(), &nn, tau->GetPr(), wk->GetPr(), &nw, &info );


     *Q = identity(mm);


     for(i=1; i<=nn; i++ ) {

        v->FillWithZeros();

        (*v)(i) = 1.0;

        for(j=i+1; j<=mm; j++ ) {
                      (*v)(j) = (*Alq)(i,j);
                      (*Alq)(i,j) = 0;
        }

        *H  = ProductT(*v,*v);
        *H *= -(*tau)(i);
        *H += identity(mm);

        (*Q) *= (*H);

     }

     *Q = -tra(*Q);

     *Alq = - (*Alq);

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *Alq;

}

#endif




#ifdef LAPACK

DMatrix& SVD( const DMatrix& A, DMatrix* U, DMatrix* V  )
{

      // This function computes the singular values of A
      // It calls dgesvd_(), which is part of CLAPACK, a
      // public domain numerical linear algebra library

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     char JOBU;
     char JOBVT;

     integer M = A.GetNoRows();
     integer N = A.GetNoCols();

     integer LDA = M;
     integer LDU = M;
     integer LDVT = N;
     integer LWORK;
     integer INFO;

     int nrows = A.GetNoRows();
     int ncols = A.GetNoCols();

     double* upr;

     double* vtpr;

     DMatrix * a;
     DMatrix * s;
     DMatrix * wk;


     s     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     a     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     wk    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     s   ->Resize( MIN( nrows, ncols) , 1             );
     a   ->Resize( nrows, ncols   );


     if ( U == NULL ) {

       JOBU  = 'N';
       upr   = NULL;

     }

     else {

       JOBU = 'A';

       U->Resize( LDU , M );

       upr = U->GetPr();

     }

     if ( V == NULL ) {

       JOBVT = 'N';

       vtpr  = NULL;

     }

     else {

       JOBVT = 'A';
       V->Resize( LDVT , N  );
       vtpr = V->GetPr();

     }

     memcpy( a->a, A.a, A.n*A.m*sizeof(double) );


     LWORK = DMatrix::GetDimAux();


/* Subroutine int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
        doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
        ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
        integer *info)
*/

     dgesvd_( &JOBU, &JOBVT, &M, &N,
          a->GetPr(), &LDA, s->GetPr(), upr,
          &LDU, vtpr , &LDVT, wk->GetPr(), &LWORK, &INFO);

     if ( V != NULL ) V->Transpose();

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *s;


}

#endif /* LAPACK */

#ifdef LAPACK
double cond( const DMatrix& A )
{
     /* Condition number of a matrix */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix* ww;

     double sigmaMax;

     double sigmaMin;


     ww     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     ww->Resize( MIN( A.n, A.m ) , 1 );

     *ww = SVD( A );

     sigmaMax = Max( *ww );

     sigmaMin = Min( *ww );

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return (sigmaMax/sigmaMin);

}
#endif /* LAPACK */

#ifdef  LAPACK
int rank_matrix( const DMatrix& A )
{
     /* rank of a matrix */

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix* ww;

     double tol;
     double eps = DMatrix::GetEPS();

     int i, sum;

     ww     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     ww->Resize( MIN( A.GetNoCols(), A.GetNoRows() ) , 1 );

     *ww = SVD( A );

     sum = 0;

     tol = MAX( A.n, A.m )* Max( *ww ) * eps;

     for ( i = 1; i <= MIN(A.n, A.m) ; i++)

         sum += (*ww)(i) > tol ? (1) : (0) ;

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return  sum;

}
#endif /* LAPACK */

#ifdef LAPACK
DMatrix& pinv(const DMatrix& A)
{

     int localFlag = 0;

     double tol;

     double eps = DMatrix::GetEPS();

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );


     int nrows = A.GetNoRows();
     int ncols = A.GetNoCols();

     DMatrix * ww;
     DMatrix * uu;
     DMatrix * vv;
     DMatrix * dd;
     DMatrix * rr;

     ww     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     uu     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     vv     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     dd     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     rr     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     ww   ->Resize( MIN(nrows,ncols), 1             );
     uu   ->Resize( nrows, nrows         );
     vv   ->Resize( ncols, ncols         );
     dd   ->Resize( MIN(nrows,ncols), MIN(nrows,ncols)  );
     rr   ->Resize( ncols, nrows         );

     dd->FillWithZeros();

     *ww = SVD( A, uu , vv  );


     tol = MAX( nrows, ncols )* Max( *ww )* eps;


     for( int i=1; i<= ww->GetNoRows(); i++ ) {

       if( (*ww)(i) > tol ) (*dd)(i,i) = 1.0/(*ww)(i);

       else  (*dd)(i,i) = 0.0;

     }


     *uu = uu->sub_matrix(1,nrows, 1, MIN(nrows,ncols) );
     *vv = vv->sub_matrix(1,ncols, 1, MIN(nrows,ncols) );
     (*rr) = (*vv)*ProductT( *dd , *uu );

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return (*rr);

}
#endif /* LAPACK */




#ifdef LAPACK
DMatrix& orth( const DMatrix& A  )
{
     // Q = orth( A ) is the orthonormal basis for the range of A
     // Q*Q' == I, the number of columns of Q is the rank of A

     int localFlag = 0;

     double tol;

     double eps = DMatrix::GetEPS();

     int r;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );


     int nrows = A.GetNoRows();
     int ncols = A.GetNoCols();

     DMatrix * ww;
     DMatrix * uu;
     DMatrix * vv;

     ww     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     uu     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     vv     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     ww   ->Resize( MIN( nrows, ncols ) , 1             );
     uu   ->Resize( nrows, ncols         );
     vv   ->Resize( ncols, ncols         );

     *ww = SVD( A, uu , vv  );


     tol = MAX( nrows, ncols )* Max( *ww )* eps;

     r = (int) ( sum( *ww > tol ) )(1,1);

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return  uu->sub_matrix( 1, uu->GetNoRows(), 1, r );


}
#endif /* LAPACK */


#ifdef LAPACK
DMatrix& null( const DMatrix& A  )
{
     // Z = null( A ) is the orthonormal basis of the null space of A
     // Z*Z' == I, A*Z == 0, the number of columns of Z is the
     // nullity of  A

     int localFlag = 0;

     int count = 0;

     double tol;

     double eps = DMatrix::GetEPS();

     int r, i;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );


     int nrows = A.GetNoRows();
     int ncols = A.GetNoCols();

     DMatrix * ww;
     DMatrix * uu;
     DMatrix * vv;
     DMatrix * Temp;

     ww     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     uu     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     vv     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     Temp   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     ww   ->Resize( MIN(ncols,nrows) , 1             );
     uu   ->Resize( nrows, ncols         );
     vv   ->Resize( ncols, ncols         );

     *ww = SVD( A, uu , vv  );

    tol = MAX( nrows, ncols )* Max( *ww )* eps;

    r = (int) ( sum( *ww > tol ) )(1,1);


    if ( r == ncols ) {

       if(!localFlag) DMatrix::SetMemberFlag( 0 );

       vv->Resize( 0, 0 );

       return *vv;

    }

    else {

      for( i= 1; i<= ncols; i++ ) {

         Temp->Resize( ncols, r );

         if ( (*ww)(i) < tol ) {

           count ++;
           Temp->SetColumn( vv->Column(i), count );

         }
      }

      if(!localFlag) DMatrix::SetMemberFlag( 0 );

      return  *Temp;

    }

}
#endif /* LAPACK */


#ifdef LAPACK
DMatrix& SVDSolve( const DMatrix& A, const DMatrix& B )
{
/*
 *   This function calls CLAPACK's routine DGELSS to compute
 *   the minimum norm solution to a real linear least
 *   squares problem:
 *
 *   Minimize 2-norm(| b - A*x |).
 *
 *   using the singular value decomposition (SVD) of A. A is an M-by-N
 *   matrix which may be rank-deficient.
 */

     int localFlag = 0;

     double eps = DMatrix::GetEPS();

     integer N = A.GetNoCols();

     integer M = A.GetNoRows();

     integer NRHS = B.GetNoCols();

     integer LDA  = A.GetNoRows();

     integer LDB  = B.GetNoRows();

     integer LWORK;

     integer RANK;

     integer INFO;

     double  RCOND;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     DMatrix * ww;
     DMatrix * aa;
     DMatrix * bb;
     DMatrix * wk;
     DMatrix * xx;

     ww     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     aa     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     bb     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     wk     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     xx     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));


     ww   ->Resize( MIN(M,N), 1             );
     aa   ->Resize( M, N         );
     bb   ->Resize( B.n, B.m         );
     xx   ->Resize( N, NRHS         );

     RCOND = MAX( N, M )* eps;
     LWORK = DMatrix::GetDimAux();


     memcpy( aa->a, A.a, A.n*A.m*sizeof(double)  );

     memcpy( bb->a, B.a, B.n*B.m*sizeof(double ) );


	 dgelss_( &M, &N, &NRHS, aa->a, &LDA, bb->a, &LDB, ww->a,

	          &RCOND, &RANK, wk->a, &LWORK, &INFO );

     *xx = bb->sub_matrix( 1, N, 1, NRHS );

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return  *xx;

}
#endif /*LAPACK*/


#ifdef LAPACK
DMatrix& schur( const DMatrix& A, DMatrix* U  )
{

      // This function computes the Schur decomposition of A
      // It calls dgees_(), which is part of CLAPACK, a
      // public domain numerical linear algebra library

     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     long nrows = A.GetNoRows();
     long ncols = A.GetNoCols();

     char JOBVS;
     char SORT  = 'N';
     integer N  = A.GetNoCols();
     integer LDA = A.GetNoRows();
     integer SDIM = 0;
     integer LDVS = N;
     integer LWORK = DMatrix::GetDimAux();
     integer INFO;

     DMatrix * a;
     DMatrix * wr;
     DMatrix * wi;
     DMatrix * wk;

     double *vs;

     if ( nrows != ncols ) error_message("Matrix must be square in schur()");


     a     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     wk    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     wr    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     wi    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     a   ->Resize( nrows, ncols   );
     wr  ->Resize( nrows,    1    );
     wi  ->Resize( nrows,    1    );

     if ( U == NULL )
     {

      JOBVS = 'N';
      vs = NULL;

     }

     else {

       JOBVS = 'V';

       U   ->Resize( nrows, ncols   );

       vs = U->GetPr();

     }

     memcpy( a->a, A.a, A.n*A.m*sizeof(double) );


      dgees_(&JOBVS, &SORT, NULL, &N,
        a->GetPr() , &LDA, &SDIM, wr->GetPr(),
        wi->GetPr(), vs   , &LDVS, wk->GetPr(),
        &LWORK, NULL, &INFO);

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     return *a;

}

#endif /* LAPACK */


#ifdef LAPACK

DMatrix& LSMNSolve( const DMatrix& A, const DMatrix& B  )
{


     int localFlag = 0;

     if (DMatrix::GetMemberFlag()) localFlag = 1;

     DMatrix::SetMemberFlag( 1 );

     char  TRANS = 'N';
     integer N =    A.GetNoCols();
     integer M =    A.GetNoRows();
     integer NRHS = B.GetNoCols();
     integer LDA =  M;
     integer LDB;
     integer LWORK = DMatrix::GetDimAux();
     integer INFO;

     integer MN = MAX( M, N );

     LDB = MAX( 1, MN );

     DMatrix * a;
     DMatrix * b;
     DMatrix * wk;


     a     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     b     = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     wk    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

     a   ->Resize( M, N     );
     b   ->Resize( LDB, NRHS  );

     b->FillWithZeros();

     memcpy( a->a, A.a, A.n*A.m*sizeof(double) );

     b->SetSubMatrix( 1,1, B );



     dgels_(&TRANS, &M, &N,
        &NRHS, a->GetPr(), &LDA, b->GetPr(), &LDB,
        wk->GetPr(), &LWORK, &INFO );

     if(!localFlag) DMatrix::SetMemberFlag( 0 );

     *b = b->sub_matrix(1,N,1,NRHS);

     return *b;
}

#endif /* LAPACK */


int isSymmetric( const DMatrix & A )
{

    int i, j;

    int isflag     = 1;

    if ( A.GetNoRows() != A.GetNoCols() ) {

       error_message( "Matrix must be square in isSymmetric()" );

    }

    for ( i = 1; i<= A.GetNoRows( ) && isflag!=0 ; i++ ) {

       for ( j = i; j<= A.GetNoRows( ) && isflag!=0; j++ ) {


          if ( A(i,j) == A(j,i) ) isflag = 1;
          else isflag = 0;

       }

    }

    return isflag;

}



#ifdef LAPACK

DMatrix& eig( const DMatrix & A, DMatrix* V  )
{

   int localFlag = 0;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix::SetMemberFlag( 1 );

   DMatrix* a;

   DMatrix* wr;

   DMatrix* wi;

   DMatrix* wk;

   DMatrix* vx;

   DMatrix* vz;

   int nn;

   int i,j;

   nn = A.GetNoRows();

   char JOBZ;
   char UPLO = 'U';
   integer LDA = nn;
   integer INFO;
   char JVL= 'N';
   char JVR;
   integer N   = nn;
   integer LDVL = nn;
   integer LDVR = nn;
   integer LWORK = DMatrix::GetDimAux();
   double   VL;
   double* vpr;





   a   = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
   wr  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
   wi  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
   wk  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
   vz  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   a    ->Resize( nn , nn            );
   wr   ->Resize( nn , 1             );
   wi   ->Resize( nn , 1             );
   wk   ->Resize( 5*nn, 1            );
   vz   ->Resize( nn ,  nn           );

   if (V!=NULL) {
     V    ->Resize( nn , 2*nn          );
     vpr  = vz->GetPr();
   }

   else {
     vx  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
     vpr = vx->GetPr();
   }


   memcpy( a->GetPr(), A.GetConstPr(), nn*nn*sizeof(double) );

   wi->FillWithZeros();

   if ( isSymmetric( A ) ) {


       if ( V==NULL )  JOBZ = 'N';

       else JOBZ = 'V';

/*       dsyev_(char *jobz, char *uplo, integer *n, doublereal *a,
         integer *lda, doublereal *w, doublereal *work, integer *lwork,
        integer *info)
*/       dsyev_(&JOBZ, &UPLO, &N, a->GetPr(),
         &LDA, wr->GetPr(), wk->GetPr(), &LWORK,
         &INFO);

//         memcpy( vpr , a->GetPr(), nn*nn*sizeof(double) );
           memcpy( V->GetPr() , a->GetPr(), nn*nn*sizeof(double) );

   }

   else {

   if ( V==NULL  ) { JVR = 'N'; }

   else { JVR = 'V'; }

//    DGEEV computes for an N-by-N real nonsymmetric matrix A, the
//    eigenvalues and, optionally, the left and/or right eigenvectors.
//    DGEEV is part of CLAPACK a public domain, numerical linear
//    algebra package

   dgeev_( &JVL, &JVR, &N , a->GetPr(),
        &LDA, wr->GetPr(), wi->GetPr(), &VL,
        &LDVL , vpr , &LDVR, wk->GetPr(),
        &LWORK, &INFO);

      if ( V!= NULL )	{
	for ( j=1; j<= nn; j++) {

	   if( (*wi)(j) == 0.0 ) {

	      V->SetColumn( vz->Column(j), j );
	      V->SetColumn( zeros(nn,1), j+nn );

	   }

	   else {

	      V->SetColumn( vz->Column(j) , j );
	      V->SetColumn(vz->Column(j+1), j+nn );
	      V->SetColumn( vz->Column(j), j+1 );
	      V->SetColumn( (-1.0)*(vz->Column(j+1)), j+nn+1 );
	      j++;

	   } /* End else */

        } /* End for */

      } /* End if */

   }


   if(!localFlag) DMatrix::SetMemberFlag( 0 );

   return (*wr || *wi);

}

#endif /* LAPACK */


#ifdef LAPACK

double rcond( const DMatrix& A )
{

   int localFlag = 0;

   if (DMatrix::GetMemberFlag()) localFlag = 1;

   DMatrix*   atmp;
   DMatrix*   wk;
   integer    N       = A.GetNoCols();
   integer    M       = A.GetNoRows();
   integer    LDA     = A.GetNoRows();
   integer    INFO;
//   integer   *IPIV    = new integer [ MIN( N, M ) ];
   integer    *IPIV   = (integer*) my_calloc( MIN(N,M), sizeof(integer) );
   char       NORM    = '1';
   double     ANORM;
   double     RCOND;

   atmp  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
   wk    = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   if ( N!= M ) error_message( "\n Matrix must be square in rcond()");

   ANORM = InfNorm( tra( A ) );

   atmp ->Resize( N, N );
   memcpy( atmp->a, A.a, N*N*sizeof(double) );

   dgetrf_(&M, &N, atmp->a, &LDA, IPIV, &INFO );

   dgecon_(&NORM,&N,atmp->a,&LDA,&ANORM,&RCOND,wk->a,IPIV,&INFO );

   if(!localFlag) DMatrix::SetMemberFlag( 0 );

   delete[] IPIV;

   return (RCOND);

}
#endif /* LAPACK */


double MachEps(void)
{
   // Function to compute machine precision
   double meps = 1.,eps = 2.;
   while(eps != 1.)
   {
      meps /= 2.;
      eps = 1. + meps;
   }

   return 2.0*meps;
}



double ArrayNorm(long n,double *x)
{
// Euclidean norm implementation with over or underflow protection
//
   const double maxDouble =                      1.7e+308;
   const double minDouble =                       8.e-307;
   double aux, xmax = 0.,xmin = maxDouble;

   for(long j = 0;j < n;j++)
   {
      aux = fabs(x[j]);
      if(xmax < aux)xmax = aux;
      if(xmin > aux)xmin = aux;
   }

   if(xmax == 0.)return xmax;

   if (xmin == 0.)xmin = minDouble;

   long double longaux =
      (long double)xmax/(long double)xmin;

   aux = sqrt(maxDouble/((double)n));

      if(xmax < aux &&
      xmax > minDouble/DMatrix::GetEPS() &&
      longaux < 1./DMatrix::GetEPS()  )

      {
      double norm = 0.;
      for(long i = 0;i < n;i++)
      {
        aux = x[i];
        norm += aux*aux;
      }
      return sqrt(norm);
   }
   else
   {
      long double norm = 0.;

      for(long i = 0;i < n;i++)
      {
        longaux = x[i];
        norm += longaux*longaux;
      }

      if(norm < maxDouble && norm > minDouble)
         return sqrt(norm);

      longaux = (long double)xmax*(long double)n;

      norm /= longaux; // avoids overflow

      norm /= longaux;

      norm = longaux*sqrt(norm); // renormalises

           // avoids overflow

      if(norm > maxDouble) norm = maxDouble;

      return norm;

   }
}


// Indexing functions

DMatrix& DMatrix::operator() (const DMatrix& RowIndx, const DMatrix& ColIndx )
{

      long nrows;
      long ncols;

      long ii;
      long jj;

      int rflag = 0;
      int cflag = 0;

      DMatrix* RtVal;

      if ( RowIndx(1) == -1 ) rflag = 1;
      if ( ColIndx(1) == -1 ) cflag = 1;


      if( !rflag ) nrows = RowIndx.GetNoRows();
      else nrows = GetNoRows();
      if( !cflag ) ncols = ColIndx.GetNoRows();
      else ncols = GetNoCols();

      RtVal  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
      RtVal->SetMType( 1 );
      RtVal->SetReferencedDMatrixPointer( this );
      RtVal->SetRowIndexPointer( &RowIndx );
      RtVal->SetColIndexPointer( &ColIndx );
      RtVal->Resize( nrows, ncols );


      for (long  i = 1; i<=nrows; i++ ) {

        if( !rflag ) ii = (long) RowIndx(i);
	else ii = i;

        for (long  j = 1; j<=ncols; j++ ) {

	   if (!cflag) jj = (long) ColIndx(j);
	   else jj=j;

	   RtVal->elem(i,j) = elem(ii,jj);

	}

      }

      return (*RtVal);
}

DMatrix& DMatrix::operator() (const DMatrix& RowIndx )
{

      long nrows;

      long ii;

      int flag = 0;

      DMatrix* RtVal;

      if (RowIndx(1)==-1.0) flag = 1;

      if (!flag) nrows = RowIndx.GetNoRows();
      else nrows = GetNoRows()*GetNoCols();


      RtVal  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
      RtVal->SetMType( 1 );
      RtVal->SetReferencedDMatrixPointer( this );
      RtVal->SetRowIndexPointer( &RowIndx );
      RtVal->SetColIndexPointer( NULL );
      RtVal->Resize( nrows, 1 );

      for ( long i = 1; i<=nrows; i++ ) {

           if( !flag ) ii= (long) RowIndx(i);
	   else ii = i;

	   RtVal->a[i-1] = a[ii-1];

      }

      return (*RtVal);

}

DMatrix& DMatrix::operator() (const DMatrix& RowIndx, long col )
{

      DMatrix* ColIndx;
      
      if (col > this->GetNoCols()) {
            this->Resize(this->GetNoRows(), col);
      }


      ColIndx  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
      ColIndx->SetMType( 0 );
      ColIndx->Resize( 1 , 1 );

      ColIndx->elem( 1,1 ) = (double) col;

      return (this->operator()(RowIndx,*ColIndx) );

}

DMatrix& DMatrix::operator() (const DMatrix& RowIndx, const char* end )
{
      long col = this->getm();

      DMatrix* ColIndx;

      ColIndx  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
      ColIndx->SetMType( 0 );
      ColIndx->Resize( 1 , 1 );

      ColIndx->elem( 1,1 ) = (double) col;

      return (this->operator()(RowIndx,*ColIndx) );

}


DMatrix& DMatrix::operator() (const char* end, const DMatrix& ColIndx )
{
      long row = this->getn();

      DMatrix* RowIndx;

      RowIndx  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
      RowIndx->SetMType( 0 );
      RowIndx->Resize( 1 , 1 );

      RowIndx->elem( 1,1 ) = (double) row;

      return (this->operator()(*RowIndx,ColIndx) );

}

DMatrix& DMatrix::operator() (long row, const DMatrix& ColIndx )
{
      DMatrix* RowIndx;
      
      if (row > this->GetNoRows()) {
            this->Resize(row, this->GetNoRows());
      }

      RowIndx  = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
      RowIndx->SetMType( 0 );
      RowIndx->Resize( 1 , 1 );

      RowIndx->elem( 1,1 ) = (double) row;

      return (this->operator()(*RowIndx,ColIndx) );

}


DMatrix& DMatrix::operator = ( double val )
{


    if ( GetMType()==1 ) {

       return ( AssignmentToColonReference( val ) );

    }

    else {

//		if ( n==0 || m == 0)
//		{
//			Resize(1,1);
//		}
//      long nelem = GetNoRows()*GetNoCols();
//      for (long i=0; i< nelem; i++ ) {
//         a[ i ] = val;
//      }
        Resize(1,1);
        a[0] = val;

      return (*this);

    }

}

DMatrix& DMatrix::AssignmentToColonReference(double val )
{
      long ii, jj;

      int rflag = 0;
      int cflag = 0;

      if ( (*rowIndx)(1) == -1.0 ) rflag = 1;
      if ( (*colIndx)(1) == -1.0 ) cflag = 1;


      for (long  i = 1; i<=GetNoRows(); i++ ) {

        if ( !rflag ) ii = (long) (*rowIndx)(i);

	else  ii = i;

        for (long  j = 1; j<=GetNoCols(); j++ ) {

	   if (!cflag)  jj = (long) (*colIndx)(j);
	   else jj= j;

	   mt->elem(ii,jj) = val;

	}

      }

      if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

      return (*this);

}


DMatrix& DMatrix::AssignmentToColonReference( const DMatrix& RightMt )
{

      long ii,jj;

      int rflag = 0;
      int cflag = 0;

      if ( (*rowIndx)(1) == -1.0 ) rflag = 1;
      if (colIndx!=NULL) {
        if ( (*colIndx)(1) == -1.0 ) cflag = 1;
      }

      if (GetNoRows() != RightMt.GetNoRows() ||
          GetNoCols() != RightMt.GetNoCols()    )
      {
	     ERROR_MESSAGE("Dimensions error in AssignmentToColonReference()");
      }

	  if ( Max(*rowIndx)>mt->GetNoRows()   )
	  {
		mt->Resize( (long) Max(*rowIndx) , mt->GetNoCols() );
	  }

          if (colIndx!=NULL) {
	    if ( Max(*colIndx)>mt->GetNoCols()   )
	    {
		mt->Resize( mt->GetNoRows() , (long) Max(*colIndx) );
	    }
          }


      for (long  i = 1; i<=GetNoRows(); i++ ) {

        if (!rflag) ii = (long) (*rowIndx)(i);
	else ii = i;

        for (long  j = 1; j<=GetNoCols(); j++ ) {
	   if (colIndx!=NULL) {
	     if (!cflag) jj = (long) (*colIndx)(j);
	     else jj = j;
           } else {jj=1;}
	   mt->elem(ii,jj) = RightMt.elem(i,j);

	}


      }

      if( !DMatrix::GetMemberFlag() ) DMatrix::SetAuxIndx( -1 );

      return (*this);

}


DMatrix& colon( int i1, int increment, int i2)
{
      DMatrix* Temp;

      long val;

      int nrows = abs( ( (i2-i1)/increment+1 ) );

      Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

      Temp->Resize( nrows, 1 );

      val = i1;

      for( long i=1; i<= nrows; i++) {

        Temp->elem(i,1) = val;

	val += increment;

      }

      return (*Temp);

}

DMatrix& colon( int i1, int i2 )
{
      DMatrix* Temp;

      long val;

      long increment = 1;

      long nrows = abs( ( (i2-i1)/increment+1 ) );

      Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

      Temp->Resize( nrows, 1 );

      val = i1;

      for( long i=1; i<= nrows; i++) {

        Temp->elem(i,1) = val;

	val += increment;

      }

      return (*Temp);

}

DMatrix& colon( double i1, double i2)
{

      DMatrix* Temp;

      double val;

      double increment = 1.0;

      long nrows = (int) fabs( ( (i2-i1)/increment )+1.0 );

      Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

      Temp->Resize( nrows, 1 );

      val = i1;

      for( long i=1; i<= nrows; i++) {

        Temp->elem(i,1) = val;

	val += increment;

      }


      return (*Temp);

}

DMatrix& colon( double i1, double increment, double i2)
{

      DMatrix* Temp;

      double val;

      long nrows =  (int) fabs( ( (i2-i1)/increment+1.0 ) );

      Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

      Temp->Resize( nrows, 1 );

      val = i1;

      for( long i=1; i<= nrows; i++) {

        Temp->elem(i,1) = val;

	val += increment;

      }

      return (*Temp);
}


DMatrix& colon( void )
{
      DMatrix* Temp;

      Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

      Temp->Resize( 1, 1 );

      Temp->elem(1,1) = -1.0;

      return (*Temp);

}



DMatrix& Sqrt( const DMatrix& A) {

   // Elementwise square root of a matrix

   DMatrix* Temp;

   Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize( A.n, A.m );

   for (long i= 0; i< A.n*A.m; i++ ) {

        if (A.a[i]<0) { error_message("sqrt(): Attempt to find square root of negative number\n...Complex numbers are not supported in this version of DMatrix");
        }
        Temp->a[i] = sqrt(A.a[i]);

   }

   return (*Temp);

}

DMatrix& triu( const DMatrix& A) {

   // Upper triangular part of a matrix

   DMatrix* Temp;

   Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

   Temp->Resize( A.n, A.m );

   Temp->FillWithZeros();

   for (long i= 1; i<= A.n; i++ ) {
      for (long j= i; j<= A.m; j++ ) {
        Temp->elem(i,j) = A.elem(i,j);
      }
   }

   return (*Temp);


}


double factorial(long j)
{

	double retval = 1;

        int k;

        if( j>1 ) {

	    for(k=2;k<=j;k++)
	       retval=retval*k;

        }

        return (retval);

}


DMatrix & eye(long n)
{
    return identity(n);
}

DMatrix & eye(long n, long m)
{
    return identity(n,m);
}

 long length(const DMatrix& A)
 {
	int retval;
	retval = A.GetNoRows()*A.GetNoCols();
	return retval;
 }

DMatrix&  find(const DMatrix& A)
{
// Returns a vector with the linear indices of the non-zero elements
// of the input matrix A.
    int i;
    int icount = 0;
    DMatrix* Temp;
    Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
    Temp->Resize( length(A), 1 );
    Temp->FillWithZeros();

    for (i=1; i<=length(A); i++)
    {
		if (A(i) != 0.0)
		{
			icount++;
			Temp->elem(icount,1) = i;
		}
    }
    Temp->Resize(icount,1);
    return (*Temp);
}

DMatrix& DMatrix::find(DMatrix& I, DMatrix& J) const
{
// Finds vector indices I and J of the non-zero elements of the calling DMatrix object, and
// returns a DMatrix object as a column vector with the non-zero elements of the calling object.
    int i, j;
    int icount = 0;
    DMatrix* Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
    Temp->Resize(n*m,1);
    Temp->FillWithZeros();
    I.Resize( n, 1 );
    J.Resize( m, 1 );
    I.FillWithZeros();
    J.FillWithZeros();

	for (i=1; i<=n; i++)
	{
		for (j=1; j<=m; j++)
		{
			if (this->elem(i,j) != 0.0)
			{
				icount++;
				I(icount) = i;
                                J(icount) = j;
                                (*Temp)(icount)= this->elem(i,j);
			}
		}
	}
	I.Resize(icount,1);
        J.Resize(icount,1);
        Temp->Resize(icount,1);
	return (*Temp);
}

DMatrix& DMatrix::find(int* I, int* J) const
{
// Finds array indices I and J of the non-zero elements of the calling DMatrix object, and
// returns a DMatrix object as a column vector with the non-zero elements of the calling object.
// Arrays I and J should be preallocated with size A.n*A.m


    int i, j;
    int icount = 0;
    DMatrix* Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
    Temp->Resize(n*m,1);
    Temp->FillWithZeros();

	for (i=1; i<=n; i++)
	{
		for (j=1; j<=m; j++)
		{
			if (this->elem(i,j) != 0.0)
			{
				I[icount] = i;
                                J[icount] = j;
                                icount++;
                                (*Temp)(icount)= this->elem(i,j);
			}
		}
	}
        Temp->Resize(icount,1);
	return (*Temp);
}


double DMatrix::random_uniform(void)

// returns a pseudo-random real number uniformly distributed
// between 0.0 and 1.0.

{
  const long Q = RAND_MODULUS / RAND_MULTIPLIER;
  const long R = RAND_MODULUS % RAND_MULTIPLIER;
        long t;

  t = RAND_MULTIPLIER * (seed[stream] % Q) - R * (seed[stream] / Q);
  if (t > 0)
    seed[stream] = t;
  else
    seed[stream] = t + RAND_MODULUS;
  return ((double) seed[stream] / RAND_MODULUS);
}


double DMatrix::random_gaussian(void)
// Returns a Gaussian random number with zero mean and unit variance
// The method uses the approximation of the normal idf due to Odeh and Evans,
// J. Applied Statistics, 1974, vol 23, pp 96-97.

{
  double p0 = 0.322232431088;
  double q0 = 0.099348462606;
  double p1 = 1.0;
  double q1 = 0.588581570495;
  double p2 = 0.342242088547;
  double q2 = 0.531103462366;
  double p3 = 0.204231210245e-1;
  double q3 = 0.103537752850;
  double p4 = 0.453642210148e-4;
  double q4 = 0.385607006340e-2;
  double u1, t, p, q, out;

  u1   = DMatrix::random_uniform();
  if (u1 < 0.5)
    t = sqrt(-2.0 * log(u1));
  else
    t = sqrt(-2.0 * log(1.0 - u1));
  p   = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
  q   = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
  if (u1 < 0.5)
    out = (p / q) - t;
  else
    out = t - (p / q);
  return out;
}



DMatrix& randu(long n, long m)
// Returns an nxm matrix where each element is a uniform random
// number in the range (0,1).
//
{

  long i;

  double *pr;

  DMatrix* Temp;

  double dummy = DMatrix::random_uniform();

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  pr   = Temp->GetPr();

  for (i=0; i< n*m; i++) {

      pr[i] = DMatrix::random_uniform();

  }

  return *Temp;
}

DMatrix& randn(long n, long m)
// Returns an nxm matrix where each element is a Gaussian random
// number with zero mean and unit variance.
//
{

  long i;

  double *pr;

  double dummy = DMatrix::random_gaussian();

  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( n, m );

  pr   = Temp->GetPr();

  for (i=0; i< n*m; i++) {

      pr[i] = DMatrix::random_gaussian();

  }

  return *Temp;
}

void tic(void)
{
// Chronometer starts
   DMatrix::SetStartTicks(clock());

}

double toc()
{
// Chronometer stops

   clock_t  stop_ticks     = clock();
   clock_t  elapsed_ticks;
   double elapsed_time;

   elapsed_ticks= stop_ticks-DMatrix::GetStartTicks();

   elapsed_time = ((double) elapsed_ticks/CLOCKS_PER_SEC);

   if (1) {

       if (DMatrix::PrintLevel() ) {
       	fprintf(stderr,"\n Elapsed time is %e seconds\n", elapsed_time );
       }

   }

   return (elapsed_time);

}

DMatrix& linspace(double X1, double X2, long N)
{
  long i;
  DMatrix* Temp;

  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));

  Temp->Resize( 1, N );

  if (N<1) {
    error_message("Length error in linspace()");
  }

  for (i=1; i<= N; i++) {

      (*Temp)(i) = X1 + ((double) (i-1))*(X2-X1)/( (double) (N-1) );

  }

  return *Temp;

}


void* my_calloc(size_t num, size_t size )
{
    void* pointer;
    pointer = mxCalloc(num, size);
    if (pointer==NULL) {
	error_message("Memory allocation error in my_calloc()");
    }
    return pointer;
}


// =====================================SPARSE MATRIX CLASS IMPLEMENTATION ==============================


SparseMatrix::SparseMatrix(void)
{
// Default constructor
   asize = 0;
   a     = 0;
   n     = 0;
   m     = 0;
   nz    = 0;

   a = NULL;
   RowIndx = NULL;
   ColIndx = NULL;
}


SparseMatrix::SparseMatrix( int nn, int mm, int nzz)
{
    n = nn;
    m = mm;
    nz= nzz;
    asize = nz;
    a       = new double[nz];
    RowIndx = new int[nz];
    ColIndx = new int[nz];

    for(int k=0;k<nz;k++)  {
       RowIndx[k]=0;
       ColIndx[k]=0;
       a[k]      = 0.0;
    }

}

SparseMatrix::SparseMatrix(double* aArg, int nArg, int mArg, int nzArg, int* RowIndxArg, int* ColIndxArg)
// Constructor using arrays
{
   int i;
   n     = nArg;
   m     = mArg;
   nz    = nzArg;
   asize = nz;

   a       = new double[nz];
   RowIndx = new int[nz];
   ColIndx = new int[nz];


   for(i=0;i<nz;i++)
   {
         a[i]       = aArg[i];
         ColIndx[i] = ColIndxArg[i];
         RowIndx[i] = RowIndxArg[i];
   }

//   memcpy( a, aArg, nz*sizeof(double) );
//   memcpy( ColIndx, ColIndxArg, nz*sizeof(int));
//   memcpy( RowIndx, RowIndxArg, nz*sizeof(int));
}



void SparseMatrix::Resize(int nnew, int mnew, int nznew)
{

    double *anew;
    int    *RowIndxNew;
    int    *ColIndxNew;
    int    dmax;

    anew = new double[nznew];
    RowIndxNew = new int[nznew];
    ColIndxNew = new int[nznew];

    dmax = MIN( nznew, nz );

    if (a!= NULL)
    {
        memcpy(anew, a, dmax*sizeof(double) );
        delete a;
    }

    if (RowIndx!=NULL) {
    	memcpy(RowIndxNew, RowIndx, dmax*sizeof(int) );
	delete RowIndx;
    }

    if (ColIndx!=NULL) {
    	memcpy(ColIndxNew, ColIndx, dmax*sizeof(int) );
        delete ColIndx;
    }

    a       = anew;
    RowIndx = RowIndxNew;
    ColIndx = ColIndxNew;

    n  = nnew;
    m  = mnew;
    nz = nznew;

}


SparseMatrix::SparseMatrix( const DMatrix& A ) // Constructor using DMatrix object
{
     DMatrix Anz;
     int* I = new int[A.GetNoRows()*A.GetNoCols()];
     int* J = new int[A.GetNoRows()*A.GetNoCols()];

     Anz = A.find( I, J);

     n  = A.GetNoRows();
     m  = A.GetNoCols();
     nz = length(Anz);

     a = new  double[nz];
     RowIndx = new int [nz];
     ColIndx = new int [nz];

     memcpy(a, Anz.GetPr(), nz*sizeof(double) );
     memcpy(RowIndx, I, nz*sizeof(int) );
     memcpy(ColIndx, J, nz*sizeof(int) );

    delete[] I;
    delete[] J;

}


SparseMatrix::SparseMatrix( const SparseMatrix& A) // copy constructor
{
   asize = A.nz;
   a     = A.a;
   n     = A.n;
   m     = A.m;
   nz    = A.nz;

   a = new double[nz];
   RowIndx = new int[nz];
   ColIndx = new int[nz];

   memcpy(a, A.a, nz*sizeof(double) );
   memcpy(ColIndx, A.ColIndx, nz*sizeof(int) );
   memcpy(RowIndx, A.RowIndx, nz*sizeof(int) );
}

// Destructor
SparseMatrix::~SparseMatrix()
{
    if (a!=NULL)
 	delete a;
    if (RowIndx!= NULL)
        delete RowIndx;
    if (ColIndx!= NULL)
        delete ColIndx;

}


DMatrix& full(const SparseMatrix& A)
{
    DMatrix* R;

    int k;

    int *rowA = A.GetRowIndx_C_Array();

    int *colA = A.GetColIndx_C_Array();

    double* Aij = A.GetPr();

    R = new DMatrix( A.GetNoRows(), A.GetNoCols() );

    R->FillWithZeros();

    for (k=0;  k< A.GetNonZero(); k++ )
    {
          R->elem( rowA[k], colA[k] ) = Aij[k];
    }

    return (*R);

}
