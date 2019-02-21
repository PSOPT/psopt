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

/*! \mainpage DMatrix and SparseMatrix classes
    \section intro Introduction
    The author developed the main features of the DMatrix class between 1994 and 1999. The class makes extensive use of operator overloading in order to facilitate the implementation in C++ of complicated matrix expressions, and it has intefaces to a number of LAPACK routines. In 2008, the class has been tested with current compilers, its functionality was expanded, and the code was published under the GNU Lesser General Public License. The class at present is restricted to dense and real matrices. In 2008, the SparseMatrix class was added to the library to incorporate basic sparse matrix functionality. The SparseMatrix class offers interfaces to some functions available in the CXSparse and LUSOL libraries.


The library should compile without problems with the following C++ compilers: GNU C++ version 4.X and Microsoft Visual Studio 2005.

    \section license License

This work is copyright (c) Victor M. Becerra (2009)

This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA, or visit http://www.gnu.org/licenses/

Author:    Dr. Victor M. Becerra,  University of Reading, School of Systems Engineering, P.O. Box 225, Reading RG6 6AY, United Kingdom, e-mail: v.m.becerra@ieee.org.

  \section install Installing the library
          See the INSTALL file.
  \section examples Examples of use
          See the source code in the examples directory.
*/


#ifndef  DMATRIX_H
#define  DMATRIX_H

#ifndef bool
// typedef int bool;
#endif

#ifndef N_TEMP_OBJECTS
#define N_TEMP_OBJECTS  (40)
#endif

#ifndef D_TEMP_OBJECTS
#define D_TEMP_OBJECTS (1000000)
#endif

#ifndef OUTPUT_STREAM
#define OUTPUT_STREAM stderr
#endif


#ifdef SPARSE_MATRIX

#include "cs.h"

#endif

#ifdef UNIX


extern "C" {


#include <stdio.h>
#include <stdarg.h>


}


#else


#include <stdio.h>
#include <stdarg.h>

#ifdef WIN32
#include <stdlib.h>
#include <malloc.h>
#endif

#ifndef WIN32
#include <mem.h> static long start_clock;
#endif

#endif

#ifdef WIN32
#define DEC_THREAD __declspec( thread )
#else
#define DEC_THREAD __thread
#endif


#ifndef MATLAB_MEX_FILE

#define mxCalloc calloc
#define mxFree   free


#endif /* MATLAB_MEX_FILE */


#include <time.h>

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

#include <string>
using std::string;


#define RAND_MODULUS    2147483647 /* DO NOT CHANGE THIS VALUE                  */
#define RAND_MULTIPLIER 48271      /* DO NOT CHANGE THIS VALUE                  */
#define RAND_CHECK      399268537  /* DO NOT CHANGE THIS VALUE                  */
#define RAND_STREAMS    256        /* # of streams, DO NOT CHANGE THIS VALUE    */
#define RAND_A256       22925      /* jump multiplier, DO NOT CHANGE THIS VALUE */
#define RAND_DEFAULT    123456789  /* initial seed, use 0 < RAND_DEFAULT < RAND_MODULUS  */


class DMatrix;
class SparseMatrix;

//! DMatrix class

/**
    A C++ class for dense and real matrix and vector computations with interfaces to
    a number of LAPACK functions
*/


class DMatrix {         // define DMatrix class

protected:
//! Array of doubles to store matrix elements using column major storage
   double *a;
  //! Number of matrix rows
   long n;
  //! Number of matrix columns
   long m;
//! Number of allocated elements in a
   long asize;
/*! Flag to indicate type of allocation.
    type = 0 : allocated matrix
    type = 1 : non-allocated matrix, uses predefined array for storage
*/
   int atype;
/*! Flag to indicate type of matrix
    mtype = 0 : normal matrix
    mtype = 1 : colon - reference matrix
*/
   int mtype;
 //! Flag to indicate auxiliary (temporary) matrix flag = 1 if temporary matrix, 0 otherwise.
    int auxFlag;
 //! Flag to indicate that element storage has been allocated
    int allocated;
//! Referenced matrix pointer
    DMatrix* mt;
//! Row indices
    const DMatrix* rowIndx;
//! Column indices
    const DMatrix* colIndx;

   // Protected static members
//! Array of Temporary matrices
   static DEC_THREAD DMatrix* auxPr;
//! Number of auxiliary matrices
   static DEC_THREAD int      noAuxArr;
//! Dimension of auxiliary arrays
   static DEC_THREAD long      dimAux;
//! Index of used auxiliary matrices
   static DEC_THREAD int      auxIndx;
//! Member function flag to control resetting of auxIndx
   static DEC_THREAD int      memberFlag;
//! Flag to indicate aux arrays allocation
   static DEC_THREAD int      initFlag;
//! Machine precision constant
   static DEC_THREAD double   MACH_EPS;
//! Flag to indicate error condition
   static DEC_THREAD int      errorFlag;
//! Print level flag,  1: output sent to sderr, 0: no output sent
   static DEC_THREAD int      print_level;

//! current state of each stream
   static DEC_THREAD long *seed;
//! stream index for pseudo-ramdon number generator
   static DEC_THREAD int  stream;
//! variable to store start time after tic() call.
   static DEC_THREAD time_t  start_time;
 //! clock_t variable
   static DEC_THREAD clock_t start_clock;
//! returns pointer to the i-th temporary object
   static DMatrix*  GetTempPr(int i);
//! Gets the memberFlag value from the object
   static int    GetMemberFlag() { return memberFlag; }
//! Sets the memberFlag value
   static void   SetMemberFlag(int arg) { memberFlag=arg; }
//! Checks the auxiliary objects
   static void   ChkAuxArrays();
//! Gets the number of temporary objects
   static int    GetNoAuxArr() { return noAuxArr; }
//! Gets the dimensions of the array of temporary objects
   static long    GetDimAux() { return dimAux; }
//! Gets the value of initFlag
   static int    GetInitFlag() { return initFlag; }
//! Sets the value of initFlag
   static void   SetInitFlag(int arg) { initFlag = arg ; }
//! Gets the current index of temporary objects
   static int    GetAuxIndx()       { return auxIndx; }
//! Increments the index of temporary objects
   static int    IncrementAuxIndx() {
#ifdef CHECK_AUX
		fprintf(stderr,"\nauxIndx=%d",auxIndx);
#endif
		return ++auxIndx; }
//! Decrements the index of temporary objects
   static int    DecrementAuxIndx() { return --auxIndx ; }
//! Sets the index of temporary objects
   static void   SetAuxIndx( int i ) {  auxIndx = i; }
//! Sets the dimension of each object in the array of temporary objects
   static void   SetDimAux( long dd ) { dimAux =  dd; }
//! Sets the number of auxiliary objects
   static void   SetNoAuxArr( int nn ) { noAuxArr = nn ; }
//! Sets the errorFlag member to true
   static void   RiseErrorFlag()   { errorFlag = true; }
//! Sets the value of auxFlag
   void   SetAuxFlag(int arg) { auxFlag = arg; }
//! Gets the value of auxFlag
   int   GetAuxFlag() { return auxFlag; }
 //! Gets the value of start_clock member
   static clock_t GetStartTicks(void) { return start_clock; }
//! Sets the value of start_clock member
   static void SetStartTicks(clock_t st) { start_clock=st; }

#ifdef DECLARED_TEMPS

   static DEC_THREAD double axp[N_TEMP_OBJECTS][D_TEMP_OBJECTS];

#endif


#define MC_EPSILON 2.221e-16

  //! Cholesky decomposition of a matrix
  /**
      \param A is a DMatrix object
      \param n is the number of columns of the input matrix A.
      \param pM (modified), is a DMatrix object which on output contains the Cholesky decomposition of input matrix a.
      \return void
      \sa CholeskySolution()
 */

   friend void CholeskyDecomp(DMatrix& A, int n, DMatrix& pM);

  //! Cholesky solution using the Cholesky decomposition of a matrix
  /**
      \param A  DMatrix object
      \param n is the number of columns of the input matrix A.
      \param pM DMatrix object which on output contains the Cholesky decomposition of input matrix a.
      \param bM DMatrix object with a vector of right-hand-side values
      \param xM (modified) DMatrix object which on return contains the solution to the system of equations
      \return void
      \sa CholeskyDecomp()
 */

   friend void CholeskySolution(const DMatrix& A, int n, const DMatrix& pM,
                         const DMatrix& bM, DMatrix& xM);

  //! Allocate memory to store matrix elements
  /**
      \param size  number of double elements to allocate
      \return void
      \sa DeAllocate()
 */

   void Allocate(long size);

  //! De-allocate memory previously allocated with DMatrix::Allocate().
  /**
      \return void
      \sa Allocate()
 */

   void DeAllocate(); // De-allocate memory


  //! Elementwise comparison of matrix elements
  /**
      \param m2  DMatrix object to be compared with the calling object.
      \param op  (char) indicates type of operator. Use 1 for >, 2 for >=,. 3 for <,  4 for <=, 5 for ==, 6 for != comparisons
      \return A reference to a DMatrix object such that its elements are 1 if the comparison is true, 0 otherwise.
  */
   DMatrix& compMat( const DMatrix& m2, char op ) const;

  //! Sets the value of the referenced matrix pointer
  /**
      \param arg  Pointer to DMatrix object
      \return void
  */
   void SetReferencedDMatrixPointer( DMatrix* arg ) { mt = arg; }
  //! Gets the value of the referenced matrix pointer
  /**
      \return pointer to DMatrix object
  */
   DMatrix* GetReferencedDMatrixPointer() { return mt; }
  //! Sets the row index pointer
  /**
      \param arg  Pointer to DMatrix object
      \return void
  */
   void SetRowIndexPointer(const  DMatrix* arg ) { rowIndx = arg; }
  //! Sets the column index pointer
  /**
      \param arg  Pointer to DMatrix object
      \return void
  */
   void SetColIndexPointer(const  DMatrix* arg ) { colIndx = arg; }
  //! Sets the type of matrix
  /**
      \param arg  (int) should be 0 for a normal matrix, or 1 for a colon reference matrix
      \return void
  */
   void SetMType( int arg );
  //! Gets the type of matrix
  /**
      \return int value
  */
   int GetMType()  { return mtype; }
  //! Assigns a matrix to the values pointed to by a colon reference matrix
  /**
      \param  A DMatrix object
      \return reference to DMatrix object
  */
   DMatrix& AssignmentToColonReference( const DMatrix& A );
  //! Assigns a double value to each value pointed to by a colon reference matrix
  /**
      \param  arg is a double value to be assigned.
      \return reference to DMatrix object
  */
   DMatrix& AssignmentToColonReference( double arg );

public:


#ifdef SPARSE_MATRIX
   friend class SparseMatrix;
#endif

   // Public methods
  //! Allocates the array of auxiliary (temporary) objects used by the class
  /** Allocates a DMatrix array of size N_TEMP_OBJECTS. Each element is allocated a storage of
      size D_TEMP_OBJECTS. These two macros are given default values but may be changed by
      the user at compilation time. The purpose of the array of temporary objects is to
      to store the intermediate objects resulting from single lines of code that
      call various operators and functions returning DMatrix objects. A simple example is as
      follows. Consider the C++ statement
        D= A*B + C;
      where A,B C and D are DMatrix objects. This statement involves two temporary objects:
      one to store the result of A*B, and another one to store the result of (A*B)+C;
      \return void
  */
   static void   AllocateAuxArr( void );
  //! De-allocates the array of auxiliary (temporary) objects previously allocated by AllocateAuxArr()
  /**
      \return void
  */
   static void   DeAllocateAuxArr( void );
  //! Returns a pseudo-random uniformly distributed number in the range [0,1]
  /**
      \return double pseudo-random value in the range [0,1]
  */
   static double random_uniform(void);
  //! Returns a pseudo-random Gaussian distributed number with zero mean and unit variance
  /**
      \return double Gaussian pseudo-random value
  */
   static double random_gaussian(void);
  //! Gets a pointer to the array of auxiliary objects
  /**
      \return DMatrix** pointer
  */
   static DMatrix**  GetAuxPr(void)   { return &auxPr; }
  //! Checks if the error flag has been raised. If so, a 1 is returned, 0 otherwise
  /**
      \return int value
  */
   static int isThereError(void)      { return errorFlag; }

   // Input / Output functions

  //! Allows the user to enter the elements of a matrix using command line prompts
  /**
      \return void
  */
   void input_matrix();

  //! Prints the elements of a DMatrix object
  /**
      \param  text is a string that serves as a prompt for the matrix being entered.
      \return void
  */
   void Print(const char *text) const;

  //! Prints information about a DMatrix object
  /**
      \param  text is a a string that serves as a label for the matrix being printed.
      \return void
  */
   void PrintInfo(const char *text) const;
  //! Reads the elements of a matrix from a file
  /**
      The file in question should contain the elements of the matrix row by row
      The elements should be separated by spaces and each row should be separated by a new line.
      The calling object should have the appropriate number of rows and columns.
      \param  filex pointer to a file already opened using "fopen()".
      \return void
  */
   void Read(FILE *filex);
  //! Saves the elements of a matrix to a file
  /**
      \param  FileName is a string with the desired file name
      \return void
  */
   void Save( const char * FileName );
  //! Reads the elements of a matrix from a file
  /**
      The calling object should have the appropriate number of rows and columns.
      \param  FileName is a string with the file name where the matrix elements are stored.
      \return void
  */
   void Load( const char * FileName );

  //! Prints the elements of a matrix elements matrix to a file
  /**
      \param  filex is a pointer to a file already opened using "fopen()".
      \return void
  */
   void Fprint( FILE *filex );

  //! Sets the print level
  /**
      \param  plevel desired print level
      \return void
  */
   static void SetPrintLevel( int plevel );

   static int PrintLevel();

   // Modifying member functions
  //! Assigns a zero value to each element of a matrix.
  /**
      \return void
  */
   void FillWithZeros( void );
  //! Swaps two rows of a matrix
  /**
      \param i: first row index
      \param j: second row index.
      \return void
  */
   void SwapRows( int i, int j );
  //! Swaps two columns of a matrix
  /**
      \param i: first column index
      \param j: second column index.
      \return void
  */
   void SwapColumns( int i, int j );
  //! Transposes a matrix
  /**
      \return void
  */
   void Transpose(void);
  //! Assigns values to the diagonal elements of a matrix, while all off-diagonal elements are set to zero.
  /**
      The dimensions of dd should be consistant with the dimensions of the calling object
      \param  dd: reference to constant DMatrix object which should contain a vector with the desired diagonal elements
      \return void
  */
   void diag( const DMatrix& dd );
  //! Assigns values to a column of a matrix, while other columns are left untouched.
  /**
      \param  Col: reference to constant DMatrix object which should contain a vector with the desired column values
      \param  icol: index to the column that is to be changed
      \return void
  */
   void SetColumn( const DMatrix& Col, int icol );
  //! Assigns values to a row of a matrix, while other columns are left untouched.
  /**
      \param  Row: reference to constant DMatrix object which should contain a vector with the desired row values
      \param  irow: index to the row that is to be changed
      \return void
  */
   void SetRow(    const DMatrix& Row, int irow );
  //! Multiples the elements of a specified column of a matrix by a constant scalar value
  /**
      \param  c:  index to the column that is to be changed
      \param  x:  scalar value
      \return void
  */
   void colMult(long c, double x  );
  //! Multiples the elements of a specified row of a matrix by a constant scalar value
  /**
      \param  r:  index to the row that is to be changed
      \param  x:  scalar value
      \return void
  */
   void rowMult(long r, double x  );
  //! Finds non-zero values of a matrix
  /**
      \param  I:  DMatrix object with the row index of each non-zero element
      \param  J:  DMatrix object with the column index of each non-zero element
      \return DMatrix object with the same dimensions as the calling object, and with elements which are 0 if the corresponding element of the calling object is 0, 1 otherwise.
  */
   DMatrix& find(DMatrix& I, DMatrix& J) const;
  //! Finds non-zero values of a matrix
  /**
      \param  I:  C++ double array with the row index of each non-zero element
      \param  J:  C++ double array with the column index of each non-zero element
      \return DMatrix object with the same dimensions as the calling object, and with elements which are 0 if the corresponding element of the calling object is 0, 1 otherwise.
  */
   DMatrix& find(int* I, int* J) const;


   // Sub-matrix functions
  //! Extracts a specified sub-matrix from a matrix
  /**
      \param  r1: start of row range
      \param  r2: end of row range
      \param  c1: start of column range
      \param  c2: end of column range
      \return DMatrix object with the specified sub-matrix
  */
   DMatrix& sub_matrix(long r1, long r2, long c1, long c2) const ;
  //! Assigns the elements of a matrix object to a section of the calling object.
  /**
      \param  row: start of row range
      \param  col: start of column range
      \param  A: DMatrix object whose element values are to be copied into the calling object
      \return void
  */
   void SetSubMatrix(long row, long col, const DMatrix& A);

   // Interface functions
  //! Gets the number of rows from the calling object
  /**
      \return long value with the number of rows.
  */
   long getn() {return n;}
  //! Gets the number of columns from the calling object
  /**
      \return long value with the number of columns
  */
   long getm() {return m;}
  //! Gets the number of rows from the calling object
  /**
      \return long value with the number of rows.
  */
   long getn() const {return n;}
  //! Gets the number of columns from the calling object
  /**
      \return long value with the number of columns
  */
   long getm() const {return m;}
  //! Gets the number of rows from the calling object
  /**
      \return long value with the number of rows.
  */
   long GetNoRows() { return n; }
  //! Gets the number of rows from the calling object
  /**
      \return long value with the number of rows.
  */
   long GetNoRows() const  { return n; }
  //! Gets the number of columns from the calling object
  /**
      \return long value with the number of columns
  */
   long GetNoCols() { return m; }
  //! Gets the number of columns from the calling object
  /**
      \return long value with the number of columns
  */
   long GetNoCols() const { return m; }
  //! Gets the pointer to the array where the elements of the matrix are stored
  /**
      \return double pointer
  */
   double *GetPr()    { return a; }
  //! Gets a pointer to the array where the elements of the matrix are stored (for const objects)
  /**
      \return double pointer
  */
   double *GetConstPr() const { return a; }
  //! Gets a pointer to the array where the elements of the matrix are stored
  /**
      \return double pointer
  */
   double* geta() {return a;}
  //! Gets the type of element storage of a matrix
  /**
      \return int value: 0 if the array is allocated, 1 if the storage is done using a previously declared array of doubles.
  */
   int getatype() { return atype; }
  //! Check if a matrix object is empty, if other words this method checks if the calling object has zero elements
  /**
      \return int value: 0 if calling object is not empty, 1 otherwise
  */
   int isEmpty() { if (n*m == 0) return 1; else return 0; }

  //! Check if a matrix object contains a row or column vector
  /**
      \return int value: 1 if calling object is a vector, 0 otherwise
  */
   int isVector() const { if (n==1 || m==1) return 1; else return 0; }

   // methods to access to matrix elements, columns and rows
  //! Returns the value of a specified element of a matrix
  /**
      \param i: row index (starting from 1)
      \param j: column index (starting from 1)
      \return double value
  */
   double element(long i, long j );
  //! Returns a reference to the specified element of a matrix
  /**
      \param i: row index (starting from 1)
      \param j: column index (starting from 1)
      \return double reference
  */
   double&   elem( long i, long j ) { return a[ (j-1)*n + i-1 ]; }
  //! Returns the value of a specified element of a matrix
  /**
      \param i: row index (starting from 1)
      \param j: column index (starting from 1)
      \return double value
  */
   double    elem( long i, long j ) const { return a[(j-1)*n+i-1]; }
  //! Returns a DMatrix object containing a speficied column of the calling object
  /**
      \param icol: column index (starting from 1)
      \return DMatrix object with the specified column
  */
   DMatrix& Column(long icol ) const;
  //! Returns a DMatrix object containing a speficied row of the calling object
  /**
      \param irow: row index (starting from 1)
      \return DMatrix object with the specified row
  */
   DMatrix& Row(long irow ) const;

  //! Computes and returns the integer power of a matrix.
  /**
      \param p: integer value, starting from 0
      \return DMatrix object with the result of the calculation
  */
   DMatrix& mpow(int p);

   void initVars();
   // Constructors

  //! Default constructor. Creates an empty matrix.
   DMatrix(void); // Default constructor
  //! Constructor with dimensions. Creates a matrix with allocated storage given specified numbers of rows and columns.
  /**
      \param Initn: number of rows
      \param Initm: number of columns
  */
   DMatrix( long Initn, long Initm);
  //! Constructor with dimensions using pre-allocated storage.
  /**
      \param vDim:  Allocated length of array v
      \param v:     double array to be used as storage by the DMatrix object
      \param Initn: number of rows
      \param Initm: number of columns
  */
   DMatrix( long vDim, double* v, long Initn, long Initm );
  //! Constructor with a single dimension, creates a column vector.
  /**
      \param Initn: number of rows
  */
   DMatrix( long Initn); // Column vector constructor using memory allocation
  //! Constructor using a variable list of element values
  /**
      \param Initn:     number of rows
      \param Initm:     number of columns
      \param a11:       first element of list doubles, values are entered column by column
  */
   DMatrix(long Initn,long Initm,double a11,...);
  //! Copy constructor. Creates a new DMatrix object with the same dimensions and element values as a given DMatrix object.
  /**
      \param A:     DMatrix object to be copied
  */
   DMatrix( const DMatrix& A); // copy constructor


  //! Destructor. Destroys a previously created DMatrix object and frees any allocated memory.
   ~DMatrix();
   // Resizing and assignment
  //! Changes the number of rows and columns of an existing matrix. Allocates new memory if necessary. If the calling object uses preallocated memory and the requested size would exceed that memory, an error is thrown.
  /**
      \param nnrow: new number of rows
      \param nncol: new number of columns
  */
   void Resize(long nnrow, long nncol );
  //! Assigns values to the elements of an existing DMatrix object using a variable list of values. Resizes the matrix if necessary.
  /**
      \param rows:    number of rows
      \param columns: number of columns
      \param a11: first element of the list of double arguments.
  */
   void assign(long rows,long columns,double a11,...); // Assignment using variable list
  //! Copies values to the elements of a DMatrix object from an existing array. The number of elements of the DMatrix object is assumed to be the number of values to be copied.
  /**
      \param aptr:    pointer to the start of the array of doubles to be copied.
  */
   void MemCpyArray(double * aptr );

   // Operators
  //! Matrix addition operator. The sizes of the matrices being added should be the same, otherwise an error is thrown.
  /**
      \param rval:  matrix located at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator+ (const DMatrix& rval) const;
  //! Matrix addition and substitution operator. The sizes of the matrices being added should be the same, otherwise an error is thrown.The left hand side object elements are replaced with the result of the operation.
  /**
      \param rval:  matrix located right hand side of the operator.
      \return Reference to the calling object
  */
   DMatrix& operator += (const DMatrix &rval);
  //! Adds a scalar real value to each element of the matrix. .
  /**
      \param x:  double value to be added
  */
   DMatrix& operator+ (double x) const;
  //! Matrix subtraction and substitution operator. The sizes of the matrices being subtracted should be the same, otherwise an error is thrown.The left hand side object elements are replaced with the result of the operation.
  /**
      \param rval:  matrix located right hand side of the operator.
      \return Reference to the calling object
  */
   DMatrix& operator -= (const DMatrix &rval);
  //! Matrix subtraction operator. The sizes of the matrices being subtracted should be the same, otherwise an error is thrown.
  /**
      \param rval:  matrix located at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator- (const DMatrix& rval) const;
  //! Subtracts a scalar real value from each element of the matrix. .
  /**
      \param x:  double value to be subtracted
  */
   DMatrix& operator- (double x) const;
  //! Matrix unary minus operator. Returns an object of the same dimensions as A but with changed element signs.
  /**
      \param A:  matrix located at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& operator- (const DMatrix& A);
  //! Matrix product operator. Returns the result of the matrix product of the calling object (left hand side of the operator) and the right hand side object. The inner dimensions of the objects being multiplied should be consistent, otherwise an error will be thrown.
  /**
      \param rval:  matrix located at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator* (const DMatrix& rval) const;
  //! Matrix product operator with substitution. Computes the matrix product of the calling object (left hand side of the operator) and the right hand side object. The inner dimensions of the objects being multiplied should be consistent, otherwise an error will be thrown. The calling object is modified to store the results of the operation.
  /**
      \param rval:  matrix located at the right hand side of the operator.
      \return Reference to the calling DMatrix object
  */
   DMatrix& operator *= (const DMatrix &rval);
  //! Computes the product of a matrix (left hand side of the operator) times a real scalar (right hand side value) and replaces the left hand side object with the result of the operation.
  /**
      \param Arg: double value that will multiply each element of the matrix.
      \return Reference the calling DMatrix object
  */
   DMatrix& operator* (double Arg) const;
  //! Computes the product of a matrix (left hand side of the operator) times a real scalar (right hand side value), and modifies the calling object to store the result of the operation.
  /**
      \param Arg: double value that will multiply each element of the matrix.
      \return Reference to the calling DMatrix object
  */
   DMatrix& operator *= (double Arg);
  //! Computes the division of a matrix (left hand side of the operator) by a real scalar (right hand side value).
  /**
      \param Arg: double value that will divide each element of the matrix.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator/ (double Arg) const;
  //! Computes the right division of a matrix (left hand side of the operator) by another matrix (right hand side value). This is conceptually equivalent to multiplying the left object by the inverse of the right hand side object but it is computed in a more efficient way. The dimensions of the matrices must be consistent, otherwise an error is returned. The right hand side object must be a square matrix.
  /**
      \param rval: DMatrix object at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator/ (const DMatrix& rval) const;
  //! Computes the left division of a matrix (left hand side of the operator) by another matrix (right hand side value). This is conceptually equivalent to multiplying the inverse of the left object by the right hand side object but it is computed in a more efficient way. The dimensions of the matrices must be consistent, otherwise an error is returned. The left hand side object must be a square matrix.
  /**
      \param rval: DMatrix object at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator% (const DMatrix& rval) const;
  //! Computes the division of a matrix (left hand side of the operator) by a real scalar (right hand side value) and modifies the left hand side object with the result of the operation.
  /**
      \param Arg: double value that will divide each element of the matrix.
      \return Reference to the calling object
  */
   DMatrix& operator /= (double Arg);
  //! Matrix assignment. The size of the left hand side object is modified if necessary, and the values of all real elements of the right hand side object and copied to the left hand side object.
  /**
      \param rval: DMatrix object at the right hand side of the operator
      \return Reference to the calling object
  */
   DMatrix& operator= (const DMatrix& rval);
  //! Matrix assignment to a scalar. The size of the left hand side object is modified to one row by one column if necessary, and the value of the right hand side argument is copied to the single element of the matrix. If the calling object is a "colon reference" matrix, then the right hand side value is copied to each element of the referenced array elements.
  /**
      \param val: double value at the right hand side of the operator
      \return Reference to the calling object
  */
   DMatrix& operator = (double val);
  //! Matrix assignment to a constant matrix defined as a character string using the bracket notation used in Matlab and Octave. The size of the left hand side object is modified if necessary. For example, the identity matrix of size two by two would be entered as "[1.0  0.0;0.0  1.0]"
  /**
      \param str: Character string containing the constant matrix defined using Matlab/Octave bracket notation.
      \return Reference to the calling object
  */
   DMatrix& operator = (const char* str);
  //! Concatenates two matrices side by side. The dimensions number of rows of the matrices involved must be the same, otherwise an error is thrown. The number of columns of the resulting matrix is the addition of the number of columns of both matrices involved.
  /**
      \param B: DMatrix object at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator ||(const DMatrix& B) const;
  //! Stacks the right hand side matrix below the left hand side matrix. The dimensions number of columns of the matrices involved must be the same, otherwise an error is thrown. The number of rows of the resulting matrix is the addition of the number of rows of both matrices involved.
  /**
      \param B: DMatrix object at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator &&(const DMatrix& B) const;
  //! Elementwise power operator. Returns a DMatrix object with the same dimensions of the calling object and each of its elements is computed as the corresponding element of the calling object to the power of the right hand side argument. Care must be taken when using this operator as the associations do not work in the same way as with the * operator. It is highly recommended to use parenthesis every time this operator is used. For example use it as follows:   (A^x).
  /**
      \param x: double argument at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator^(double x);
  //! Elementwise product operator. Returns a DMatrix object with the same dimensions of the calling objects and each of its elements is computed as the product of the corresponding elements of the calling object and the right hand side object. The dimensions of the calling objects must be the same, otherwise an error is thrown. Care must be taken when using this operator as the associations do not work in the same way as with the * operator. It is highly recommended to use parenthesis every time this operator is used. For example use it as follows:   (A&B).
  /**
      \param B: DMatrix object at the right hand side of the operator
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator &(const DMatrix& B) const;
  //! Elementwise division operator. Returns a DMatrix object with the same dimensions of the calling objects and each of its elements is computed as the of the corresponding element of the calling object by the corresponding element of the right hand side object. The dimensions of the calling objects must be the same, otherwise an error is thrown. Care must be taken when using this operator as the associations do not work in the same way as with the / operator. It is highly recommended to use parenthesis every time this operator is used. For example use it as follows:   (A|B).
  /**
      \param B: DMatrix object at the right hand side of the operator
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator |(const DMatrix& B) const;
  //! Matrix indexing. Returns a reference to the matrix element located at the position indicated by the row and column indices. Indices start from 1.
  /**
      \param row: Row index starting from 1.
      \param col: Column index starting from 1.
      \return Reference to the indexed matrix element.
  */
   double& operator() (long row,long col);
  //! Matrix indexing. Returns a reference to the matrix element located at the position indicated by the row index and the last column. Indices start from 1. An error in thrown in case of zero or negative indices. The matrix is resized if necessary.
  /**
      \param row: Row index starting from 1.
      \param end: Character string containing the word "end".
      \return Reference to the indexed matrix element.
  */
   double& operator() (long row,const char* end);
  //! Matrix indexing. Returns a reference to the matrix element located at the position indicated by the column index and the last row. Indices start from 1. An error in thrown in case of a range violation.
  /**
      \param end: Character string containing the word "end".
      \param col: Column index starting from 1.
      \return Reference to the indexed matrix element.
  */
   double& operator() (const char* end,long col);
  //! Matrix indexing. Returns the value of the matrix element located at the position indicated by the row and column indices. Indices start from 1. An error in thrown in case of a range violation.
  /**
      \param row: Row index starting from 1.
      \param col: Column index starting from 1.
      \return double value of the indexed element.
  */
   double operator() (long row,long col) const;
  //! Matrix indexing. Returns the value of the matrix element located at the position indicated by the row index and the last column. Indices start from 1. An error in thrown in case of a range violation.
  /**
      \param row: Row index starting from 1.
      \param end: Character string containing the word "end".
      \return double value of the indexed matrix element.
  */
   double operator() (long row,const char* end) const;
  //! Matrix indexing. Returns the value of the matrix element located at the position indicated by the column index and the last row. Indices start from 1. An error in thrown in case of a range violation.
  /**
      \param end: Character string containing the word "end".
      \param col: Column index starting from 1.
      \return double value of the indexed matrix element.
  */
   double operator() (const char* end,long col) const;
  //! Single index matrix indexing. Returns a reference to the matrix element located at the linear position indicated by the index, assuming column major storage. The indexs start from 1. The matrix is resized if necessary. An error is thrown in case of zero or negative indices.
  /**
      \param index: index value
      \return reference to indexed matrix element.
  */
   double& operator() (long index );
  //! Access to last linear element. Returns a reference to the last linear matrix element, assuming column major storage.
  /**
      \param end: Character string containing the word "end".
      \return reference to indexed matrix element.
  */
   double& operator() (const char* end);
  //! Single index matrix indexing. Returns the value of the matrix element located at the linear position indicated by the index, assuming column major storage. The index starts from 1. An error is thrown in case of range error.
  /**
      \param k: index value
      \return Value of indexed matrix element.
  */
   double operator() (long k ) const;
  //! Access to last linear element. Returns the value of the last linear matrix element, assuming column major storage.
  /**
      \param end: Character string containing the word "end".
      \return value of the last matrix element.
  */
   double operator() (const char* end) const;
  //! Checks if each element of the DMatrix object on the left hand side is greater than the right hand side value. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise.
  /**
      \param val: right hand side value
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator >  ( double val ) const ;
  //! Checks if each element of the DMatrix object on the left hand side is lower than the right hand side value. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise.
  /**
      \param val: right hand side value
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator <  ( double val ) const ;
  //! Checks if each element of the DMatrix object on the left hand side is greater or equal than the right hand side value. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise.
  /**
      \param val: right hand side value
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator >= ( double val ) const ;
  //! Checks if each element of the DMatrix object on the left hand side is lower or equal than the right hand side value. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise.
  /**
      \param val: right hand side value
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator <= ( double val ) const ;
  //! Checks if each element of the DMatrix object on the left hand side is equal to the right hand side value. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise.
  /**
      \param val: right hand side value
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator == ( double val ) const ;
  //! Checks if each element of the DMatrix object on the left hand side is different from the right hand side value. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise.
  /**
      \param val: right hand side value
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator != ( double val ) const ;
  //! Elementwise matrix comparison. Checks if each element of the matrix on the left hand side is greater than the corresponding element of the right hand side matrix. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise. The dimensions of the two matrices involved must be the same, otherwise an error is thrown.
  /**
      \param val: val: right hand side object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator >  ( const DMatrix& val ) const ;
  //! Elementwise matrix comparison. Checks if each element of the matrix on the left hand side is lower than the corresponding element of the right hand side matrix. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise. The dimensions of the two matrices involved must be the same, otherwise an error is thrown.
  /**
      \param val: val: right hand side object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator <  ( const DMatrix& val ) const ;
  //! Elementwise matrix comparison. Checks if each element of the matrix on the left hand side is greater or equal than the corresponding element of the right hand side matrix. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise. The dimensions of the two matrices involved must be the same, otherwise an error is thrown.
  /**
      \param val: val: right hand side object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator >= ( const DMatrix& val ) const ;
  //! Elementwise matrix comparison. Checks if each element of the matrix on the left hand side is lower or equal than the corresponding element of the right hand side matrix. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise. The dimensions of the two matrices involved must be the same, otherwise an error is thrown.
  /**
      \param val: val: right hand side object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator <= ( const DMatrix& val ) const ;
  //! Elementwise matrix comparison. Checks if each element of the matrix on the left hand side is equal to the corresponding element of the right hand side matrix. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise. The dimensions of the two matrices involved must be the same, otherwise an error is thrown.
  /**
      \param val: val: right hand side object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator == ( const DMatrix& val ) const ;
  //! Elementwise matrix comparison. Checks if each element of the matrix on the left hand side is different from the corresponding element of the right hand side matrix. The result is a DMatrix object where each element has the value of 1 if the corresponding comparison was true, 0 otherwise. The dimensions of the two matrices involved must be the same, otherwise an error is thrown.
  /**
      \param val: val: right hand side object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   DMatrix& operator != ( const DMatrix& val ) const ;

   // Indexing functions
  //! Submatrix extraction and referencing using arrays of indices.
  /**
      \param RowIndx is a DMatrix array that contains row index values usually generated using the colon() function.
      \param ColIndx is a DMatrix array that contains column index values usually generated using the colon() function.
      \return Reference to a DMatrix mtype 1 object that maps to the referenced elements of the calling object.
  */
   DMatrix& operator() (const DMatrix& RowIndx, const DMatrix& ColIndx );
  //! Linear sub-vector extraction and referencing using an array of indices assuming column-major storage
  /**
      \param RowIndx is a DMatrix array that contains row index values usually generated using the colon() function.
      \return Reference to a DMatrix mtype 1 object that maps to the referenced elements of the calling object
  */
   DMatrix& operator() (const DMatrix& RowIndx );
  //! Sub-vector extraction and referencing using an array of row indices for a given column
  /**
      \param RowIndx is a DMatrix array that contains row index values usually generated using the colon() function.
      \param col     is a column index
      \return Reference to a DMatrix mtype 1 object that maps to the referenced elements of the calling object
  */
   DMatrix& operator() (const DMatrix& RowIndx, long col );
  //! Sub-vector extraction and referencing using an array of column indices for a given row
  /**
      \param row is a row index.
      \param ColIndx is a DMatrix array that contains column index values usually generated using the colon() function.
      \return Reference to a DMatrix mtype 1 object that maps to the referenced elements of the calling object
  */
   DMatrix& operator() (long row, const DMatrix& ColIndx);
  //! Sub-vector extraction and referencing using an array of row indices for the last column of a matrix.
  /**
      \param RowIndx is a DMatrix array that contains row index values usually generated using the colon() function.
      \param end     is a character string containing the word "end".
      \return Reference to a DMatrix mtype 1 object that maps to the referenced elements of the calling object
  */
   DMatrix& operator() (const DMatrix& RowIndx,const char* end );
  //! Sub-vector extraction and referencing using an array of column indices for the last column of a matrix
  /**
      \param end is a character string containing the word "end".
      \param ColIndx is a DMatrix array that contains column index values usually generated using the colon() function.
      \return Reference to a DMatrix mtype 1 object that maps to the referenced elements of the calling object
  */
   DMatrix& operator() (const char* end, const DMatrix& ColIndx);

   // Miscellaneous friend functions
  //! This function generates a DMatrix object with a vector starting from a given value, with given increments and ending in a given value.
  /**
      \param i1 is the first value
      \param increment is the increment
      \param i2 is the last value
      \return Reference to a temporary DMatrix object
  */
   friend DMatrix& colon( double i1, double increment, double i2 );
  //! This function generates DMatrix object with a vector starting from a given value, with given increments and ending in a given value.
  /**
      \param i1 is the first value
      \param increment is the increment
      \param i2 is the last value
      \return Reference to a temporary DMatrix object
  */
   friend DMatrix& colon( int i1, int increment, int i2);
  //! This function generates DMatrix object with a vector starting from a given value, with unit increments, and ending in a given value.
  /**
      \param i1 is the first value
      \param i2 is the last value
      \return Reference to a temporary DMatrix object
  */
   friend DMatrix& colon( int i1, int i2 );
  //! This function generates DMatrix object with a vector starting from a given value, with unit increments, and ending in a given value.
  /**
      \param i1 is the first value
      \param i2 is the last value
      \return Reference to a temporary DMatrix object
  */
   friend DMatrix& colon( double i1, double i2);
  //! This function generates a special DMatrix object with one row and one column which is understood by the indexing functions that take a DMatrix object as an argument to mean "all rows" or "all columns".
  /**
      \return Reference to a temporary DMatrix object
  */
   friend DMatrix& colon( void );

   // Miscellaneous friend functions
  //! This function returns a 1 if any alement of DMatrix object that is passed as argument is non-zero, otherwise it returns a zero.
  /**
      \param  A is a DMatrix object.
      \return an integer value which is either 1 or 0.
  */
   friend int any( const DMatrix& A );
  //! This function calculates the integer matrix power.
  /**
      \param  A is a DMatrix object.
      \param  p is an integer value which can be positive or negative
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& mpow( DMatrix& A, int p );
  //! This function multiplies a real number by a matrix
  /**
      \param  r is a double value
      \param  A is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& operator *(double r, const DMatrix& A);
  //! This function returns the transpose of a given matrix.
  /**
      \param  A is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& tra(const DMatrix& A);
  //! This function returns the inverse of a given square matrix. If the argument is not a square matrix an error is thrown.
  /**
      \param  A is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& inv(const DMatrix& A);
  //! This function returns the pseudo-inverse of a given rectangular matrix.
  /**
      \param  A is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& pinv(const DMatrix& A);
  //! This function returns the identity matrix with a given number of rows and columns.
  /**
      \param  n is the desired number of rows and columns
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& identity(long n);
  //! This function returns a truncated identity matrix with specified numbers of rows and columns.
  /**
      \param  n is the desired number of rows
      \param  m is the desired number of columns
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& identity(long n, long m);
  //! This function returns the identity matrix with a given number of rows and columns.
  /**
      \param  n is the desired number of rows and columns
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& eye(long n);
  //! This function returns a truncated identity matrix with specified numbers of rows and columns.
  /**
      \param  n is the desired number of rows
      \param  m is the desired number of columns
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& eye(long n,long m);
  //! This function returns a matrix full of zeros with specified numbers of rows and columns.
  /**
      \param  n is the desired number of rows
      \param  m is the desired number of columns
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& zeros(long n, long m);
  //! This function returns a matrix full of ones with specified numbers of rows and columns.
  /**
      \param  n is the desired number of rows
      \param  m is the desired number of columns
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& ones(long n, long m);
  //! This function returns the exponential matrix of a given square matrix
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& expm(const DMatrix& A);
  //! This function returns a matrix with the sine of each element of the input matrix.
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& sin(const DMatrix& A);
  //! This function returns a matrix with the cosine of each element of the input matrix.
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& cos(const DMatrix& A);
  //! This function returns a matrix with the tangent of each element of the input matrix.
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& tan(const DMatrix& A);
  //! This function returns a matrix with the natural exponential of each element of the input matrix.
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& exp(const DMatrix& A);
  //! This function returns a matrix with the hyperbolic sine of each element of the input matrix.
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& sinh(const DMatrix& A);
  //! This function returns a matrix with the hyperbolic cosine of each element of the input matrix.
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& cosh(const DMatrix& A);
  //! This function returns a matrix with the hyperbolic tangent of each element of the input matrix.
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& tanh(const DMatrix& A);
  //! This function returns a matrix with the natural logarithm of each element of the input matrix.
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& log(const DMatrix& A);
  //! if A is a matrix this function extracts a column vector with the diagonal values of A. If A is a vector this function returns a matrix having the elements of A in the diagonal
  /**
      \param  A is a DMatrix object
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& diag( const DMatrix& A );
  //! This function returns the product of the first matrix transposed times the second matrix. The number of rows of both matrices must be the same, otherwise an error is thrown.
  /**
      \param  A is a DMatrix object
      \param  B is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& TProduct(const DMatrix& A,const DMatrix& B);
  //! This function returns the product of the first matrix times the second matrix transposed. The number of columns of both matrices must be the same, otherwise an error is thrown.
  /**
      \param  A is a DMatrix object
      \param  B is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& ProductT(const DMatrix& A,const DMatrix& B);
  //! This function returns the product of the first matrix transposed times the second matrix transposed. The number of rows of the first matrix must be the same as the number of columns of the second matrix, otherwise an error is thrown.
  /**
      \param  A is a DMatrix object
      \param  B is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& TProductT(const DMatrix& A,const DMatrix& B);
  //! This function calculates the product of two matrices. The number of columns of the first matrix must be the same as the number of rows of the second matrix, otherwise an error is thrown.
  /**
      \param  A is a DMatrix object
      \param  B is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& Product(const DMatrix& A,const DMatrix& B);
  //! Solves the system of equations \f$ A x = b \f$ using LU factorisation.
  /**
      \param  A is a DMatrix object
      \param  b is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation (the resulting vector x)
      \sa     LUFSolve(), LU()
  */
   friend DMatrix& LUSolve(const DMatrix& A ,const DMatrix& b );
  //! Solves the system of equations \f$ A x = b \f$ using LU factorisation using a previously found LU factors.
  /**
      \param  ALU is a DMatrix object resulting from a previous call to the LU() function
      \param  b is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation (the resulting vector x)
      \sa LU(), LUSolve()
  */
   friend DMatrix& LUFSolve(const DMatrix& ALU ,const DMatrix& b );
  //! Solves the system of equations \f$ A x = b \f$ using Cholesky factorisation.
  /**
      \param  A is a DMatrix object
      \param  b is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation (the resulting vector x)
      \sa     CholFSolve(), Chol()
  */
   friend DMatrix& CholSolve( const DMatrix& A, const DMatrix &b );
  //! Solves the system of equations \f$ A x = b \f$ using Cholesky factorisation. The function uses a previously found Cholesky factorisation.
  /**
      \param  Achol is a DMatrix object resulting from a previous call to Chol().
      \param  b is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation (the resulting vector x)
      \sa     CholFSolve(), Chol()
  */
   friend DMatrix& CholFSolve( const DMatrix& Achol, const DMatrix &b );
  //! Returns the Cholesky factorisation of a matrix A, which must be a positive definite symmetric matrix.
  /**
      \param  A is a DMatrix object·
      \return Reference to a temporary DMatrix object with the result of the operation
      \sa     CholFSolve(), CholSolve()
  */
   friend DMatrix& Chol( const DMatrix& A  );
  //! Returns the Cholesky root R of a given matrix A, such that A=R'R. A must be a positive definite symmetric matrix.
  /**
      \param  A is a DMatrix object·
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& CholeskyRoot(const DMatrix& A );
  //! Solves the system of equations \f$ A x = b \f$ using QR factorisation. The number of rows of matrix A must be greater or equal than the number of columns.
  /**
      \param  A is a DMatrix object.
      \param  b is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation (the resulting vector x)
      \sa     QRFSolve(), QR()
  */
   friend DMatrix& QRSolve( const DMatrix& A, const DMatrix &b );
  //! Solves the system of equations \f$ A x = b \f$ using QR factorisation. The function uses a previously found QR factorisation.
  /**
      \param  A is a DMatrix object.
      \param  b is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation (the resulting vector x)
      \sa     QRSolve(), QR()
  */
   friend DMatrix& QRFSolve( const DMatrix& A, const DMatrix &b );
  //! Returns the QR factorisation of a matrix A. The number of rows of matrix A must be greater or equal than the number of columns.
  /**
      \param  A is a DMatrix object·
      \return Reference to a temporary DMatrix object with the result of the operation
      \sa     QRFSolve(), QRSolve()
  */
   friend DMatrix& QR( const DMatrix& A );

  //! This function solves overdetermined  or  underdetermined  real  linear  systems \f$ A x = B \f$ using a QR or LQ factorization of A.  It is assumed that matrix A has full rank. The function uses the LAPACK routine dgels.
  /**
      \param  A is a DMatrix object
      \param  B is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation (the resulting vector x)
  */
   friend DMatrix& LSMNSolve( const DMatrix& A, const DMatrix& B );

  //! Returns the LQ factorisation of a matrix A. The function uses the LAPACK routine dgelqf_().
  /**
      \param  A is a DMatrix object·
      \param  Q is a pointer to a DMatrix object, which is modified on output to contain the Q factor of the decomposition.
      \return Reference to a temporary DMatrix object with the L factor of the decomposition.
  */
   friend DMatrix& LQ( const DMatrix& A, DMatrix* Q );
  //! Returns the LU factorisation of a matrix A.
  /**
      \param  A is a DMatrix object·
      \return Reference to a temporary DMatrix object with the result of the operation
      \sa     LUFSolve(), LUSolve()
  */
   friend DMatrix& LU( const DMatrix& A  );
  //! Returns the singular value decomposition of a matrix \f$ A = U' diag(s) V \f$, where vector s contains the singular values of A. The function uses the LAPACK routine dgesvd_().
  /**
      \param  A is a DMatrix object·
      \param  U is a pointer to a DMatrix object, which is modified on output to contain the U factor of the decomposition.
      \param  V is a pointer to a DMatrix object, which is modified on output to contain the V factor of the decomposition.
      \return Reference to a temporary DMatrix object with a vector that contains the singular values of matrix A.
  */
   friend DMatrix& SVD( const DMatrix& A, DMatrix* U=NULL, DMatrix* V=NULL );
   //! This function returns Q, the orthonormal basis for the range of a matrix A, such that \f$ Q Q' = I \f$. The number of columns of Q is the rank of A.
   /**
   	\param A is a DMatrix object
   	\return Reference to a temporary DMatrix object with the result of the operation
   */
   friend DMatrix& orth( const DMatrix& A );
   //! This function returns Z, the orthonormal basis for the null space of a matrix A, such that \f$ Z Z' = I \f$ and \f$ A Z=0 \f$. The number of columns of Z is the nullity of A.
   /**
   	\param A is a DMatrix object
   	\return Reference to a temporary DMatrix object with the result of the operation
   */
   friend DMatrix& null( const DMatrix& A );
   //!  This function uses the LAPACK routine dgelss_() to compute the minimum norm solution to a real linear least squares problem: Minimize \f$ || B - A x ||_2 \f$ using the singular value decomposition (SVD) of A. A is a rectangular matrix which may be rank-deficient.
  /**
      \param  A is a DMatrix object.
      \param  B is a DMatrix object.
      \return Reference to a temporary DMatrix object with the result of the operation (the resulting vector x)
      \sa     SVD()
  */
   friend DMatrix& SVDSolve( const DMatrix& A, const DMatrix& B );
   //!  This function computes and returns the Schur decomposition of a matrix A, such that \f$ A=Q'U Q \f$, where U is an upper triangular matrix and Q is a unitary matrix. This function uses the LAPACK routine dgees_()
  /**
      \param  A is a DMatrix object.
      \param  U is a pointer to a DMatrix object.
      \return Reference to a temporary DMatrix object with the unitary matrix Q.
  */
   friend DMatrix& schur(const DMatrix& A, DMatrix* U= NULL  );
   //!  This function computes the eigenvalues and (optionally) the eigenvectors of a matrix A. This function uses the LAPACK routines dsyev_() and dgeev_.
  /**
      \param  A is a DMatrix object.
      \param  V is a pointer to a DMatrix object.
      \return Reference to a temporary DMatrix object with the real part of the eigenvalues in the first column and the complex part of the eigenvalues in the second column.
  */
   friend DMatrix& eig(const DMatrix& A, DMatrix* V= NULL  );
   //!  This function computes and return the Euclidean norm of a matrix A, which is the square root of the sum of its squared elements.
  /**
      \param  A is a DMatrix object.
      \return the value of the Euclidean norm.
  */
   friend double enorm(const DMatrix& A);
   //!  This function computes 2-norm of matrix A, which is computed as the maximum singular value of A.
  /**
      \param  A is a DMatrix object.
      \return the value of the 2-norm
  */
   friend double norm(const DMatrix& A);
   //!  This function computes infinity norm of matrix A, which is computed as the maximum absolute value row sum.
  /**
      \param  A is a DMatrix object.
      \return the value of the infinity norm
  */
   friend double InfNorm(const DMatrix& A);
   //!  This function computes Frobenius norm of matrix A.
  /**
      \param  A is a DMatrix object.
      \return the value of the Frobenius norm
  */
   friend double Fnorm(  const DMatrix& A);
   //!  This function computes and returns the element-wise absolute value of matrix A.
  /**
      \param  A is a DMatrix object.
      \return a reference to a temporary DMatrix object with the result of the operation.
  */
   friend DMatrix& Abs(const DMatrix& A);
   //!  This function finds and returns the element of matrix A with maximum value. It also returns the indices of such element. If more than one element has the same maximum value, the indices of the first element found when searching column by column is returned.
  /**
      \param  A is a DMatrix object.
      \param  rindx is an optional pointer to an integer which is modified with the row index.
      \param  cindx is an optional pointer to an integer which is modified with the column index.
      \return the value of the element with maximum value.
  */
   friend double Max(const DMatrix& A,int* rindx=NULL, int* cindx=NULL);
   //!  This function finds and returns the element of matrix A with maximum absolute value. It also returns the indices of such element. If more than one element has the same maximum absolute value, the indices of the first element found when searching column by column is returned.
  /**
      \param  A is a DMatrix object.
      \param  rindx is a pointer to an integer which is modified with the row index.
      \param  cindx is a pointer to an integer which is modified with the column index.
      \return the absolute value of the element with maximum absolute value.
  */
   friend double MaxAbs(const DMatrix& A, int* rindx=NULL, int* cindx=NULL);
   //!  This function finds and returns the element of matrix A with minimum value. It also returns the indices of such element.  If more than one element has the same minimum value, the indices of the first element found when searching column by column is returned.
  /**
      \param  A is a DMatrix object.
      \param  rindx is a pointer to an integer which is modified with the row index.
      \param  cindx is a pointer to an integer which is modified with the column index.
      \return the absolute value of the element with minimum absolute value.
  */
   friend double Min(const DMatrix& A, int* rindx=NULL, int* cindx=NULL );
   //!  This function finds and returns the element of matrix A with minimum absolute value. It also returns the indices of such element.  If more than one element has the same minimum absolute value, the indices of the first element found when searching column by column is returned.
  /**
      \param  A is a DMatrix object.
      \param  rindx is a pointer to an integer which is modified with the row index.
      \param  cindx is a pointer to an integer which is modified with the column index.
      \return the absolute value of the element with minimum absolute value.
  */
   friend double MinAbs(const DMatrix& A, int* rindx=NULL, int* cindx=NULL);
   //!  This function sorts the input vector x in ascending order. Optionally, it also returns an integer array of sorted indices. If the input object is not a vector, then an error is thrown.
  /**
      \param  x is a DMatrix object which upon input contains the unsorted vector and upon output contains the sorted vector.
      \param  indx is a pointer to the first element of the array of sorted indices.
      \return void
  */
   friend void   sort( DMatrix& x, int indx[] = NULL);
   //!  This function sorts the input vector x in ascending order. It also returns a DMatrix object with the sorted indices. If the input object is not a vector, then an error is thrown.
  /**
      \param  x is a DMatrix object which upon input contains the unsorted vector and upon output contains the sorted vector.
      \param  indx is a DMatrix object which upon output contains the values of the sorted indices.
      \return void
  */
   friend void   sort( DMatrix& x, DMatrix& indx);
   //!  This function computes the dot product of two vectors. If any of the input arguments does not contain a vector, or if the vector lengths are not equal, an error is thrown.
  /**
      \param  x is a DMatrix object.
      \param  y is a DMatrix object
      \return the value of the dot product of x and y.
  */
   friend double dotProduct( const DMatrix& x, const DMatrix& y );
   //!  This function computes the dot product of two vectors. If any of the input arguments does not contain a vector, or if the vector lengths are not equal, an error is thrown.
  /**
      \param  x is a DMatrix object.
      \param  y is a DMatrix object
      \return the value of the dot product of x and y.
  */
   friend double dot(const DMatrix& x, const DMatrix& y) { return dotProduct(x,y); }
   //!  This function computes the cross product of two vectors. If any of the input arguments does not contain a vector, or if the  length of any of the vectors is not 3, an error is thrown.
  /**
      \param  x is a DMatrix object.
      \param  y is a DMatrix object
      \return the value of the cross product of x and y.
  */
   friend DMatrix& crossProduct(const DMatrix& x, const DMatrix& y);
   //!  This function computes the cross product of two vectors. If any of the input arguments does not contain a vector, or if the  length of any of the vectors is not 3, an error is thrown.
  /**
      \param  x is a DMatrix object.
      \param  y is a DMatrix object
      \return the value of the cross product of x and y.
  */
   friend DMatrix& cross(const DMatrix& x, const DMatrix& y) { return crossProduct(x,y); }
   //!  This function checks if the input matrix is symmetric. If the input matrix is not square, an error is thrown.
  /**
      \param  A is a DMatrix object.
      \return 1 if the input matrix is symmetric, 0 otherwise.
  */
   friend int isSymmetric( const DMatrix& A );
   //!  This function calculates the 2-norm condition number of a matrix, which is the ratio of the maximum singular value to the minimum singular value of the matrix.  A large condition number indicates a nearly singular matrix.  If the input matrix is not square, an error is thrown.
  /**
      \param  A is a DMatrix object.
      \return the 2-norm condition number
  */
   friend double cond( const DMatrix& A );
   //!  This function estimates the 1-norm reciprocal condition number of a matrix. The function uses the LAPACK function dgecon. if A is well conditioned, then rcond(A) is near 1. If A is badly conditioned, then rcond(A) is close to the machine numerical precision (very small).  If the input matrix is not square, an error is thrown.
  /**
      \param  A is a DMatrix object.
      \return the reciprocal condition number estimate
  */
   friend double rcond( const DMatrix& A );
   //!  This function returns an estimate of the rank of a matrix, which is the number of linearly independent rows or columns.
  /**
      \param  A is a DMatrix object.
      \return the rank estimate.
  */
   friend int rank_matrix( const DMatrix& A );
   //!  This function returns the determinant of a square matrix. If the input matrix is not square, an error is thrown.
  /**
      \param  A is a DMatrix object.
      \return the determinant of the matrix.
  */
   friend double det(  const DMatrix& A );
   //!  This function returns the trace of a square matrix. If the input matrix is not square, an error is thrown.
  /**
      \param  A is a DMatrix object.
      \return the trace of the matrix.
  */
   friend double trace( const DMatrix& A );
   //!  This function returns a row vector with the mean values of the columns of matrix A.
  /**
      \param  A is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& mean( const DMatrix& A );
   //!  This function returns a row vector with the standard deviation of each column of matrix A. If ntype is 0 (default) the result is normalised with (n-1), where n is the number of rows of A. Otherwise, the result is normalised with n.
  /**
      \param  A is a DMatrix object.
      \param  ntype is the type of normalization, 0 (default) or 1.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& Std(  const DMatrix& A, int ntype=0 );
   //!  Computes the covariance matrix of a data matrix where the N rows correspond to samples and the M columns are variables. The result is returned as an M x M matrix. If ntype=0 (default) then the result is normalised with N-1. Otherwise, if ntype=1, the  result is normalised with N.
  /**
      \param  A is a DMatrix object.
      \param  ntype is the type of normalization, 0 (default) or 1.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& cov( const DMatrix& A, int ntype=0 );
   //!  Computes the covariance matrix of two vectors X and Y of dimension N. The result is returned as an 1 x 1 DMatrix object. If ntype=0 (default) then the result is normalised with N-1. Otherwise, if ntype=1, the  result is normalised with N.
  /**
      \param  X is a DMatrix object.
      \param  Y is a DMatrix object
      \param  ntype is the type of normalization, 0 (default) or 1.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& cov(DMatrix& X, DMatrix& Y, int ntype=0 );
   //!  This function returns a row vector with the variance of each column of matrix A. If ntype is 0 (default) the result is normalised with (n-1), where n is the number of rows of A. Otherwise, the result is normalised with n.
  /**
      \param  A is a DMatrix object.
      \param  ntype is the type of normalization, 0 (default) or 1.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& var(DMatrix& A, int ntype=0 );
   //!  This function returns a row vector with the sum of the elements of each column of matrix A.
  /**
      \param  A is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& sum(  const DMatrix& A );
   //!  This function returns a row vector with the product of the elements of each column of matrix A.
  /**
      \param  A is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& prod(  const DMatrix& A );
   //!  This function computes and returns the element-wise product of two matrices of the same dimensions. If the dimensions of the two input matrices are not the same, an error is thrown.
  /**
      \param  A is a DMatrix object.
      \param B is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& elemProduct( const DMatrix& A, const DMatrix& B );
   //!  This function computes and returns the element-wise division of two matrices of the same dimensions. If the dimensions of the two input matrices are not the same, an error is thrown. The dimensions of the returned object are the same as the dimensions of the factors.
  /**
      \param  A is a DMatrix object.
      \param B is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& elemDivision( const DMatrix& A, const DMatrix& B);
   //!  This function computes and returns the Kronecker product of two matrices. The row (column) dimension of the returned object is the product of the row (column) dimensions of both factors.
  /**
      \param  A is a DMatrix object.
      \param B is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& kronProduct( const DMatrix& A, const DMatrix& B );
   //!  This function returns a column vector made by stacking the columns of a matrix one below the other from left to right.
  /**
      \param  A is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& vec( const DMatrix& A );
   //!  This function returns a DMatrix object with the same dimensions as the input matrix such that each of its elements is 1 is the corresponding value of the input matrix is positive, -1 if the corresponding value of the input matrix is negative, and 0 if the corresponding value of the input matrix is 0.
  /**
      \param  A is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& MatrixSign( const DMatrix& A );
   //! This function returns a column vector with the linear indices of the non-zero elements of the input matrix A. The linear index is 1 for element (1,1) of the input matrix A, and length(A) for the (nrows,ncols) element of the input matrix A.
  /**
      \param  A is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix&  find(const DMatrix& A);
   //! This function returns an nxm matrix where each element is a uniform pseudo-random number in the range (0,1).
  /**
      \param  n is the desired number of rows of the returned matrix
      \param  m is the desired number of columns of the returned matrix
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& randu(long n, long m);
   //! This function returns an nxm matrix where each element is a Gaussian pseudo-random number in the range with zero mean and variance 1.
  /**
      \param  n is the desired number of rows of the returned matrix
      \param  m is the desired number of columns of the returned matrix
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& randn(long n, long m);
   //! This function returns a linearly spaced vector with N points between the values X1 and X2.
  /**
      \param  X1 is a real number
      \param  X2 is a real number
      \param N is the desired length of the returned vector
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& linspace(double X1, double X2, long N);
   //! This function returns the machine numerical precision
   static double GetEPS() { return MC_EPSILON ; }

  //! This function prints an error message and throws an exception to be handled by the ErrorHandler class.
  /**
      \param  input_text is the error message
      \return void
  */
   friend void error_message(const char *input_text);
  //! This function, which is to be used in conjunction with function toc(), starts counting elapsed CPU time.
   friend void tic(void);
  //! This function, which is to be used in conjunction with function tic(), stops counting CPU time, and it prints and returns the elapsed time in seconds since the function tic() was called.
  /**
      \return the elapsed time in seconds.
  */
   friend double toc();

#ifdef MATLAB_MEX_FILE
   friend void mxArray2DMatrix( DMatrix& A, mxArray* mx );

   friend void DMatrix2mxArray( char* Name, mxArray* mx, DMatrix& A );

#endif
   //! This function computes the square root of each element of the input matrix A, and it returns a DMatrix object with the same dimensions as the input matrix. If any element of the input matrix is negative, an error is thrown.
  /**
      \param  A is a DMatrix object.
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& Sqrt(const DMatrix& A);
   //! This function extracts and return the triangular upper part of the input matrix A. The returned object has the same dimensions as the input object.
  /**
      \param  A is a DMatrix object
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& triu(const DMatrix& A);
   //! This function returns the N-by-M matrix whose elements are taken columnwise from the input matrix A.  An error is thrown if A does not have N*M elements.
  /**
      \param  A is a DMatrix object
      \param N is the desired number of rows of the reshaped matrix
      \param M is the desired number of columns of the reshaped matrix
      \return a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& reshape(DMatrix& A, long N, long M);
   //! This function returns the number of elements of a matrix A.
   /**
      \param A is a DMatrix object
      \return the number of elements of the input matrix.
   */
   friend long length(const DMatrix& A);

};


//! InitializeDMatrixClass class
/**
   This is a dummy C++ class intended to initialise the temporary objects of the DMatrix class.
*/

class InitializeDMatrixClass{

    public:
   //! This is the default constructor which calls the function DMatrix::AllocateAuxArr().
    InitializeDMatrixClass() { DMatrix::AllocateAuxArr(); }
   //! This is the destructor which calls the function DMatrix::DeAllocateAuxArr().
    ~InitializeDMatrixClass() { DMatrix::DeAllocateAuxArr(); }

};



/* inline utility functions and prototypes */
// Swap double
inline void Swap(double *x,double *y)
{double temp = *x; *x = *y; *y = temp;}
// Swap int
inline void Swap(int *x,int *y)
{int temp = *x; *x = *y; *y = temp;}
double MachEps(void);
double ArrayNorm( long n, double* x );
double factorial(long j);
/* Declaration of all friend functions of DMatrix class  */
void CholeskyDecomp(DMatrix& a, int n, DMatrix& pM);
void CholeskySolution(const DMatrix& a, int n, const DMatrix& pM,
                         const DMatrix& bM, DMatrix& xM);
void Hessemberg(DMatrix& a );
DMatrix& operator- (const DMatrix& A);
DMatrix& colon( double i1, double increment, double i2 );
DMatrix& colon( int i1, int increment, int i2);
DMatrix& colon( int i1, int i2 );
DMatrix& colon( double i1, double i2);
DMatrix& colon( void );
int any( const DMatrix& val );
DMatrix& mpow( DMatrix& A, int p );
DMatrix& operator *(double Arg, const DMatrix& A);
DMatrix& tra(const DMatrix& A);
DMatrix& inv(const DMatrix& A);
DMatrix& pinv(const DMatrix& A);
DMatrix& identity(long n);
DMatrix& identity(long n, long m);
DMatrix& eye(long n);
DMatrix& eye(long n,long m);
DMatrix& zeros(long n, long m);
DMatrix& ones(long n, long m);
DMatrix& expm(const DMatrix& A);
DMatrix& sin(const DMatrix& A);
DMatrix& cos(const DMatrix& A);
DMatrix& tan(const DMatrix& A);
DMatrix& exp(const DMatrix& A);
DMatrix& log(const DMatrix& A);
DMatrix& diag( const DMatrix& A );
DMatrix& TProduct(const DMatrix& A,const DMatrix& B);
DMatrix& ProductT(const DMatrix& A,const DMatrix& B);
DMatrix& TProductT(const DMatrix& A,const DMatrix& B);
DMatrix& Product(const DMatrix& A,const DMatrix& B);
DMatrix& LUSolve(const DMatrix& A ,const DMatrix& b );
DMatrix& LUFSolve(const DMatrix& ALU ,const DMatrix& b );
DMatrix& CholSolve( const DMatrix& A, const DMatrix &b );
DMatrix& CholFSolve( const DMatrix& Achol, const DMatrix &b );
DMatrix& Chol( const DMatrix& A  );
DMatrix& CholeskyRoot(const DMatrix& A );
DMatrix& QRSolve( const DMatrix& A, const DMatrix &b );
DMatrix& QRFSolve( const DMatrix& A, const DMatrix &b );
DMatrix& QR( const DMatrix& A );
DMatrix& LQ( const DMatrix& A, DMatrix* Q );
DMatrix& LU( const DMatrix& A  );
DMatrix& SVD( const DMatrix& A, DMatrix* U, DMatrix* V );
DMatrix& orth( const DMatrix& A );
DMatrix& null( const DMatrix& A );
DMatrix& SVDSolve( const DMatrix& A, const DMatrix& B );
DMatrix& schur(const DMatrix& A, DMatrix* U  );
DMatrix& eig(const DMatrix& A, DMatrix* V  );
double enorm(const DMatrix& A);
double norm(const DMatrix& A);
DMatrix& Abs(const DMatrix& A);
double Max(const DMatrix& A,int* rindx, int* cindx);
double MaxAbs(const DMatrix& A, int* rindx, int* cindx);
double Min(const DMatrix& A, int* rindx, int* cindx );
double MinAbs(const DMatrix& A, int* rindx, int* cindx);
void   sort( DMatrix& A, int indx[]);
void   sort( DMatrix& A, DMatrix& indx);
double InfNorm(const DMatrix& A);
double Fnorm(  const DMatrix& A);
double dotProduct( const DMatrix& A, const DMatrix& B );
int isSymmetric( const DMatrix& A );
double cond( const DMatrix& A );
double rcond( const DMatrix& A );
int rank_matrix(const DMatrix& A );
double det(  const DMatrix& A );
double trace( const DMatrix& A );
DMatrix& mean( const DMatrix& A );
DMatrix& Std(  const DMatrix& A, int ntype );
DMatrix& cov( const DMatrix& A, int ntype );
DMatrix& cov(DMatrix& X, DMatrix& Y, int ntype );
DMatrix& var(DMatrix& A, int ntype);
DMatrix& sum(  const DMatrix& A );
DMatrix& prod(  const DMatrix& A );
DMatrix& elemProduct( const DMatrix& A, const DMatrix& B );
DMatrix& elemDivision( const DMatrix& A, const DMatrix& B);
DMatrix& kronProduct( const DMatrix& A, const DMatrix& B );
DMatrix& vec( const DMatrix& A );
DMatrix& MatrixSign( const DMatrix& A );
DMatrix&  find(const DMatrix& A);
DMatrix& randu(long n, long m);
DMatrix& randn(long n, long m);
DMatrix& linspace(double X1, double X2, long N);
DMatrix& LSMNSolve( const DMatrix& A, const DMatrix& B );
void error_message(const char *input_text);
void tic(void);
double toc();
#ifdef MATLAB_MEX_FILE
void mxArray2DMatrix( DMatrix& A, mxArray* mx );
void DMatrix2mxArray( const char* Name, mxArray* mx, DMatrix& A );
#endif
DMatrix& Sqrt(const DMatrix& A);
DMatrix& triu(const DMatrix& A);
DMatrix& reshape(DMatrix& A, long nn, long mm);
long length(const DMatrix& A);
void* my_calloc(size_t num, size_t size );


// ===========================================================



#ifdef LAPACK



extern "C" {

#include "f2c.h"

/* CLAPACK real matrix eigenvalues routine */
int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *
	a, integer *lda, doublereal *wr, doublereal *wi, doublereal *vl,
	integer *ldvl, doublereal *vr, integer *ldvr, doublereal *work,
	integer *lwork, integer *info);

/* CLAPACK symmetric real matrix eigenvalues routine */
int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a,
	 integer *lda, doublereal *w, doublereal *work, integer *lwork,
	integer *info);

/* CLAPACK  real matrix singular value decomposition routine */
int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
	doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
	integer *info);

/* CLAPACK real matrix schur decomposition routine */
int dgees_(char *jobvs, char *sort, L_fp select, integer *n,
	doublereal *a, integer *lda, integer *sdim, doublereal *wr,
	doublereal *wi, doublereal *vs, integer *ldvs, doublereal *work,
	integer *lwork, logical *bwork, integer *info);

/* CLAPACK rectangular linear system solution routine */
int dgels_(char *trans, integer *m, integer *n, integer *
	nrhs, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	doublereal *work, integer *lwork, integer *info);

/* CLAPACK reciprocal condition number estimator      */
int dgecon_(char *norm, integer *n, doublereal *a, integer *
	lda, doublereal *anorm, doublereal *rcond, doublereal *work, integer *
	iwork, integer *info);

/* CLAPACK LU factorization routine */
int dgetrf_(integer *m, integer *n, doublereal *a, integer *
	lda, integer *ipiv, integer *info);

/* CLAPACK Minimum norm solution by using singular value decomp */
int dgelss_(integer *m, integer *n, integer *nrhs,
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	s, doublereal *rcond, integer *rank, doublereal *work, integer *lwork,
	 integer *info);

  /* CLAPACK LQ factorization routine */

int dgelqf_(integer *m, integer *n, doublereal *a, integer *
        lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);




} /* End of extern "C" block */

#endif /* LAPACK */


// #define SPARSE

#ifdef SPARSE_MATRIX


//! SparseMatrix class
/**
    A C++ class for sparse numerical linear algebra with interfaces to
    a number of CXSparse and LUSOL functions
*/
class SparseMatrix {
protected:
//! Array to store matrix elements
   double *a;
//! Number of matrix rows
   int n;
//! Number of matrix columns
   int m;
//! Number of allocated elements in a
   int asize;
//! Number of non-zero elements
   int nz;
//! Array of nonzero row indices
   int *RowIndx;
//! Array of nonzero column indices
   int *ColIndx;

public:

   // Interface functions
  //! Gets the number of rows from the calling object
  /**
      \return integer value with the number of rows.
  */
   int GetNoRows() const  { return n; }
  //! Gets the number of columns from the calling object
  /**
      \return integer value with the number of columns.
  */
   int GetNoCols() const { return m; }
  //! Gets the number of non-zero elements from the calling object
  /**
      \return integer value with the number of non-zero elements.
  */
   int GetNonZero() const {return nz;}
  //! Returns a DMatrix object with the row indices for the non-zero elements
  /**
      \return Temporary DMatrix object with the row indices for the non-zero elements
  */
   DMatrix& GetRowIndxDMatrix();
  //! Returns a DMatrix object with the column indices for the non-zero elements
  /**
      \return Temporary DMatrix object with the column indices for the non-zero elements
  */
   DMatrix& GetColIndxDMatrix();
  //! Returns a int pointer to the start of an array of row indices for the non-zero elements
  /**
      \return int pointer to the start of the array
  */
   int* GetRowIndx_C_Array() const {return RowIndx;}
  //! Returns a int pointer to the start of an array of column indices for the non-zero elements
  /**
      \return int pointer to the start of the array
  */
   int* GetColIndx_C_Array() const {return ColIndx;}
  //! Returns a double pointer to the start of an array where the non-zero elements are stored
  /**
      \return double pointer to the start of the array
  */
   double *GetPr()  const  { return a; }
  //! Inserts a non-zero value at a specified location of the sparse matrix. The function re-allocates storage as required.
  /**
      \param i row index (starting from 1)
      \param j column index (starting from 1)
      \param val non-zero double value to be inserted.
      \return void
  */
   void InsertNonZero(int i, int j, double val);
  //! Saves a sparse matrix in triplet form. The first row of the saved file contains the number of rows, number of columns and the number of non-zeros of the matrix. Each subsequent row contains the row index, the column index and the corresponding nonzero value.
  /**
      \param fname name of the file to be created
      \return void
  */
   void Save(const char* fname) const;
  //! Saves a sparse sparsity pattern for the sparse matrix, such that the saved matrix contains only zeros and asterisks.
  /**
      \param fname name of the file to be created
      \return void
  */
   void SaveSparsityPattern(const char* fname) const;
  //! Loads a sparse matrix in triplet form. The first row of the specified file must contain the number of rows, number of columns and the number of non-zeros of the matrix. Each subsequent row must contain the row index, the column index and the corresponding nonzero value, separated by spaces. An error is thrown is the file cannot be opened.
  /**
      \param fname name of the file to be read
      \return void
  */
   void Load(const char* fname);
  //! This function transposes a sparse matrix and returns the transposed sparse matrix.
  /**
      \param  A  sparse matrix to be transposed
      \return a temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& tra(const SparseMatrix& A);
  //! Extracts a specified sub-matrix from a sparse matrix and returns a sparse matrix.
  /**
      \param  istart: start of row range
      \param  iend: end of row range
      \param  jstart: start of column range
      \param  jend: end of column range
      \return Temporary SparseMatrix object with the specified sub-matrix
  */
   SparseMatrix& sub_matrix(int istart, int iend, int jstart, int jend) const;
  //! Extracts a specified column from a sparse matrix and returns a DMatrix object.
  /**
      \param  j: column index
      \return Temporary DMatrix object with the specified column
  */
   DMatrix& Column(int j) const;
  //! Extracts a specified row from a sparse matrix and returns a DMatrix object.
  /**
      \param  i: row index
      \return Temporary DMatrix object with the specified row
  */
   DMatrix& Row(int i) const;
  //! Assigns the elements of a SparseMatrix object to a section of the calling object.
  /**
      \param  istart: start of row range
      \param  iend: end of row range
      \param  jstart: start of column range
      \param  jend: end of column range
      \param  B: SparseMatrix object whose element values are to be copied into the calling object
      \return void
  */
   void set_sub_matrix(const SparseMatrix& B, int istart, int iend, int jstart, int jend);

   // Display functions
  //! Prints a sparse matrix in triplet format.
  /**
      \param  text: label to identify the sparse matrix to be printed
      \return void
  */
   void Print(const char* text);

   // Constructors
  //! Default constructor. Creates an empty sparse matrix with zero rows and columns.
   SparseMatrix(void); // Default constructor
  //! Constructor using triplet arrays.
  /**
      \param  aa: array with double non-zero elements
      \param  nn: number of rows
      \param  mm: number of columns
      \param  nnz: number of non-zero elements
      \param  RowIndxArg: array of int values with the row indices of the nonzero elements (indices start from 1).
      \param  ColIndxArg: array of int values with the column indices of the nonzero elements (indices start from 1).
  */
   SparseMatrix(double* aa, int nn, int mm, int nnz, int* RowIndxArg, int* ColIndxArg); // Constructor using arrays
  //! Constructor using DMatrix object. Creates a sparse matrix with the same number of rows and columns as the argument, and copies the non-zero values of the argument into the created sparse matrix.
  /**
      \param  A: DMatrix object
  */
   SparseMatrix( const DMatrix& A ); // Constructor using DMatrix object
  //! Copy constructor. Creates a copy of the argument.
  /**
        \param  A: SparseMatrix object
  */
   SparseMatrix( const SparseMatrix& A); // copy constructor
  //! Constructor without element asignment. Creates a SparseMatrix object with storage for a specified number of non-zero values.
  /**
        \param  nn: number of rows.
        \param  mm: number of columns.
        \param  nnz: number of non-zero values.
  */
   SparseMatrix( int nn, int mm, int nnz); // Constructor without element assigment

   // Destructor
  //! Desctructor.
   ~SparseMatrix();

   // Other functions
  //! Resizes a SparseMatrix object
  /**
      \param  nnew: new number of rows.
      \param  mnew: new number of columns.
      \param  nznew: new number of non-zero values.
      \return void
  */
   void Resize(int nnew, int mnew, int nznew);
  //! Inverts a SparseMatrix object and returns a SparseMatrix object
  /**
      \param  A: SparseMatrix object to be inverted
      \return A temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& inv(SparseMatrix& A);
  //! Transposes a SparseMatrix object.
  /**
      \return void
  */
   void Transpose();
  //! Computes the product of two sparse matrices and returns a SparseMatrix object. The inner matrix dimensions must be consistent, otherwise an error is thrown.
  /**
      \param  A: SparseMatrix object
      \param  B: SparseMatrix object
      \return A temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& Product(const SparseMatrix& A,const  SparseMatrix& B);
  //! Computes the product of two sparse matrices and returns a SparseMatrix object. The left hand side matrix is transposed. The  matrix dimensions must be consistent, otherwise an error is thrown.
  /**
      \param  A: SparseMatrix object
      \param  B: SparseMatrix object
      \return A temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& TProduct(SparseMatrix& A,SparseMatrix& B);
  //! Computes the product of two sparse matrices and returns a SparseMatrix object. The right hand side matrix is transposed. The  matrix dimensions must be consistent, otherwise an error is thrown.
  /**
      \param  A: SparseMatrix object
      \param  B: SparseMatrix object
      \return A temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& ProductT(SparseMatrix& A, SparseMatrix& B);
   //!  This function computes and returns the element-wise product of two sparse matrices of the same dimensions. If the dimensions of the two input matrices are not the same, an error is thrown.
  /**
      \param  A is a SparseMatrix object.
      \param B is a SparseMatrix object.
      \return a temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& elemProduct(const SparseMatrix A, const SparseMatrix& B);
  //! Converts a DMatrix object into a SparseMatrix object.
  /**
      \param  A: DMatrix object
      \return A temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& sparse(DMatrix& A);
  //! Creates a sparse identity matrix of specified dimension.
  /**
      \param  n: number of rows and columns of the sparse matrix to be created.
      \return A temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& speye(int n);
  //! Creates a sparse matrix with ones and zeros based on the sparsity pattern of a sparse matrix.
  /**
      \param  s: SparseMatrix object
      \return A temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& spones(const SparseMatrix& s);
  //! Creates a SparseMatrix object given triplet information stored in the columns of a DMatrix object.
  /**
      \param  A: DMatrix object. The first column contains the row indices, the second column contains the column indices and the third column contains the non-zero values.
      \return A temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& spconvert(DMatrix& A);
  //! Creates a DMatrix object given a SparseMatrix object.
  /**
      \param  A: SparseMatrix object.
      \return A temporary DMatrix object with the result of the operation
  */
   friend DMatrix& full(const SparseMatrix& A);
   //! This function returns Z, the orthonormal basis for the null space of a sparse matrix A, such that \f$ Z Z' = I \f$ and \f$ A Z=0 \f$. The number of columns of Z is the nullity of A.
   /**
   	\param A is a SparseMatrix object
   	\return Reference to a temporary DMatrix object with the result of the operation
   */
   friend DMatrix& null( const SparseMatrix& A );
   //! Returns the singular value decomposition of a sparse matrix \f$ A = U' diag(s) V \f$, where vector s contains the singular values of A. The function uses the LAPACK routine dgesvd_().
  /**
      \param  A is a SparseMatrix object·
      \param  U is a pointer to a DMatrix object, which is modified on output to contain the U factor of the decomposition.
      \param  V is a pointer to a DMatrix object, which is modified on output to contain the V factor of the decomposition.
      \return Reference to a temporary DMatrix object with a vector that contains the singular values of matrix A.
  */
   friend DMatrix& SVD(const SparseMatrix& A, DMatrix* U=NULL, DMatrix* V=NULL );
  //! Returns the QR factorisation of a sparse matrix A. The number of rows of matrix A must be greater or equal than the number of columns.
  /**
      \param  A is a SparseMatrix object·
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   friend DMatrix& QR( const SparseMatrix& A );
  //! Returns the LQ factorisation of a matrix A. The function uses the LAPACK routine dgelqf_().
  /**
      \param  A is a SparseMatrix object·
      \param  Q is a pointer to a DMatrix object, which is modified on output to contain the Q factor of the decomposition.
      \return Reference to a temporary DMatrix object with the L factor of the decomposition.
  */
   friend DMatrix& LQ( const SparseMatrix& A, DMatrix* Q );
   //! This function returns Q, the orthonormal basis for the range of a sparse matrix A, such that \f$ Q Q' = I \f$. The number of columns of Q is the rank of A.
   /**
   	\param A is a SparseMatrix object
   	\return Reference to a temporary DMatrix object with the result of the operation
   */
   friend DMatrix& orth( const SparseMatrix& A );
   //!  This function computes and returns the Schur decomposition of a sparse matrix A, such that \f$ A=Q'U Q \f$, where U is an upper triangular matrix and Q is a unitary matrix. This function uses the LAPACK routine dgees_()
  /**
      \param  A is a SparseMatrix object.
      \param  U is a pointer to a DMatrix object.
      \return Reference to a temporary DMatrix object with the unitary matrix Q.
  */
   friend DMatrix& schur(const SparseMatrix& A, DMatrix* U= NULL  );
   //!  This function computes the eigenvalues and (optionally) the eigenvectors of a sparse matrix A. This function uses the LAPACK routines dsyev_() and dgeev_.
  /**
      \param  A is a SparseMatrix object.
      \param  V is a pointer to a DMatrix object.
      \return Reference to a temporary DMatrix object with the real part of the eigenvalues in the first column and the complex part of the eigenvalues in the second column.
  */
   friend DMatrix& eig(const SparseMatrix& A, DMatrix* V= NULL  );
   //!  This function computes and return the Euclidean norm of a sparse matrix A, which is the square root of the sum of its (nonzero) squared elements.
  /**
      \param  A is a SparseMatrix object.
      \return the value of the Euclidean norm.
  */
   friend double enorm(const SparseMatrix& A);
   //!  This function computes 2-norm of sparse matrix A, which is computed as the maximum singular value of A.
  /**
      \param  A is a SparseMatrix object.
      \return the value of the 2-norm
  */
   friend double norm(const SparseMatrix& A);
   //!  This function computes and returns the element-wise absolute value of a sparse matrix A.
  /**
      \param  A is a SparseMatrix object.
      \return a reference to a temporary SparseMatrix object with the result of the operation.
  */
   friend SparseMatrix& Abs(const SparseMatrix& A);
   //!  This function calculates the 2-norm condition number of a sparse matrix, which is the ratio of the maximum singular value to the minimum singular value of the matrix.  A large condition number indicates a nearly singular matrix.  If the input matrix is not square, an error is thrown.
  /**
      \param  A is a SparseMatrix object.
      \return the 2-norm condition number
  */
   friend double cond( const SparseMatrix& A );
   //!  This function estimates the 1-norm reciprocal condition number of a sparse matrix. The function uses the LAPACK function dgecon. if A is well conditioned, then rcond(A) is near 1. If A is badly conditioned, then rcond(A) is close to the machine numerical precision (very small).  If the input matrix is not square, an error is thrown.
  /**
      \param  A is a SparseMatrix object.
      \return the reciprocal condition number estimate
  */
   friend double rcond( const SparseMatrix& A );
   //!  This function returns an estimate of the rank of a sparse matrix, which is the number of linearly independent rows or columns.
  /**
      \param  A is a SparseMatrix object.
      \return the rank estimate.
  */
   friend int rank_sparse( const SparseMatrix& A );
   //! This function returns a sparse matrix with a given sparsity pattern where each non-zero element is a uniform pseudo-random number in the range (0,1).
  /**
      \param  S is a SparseMatrix object from which the sparsity pattern is taken.
      \return a temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& sprand(const SparseMatrix& S);
   //! This function returns a sparse matrix with a given density where each non-zero element is a uniform pseudo-random number in the range (0,1).
  /**
      \param  n is the number of rows
      \param  m is the number of columns
      \param  density is the desired density of the sparse matrix to be created.
      \return a temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& sprand(int n, int m, double density);
   //! This function returns a sparse matrix with a given sparsity pattern where each non-zero element is a Gaussian pseudo-random number in the range with zero mean and unit variance.
  /**
      \param  S is a SparseMatrix object from which the sparsity pattern is taken.
      \return a temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& sprandn(const SparseMatrix& S);
   //! This function returns a sparse matrix with a given density where each non-zero element is a Gaussian pseudo-random number with zero mean and unit variance.
  /**
      \param  n is the number of rows
      \param  m is the number of columns
      \param  density is the desired density of the sparse matrix to be created.
      \return a temporary SparseMatrix object with the result of the operation
  */
   friend SparseMatrix& sprandn(int n, int m, double density);

   //! This function eliminates zero elements from a sparse matrix and deletes unnecessary storage.
  /**
      \return void
  */
   void Compress();

   // Operators
  //! Sparse matrix addition and substitution operator. The sizes of the matrices being added must be the same, otherwise an error is thrown.The left hand side object is replaced with the result of the operation.
  /**
      \param rval:  SparseMatrix object located right hand side of the operator.
      \return Reference to the calling SparseMatrix object
  */
   SparseMatrix& operator += (const SparseMatrix &rval);
  //! Sparse matrix subtraction and substitution operator. The sizes of the matrices being subtracted must be the same, otherwise an error is thrown. The left hand side object is replaced with the result of the operation.
  /**
      \param rval:  SparseMatrix object located right hand side of the operator.
      \return Reference to the calling SparseMatrix object
  */
   SparseMatrix& operator -= (const SparseMatrix &rval);
  //! Sparse matrix product and substitution operator. The inner sizes of the matrices being multiplied must be the consistent, otherwise an error is thrown. The left hand side object is replaced with the result of the operation.
  /**
      \param rval:  SparseMatrix object located right hand side of the operator.
      \return Reference to the calling SparseMatrix object
  */
   SparseMatrix& operator *= (const SparseMatrix &rval);
  //! Computes the product of a sparse matrix (left hand side of the operator) times a real scalar (right hand side value) and replaces the left hand side object with the result of the operation.
  /**
      \param Arg: double value that will multiply each non-zero element of the sparse matrix.
      \return Reference the calling SparseMatrix object
  */
   SparseMatrix& operator*= (double Arg);
  //! Computes the division of a sparse matrix (left hand side of the operator) by a real scalar (right hand side value) and replaces the left hand side object with the result of the operation.
  /**
      \param Arg: double value that will divide each non-zero element of the sparse matrix.
      \return Reference the calling SparseMatrix object
  */
   SparseMatrix& operator/= (double Arg);
  //! Sparse matrix addition operator. The row and column sizes of the matrices being added must be the same, otherwise an error is thrown.
  /**
      \param rval:  sparse matrix located at the right hand side of the operator.
      \return Reference to a temporary SparseMatrix object with the result of the operation
  */
   SparseMatrix& operator+ (const SparseMatrix& rval) const;
  //! Sparse matrix subtraction operator. The row and column sizes of the matrices being subtracted must be the same, otherwise an error is thrown.
  /**
      \param rval:  sparse matrix located at the right hand side of the operator.
      \return Reference to a temporary SparseMatrix object with the result of the operation
  */
   SparseMatrix& operator - (const SparseMatrix& rval) const;
  //! Sparse matrix product operator. The inner dimensions of the matrices being multiplied must be the same, otherwise an error is thrown.
  /**
      \param rval:  sparse matrix located at the right hand side of the operator.
      \return Reference to a temporary SparseMatrix object with the result of the operation
  */
   SparseMatrix& operator* (const SparseMatrix& rval) const;
  //! This operator computes the product of a sparse matrix by a dense matrix and returns a sparse matrix object. The inner dimensions of the matrices being multiplied must be the same, otherwise an error is thrown.
  /**
      \param rval:  DMatrix object located at the right hand side of the operator.
      \return Reference to a temporary SparseMatrix object with the result of the operation
  */
   SparseMatrix& operator* (DMatrix& rval) const;
  //! Computes the product of a sparse matrix (left hand side of the operator) times a real scalar (right hand side value).
  /**
      \param Arg: double value that will multiply each non-zero element of the sparse matrix.
      \return Reference to a temporary object with the result of the operation.
  */
   SparseMatrix& operator* (double Arg) const;
  //! Computes the product of a real value (left hand side of the operator) by a sparse matrix (right hand side of the operator).
  /**
      \param Arg: double value that will multiply each non-zero element of the sparse matrix.
      \param A:  SparseMatrix object to be multiplied by a real value.
      \return Reference to a temporary object with the result of the operation.
  */
   friend SparseMatrix& operator *(double Arg, const SparseMatrix& A);
  //! Computes the division of a sparse matrix (left hand side of the operator) by a real scalar (right hand side value).
  /**
      \param Arg: double value that will divide each non-zero element of the sparse matrix.
      \return Reference to a temporary object with the result of the operation.
  */
   SparseMatrix& operator/ (double Arg) const;
  //! Computes the right division of a sparse matrix (left hand side of the operator) by another sparse matrix (right hand side value). This is conceptually equivalent to multiplying the left object by the inverse of the right hand side object but it is computed in a more efficient way. The dimensions of the matrices must be consistent, otherwise an error is returned. The right hand side object must be a square sparse matrix.
  /**
      \param rval: SparseMatrix object at the right hand side of the operator.
      \return Reference to a temporary SparseMatrix object with the result of the operation
  */
   SparseMatrix& operator/ (SparseMatrix& rval);
  //! Sparse matrix assignment. The size of the left hand side object is modified if necessary, and the values of all non-zero real elements of the right hand side object and copied to the left hand side object.
  /**
      \param rval: SparseMatrix object at the right hand side of the operator
      \return Reference to the calling SparseMatrix object
  */
   SparseMatrix& operator= (const SparseMatrix& rval);
  //! Elementwise power operator. Returns a SparseMatrix object with the same dimensions of the calling object and each of its elements is computed as the corresponding element of the calling object to the power of the right hand side argument. Care must be taken when using this operator as the associations do not work in the same way as with the * operator. It is highly recommended to use parenthesis every time this operator is used. For example use it as follows:   (A^x).
  /**
      \param x: double argument at the right hand side of the operator.
      \return Reference to a temporary DMatrix object with the result of the operation
  */
   SparseMatrix& operator^(double x) const;
  //! Concatenates two sparse matrices side by side. The number of rows of the matrices involved must be the same, otherwise an error is thrown. The number of columns of the resulting matrix is the addition of the number of columns of both matrices involved.
  /**
      \param B: SparseMatrix object at the right hand side of the operator.
      \return Reference to a temporary SparseMatrix object with the result of the operation
  */
   SparseMatrix& operator||(SparseMatrix& B) const;
  //! Stacks the right hand side matrix below the left hand side matrix. The dimensions number of columns of the matrices involved must be the same, otherwise an error is thrown. The number of rows of the resulting matrix is the addition of the number of rows of both matrices involved.
  /**
      \param B: SparseMatrix object at the right hand side of the operator.
      \return Reference to a temporary SparseMatrix object with the result of the operation
  */
   SparseMatrix& operator&&(SparseMatrix& B) const;
  //! Matrix indexing. Returns a reference to the matrix element located at the position indicated by the row and column indices. Indices start from 1. If the row and column indices refer to a zero (non-allocated) element, a new element with zero value is allocated so that it can be returned as a reference.
  /**
      \param row: Row index starting from 1.
      \param col: Column index starting from 1.
      \return Reference to the indexed matrix element.
  */
   double& operator() (int row,int col);
  //! Matrix indexing. Returns the value of the matrix element located at the position indicated by the row and column indices. Indices start from 1.
  /**
      \param row: Row index starting from 1.
      \param col: Column index starting from 1.
      \return double value of the element at the indicated position.
  */
   double operator() (int row,int col) const;
  //! Computes the left division of a sparse matrix (left hand side of the operator) by another sparse matrix (right hand side value). This is conceptually equivalent to multiplying the inverse of the left object by the right hand side object but it is computed in a more efficient way. The dimensions of the matrices must be consistent, otherwise an error is returned. The left hand side object must be a square matrix.
  /**
      \param B: SparseMatrix object at the right hand side of the operator.
      \return Reference to a temporary SparseMatrix object with the result of the operation
   */
   SparseMatrix& operator %(SparseMatrix& B);
  //! Solves a sparse system of equations \f$ Ax = b \f$ by using sparse LU decomposition. Matrix A must be a square matrix. This function uses LUSOL library functions.
  /**
      \param A: SparseMatrix object
      \param b: SparseMatrix object
      \param x: Reference to DMatrix object. Upon return this vector contains the solution to the system of equations.
      \param r: Reference to DMatrix object. Upon return this vector contains the residual of the solution, \f$ r=Ax-b \f$.
      \param argc: number of extra parameters (plus 1) to be passed to LUSOL. The default value is 1.
      \param argv: pointer to array of character strings with the extra parameters passed to LUSOL. The default is NULL. See the LUSOL documentation and source code for more details.
      \return void
   */
   friend void SparseLUSolve(SparseMatrix& A, DMatrix& b, DMatrix& x, DMatrix& r, int argc, char *argv[]);
  //! Finds the sparse LU factor of a sparse matrix A. Matrix A must be a square matrix. This function uses LUSOL library functions.
  /**
      \param A: SparseMatrix object
      \param argc: number of extra parameters (plus 1) to be passed to LUSOL. The default value is 1.
      \param argv: pointer to array of character strings with the extra parameters passed to LUSOL. The default is NULL. See the LUSOL documentation and source code for more details.
      \return void pointer.
      \sa SparseLUSolveGivenFactor()
   */
   friend void* SparseLUFactor(SparseMatrix& A, int argc, char *argv[]);
  //! Solves a sparse system of equations \f$ Ax = b \f$ by using sparse LU decomposition and a given LU factor. Matrix A must be a square matrix. This function uses LUSOL library functions.
  /**
      \param LUSOLv: void pointer resulting from a previous call to SparseLUFactor().
      \param A: SparseMatrix object
      \param b: SparseMatrix object
      \param x: Reference to DMatrix object. Upon return this vector contains the solution to the system of equations.
      \param r: Reference to DMatrix object. Upon return this vector contains the residual of the solution, \f$ r=Ax-b \f$.
      \param argc: number of extra parameters (plus 1) to be passed to LUSOL. The default value is 1.
      \param argv: pointer to array of character strings with the extra parameters passed to LUSOL. The default is NULL. See the LUSOL documentation and source code for more details.
      \return void
      \sa SparseLUFactor()
   */
   friend void SparseLUSolveGivenFactor(void* LUSOLv, SparseMatrix& A, DMatrix& b, DMatrix& x, DMatrix& r,  int argc, char *argv[]);
  //! This function prints an error message and throws an exception to be handled by the ErrorHandler class.
  /**
      \param  error_text is a character string with the error message
      \return void
  */
   friend void sp_error_message(const char *error_text);


};

// Declaration of all friend functions of SparseMatrix class

void SparseLUSolve(SparseMatrix& A, DMatrix& bb, DMatrix& xx, DMatrix& rr, int argc=1, char *argv[]=NULL);

void* SparseLUFactor(SparseMatrix& A, int argc=1, char *argv[]=NULL);

void SparseLUSolveGivenFactor(void* LUSOLv, SparseMatrix& A, DMatrix& bb, DMatrix& xx, DMatrix& rr,  int argc=1, char *argv[]=NULL);

void sp_error_message(const char *error_text);

SparseMatrix& TProduct(SparseMatrix& A,SparseMatrix& B);

SparseMatrix& ProductT(SparseMatrix& A, SparseMatrix& B);

SparseMatrix& elemProduct(const SparseMatrix A, const SparseMatrix& B);

SparseMatrix& tra(const SparseMatrix& A);
SparseMatrix& inv(SparseMatrix& A);
SparseMatrix& Product(const SparseMatrix& A,const  SparseMatrix& B);
SparseMatrix& sparse(DMatrix& A);
SparseMatrix& speye(int n);
SparseMatrix& spones(const SparseMatrix& s);
SparseMatrix& spconvert(DMatrix& A);
DMatrix& full(const SparseMatrix& A);
DMatrix& null( const SparseMatrix& A );
DMatrix& SVD(const SparseMatrix& A, DMatrix* U, DMatrix* V );
DMatrix& QR( const SparseMatrix& A );
DMatrix& LQ( const SparseMatrix& A, DMatrix* Q );
DMatrix& orth( const SparseMatrix& A );
DMatrix& schur(const SparseMatrix& A, DMatrix* U  );
DMatrix& eig(const SparseMatrix& A, DMatrix* V  );
double enorm(const SparseMatrix& A);
double norm(const SparseMatrix& A);
SparseMatrix& Abs(const SparseMatrix& A);
double cond( const SparseMatrix& A );
double rcond( const SparseMatrix& A );
int rank_sparse( const SparseMatrix& A );
SparseMatrix& sprand(const SparseMatrix& S);
SparseMatrix& sprand(int n, int m, double density);
SparseMatrix& sprandn(const SparseMatrix& S);
SparseMatrix& sprandn(int n, int m, double density);
SparseMatrix& operator *(double Arg, const SparseMatrix& A);


// ===========================================================




#endif /* SPARSE_MATRIX */


//! ErrorHandler class
/**
   This is a C++ class intended to handle error conditions.
*/
class ErrorHandler
{
    public:
        //! A string of characters which contains the error message
        string error_message;
        //! A constructor which takes the error message as an argument and assigns it to error_message
        /**
          \param m is the error message string.
          \sa function error_message().
        */
        ErrorHandler(const string m);
};



#endif /* DMATRIX_H */



