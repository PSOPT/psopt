//
//  SparseMatrix.hpp
//  dmatrix
//
//  Created by Philipp Waxweiler on 21.07.20.
//

#ifndef SparseMatrix_hpp
#define SparseMatrix_hpp

#include "dmatrix.h"
#include "cs.h"

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
  friend DMatrix& SVD(const SparseMatrix& A);
   friend DMatrix& SVD(const SparseMatrix& A, DMatrix* U, DMatrix* V );
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
  friend DMatrix& schur(const SparseMatrix& A);
   friend DMatrix& schur(const SparseMatrix& A, DMatrix* U  );
   //!  This function computes the eigenvalues and (optionally) the eigenvectors of a sparse matrix A. This function uses the LAPACK routines dsyev_() and dgeev_.
  /**
      \param  A is a SparseMatrix object.
      \param  V is a pointer to a DMatrix object.
      \return Reference to a temporary DMatrix object with the real part of the eigenvalues in the first column and the complex part of the eigenvalues in the second column.
  */
  friend DMatrix& eig(const SparseMatrix& A);
   friend DMatrix& eig(const SparseMatrix& A, DMatrix* V  );
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


#endif /* SparseMatrix_hpp */
