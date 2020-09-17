/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

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

#include "psopt.h"

// Implementation of TripletSparseMatrix functions


TripletSparseMatrix::TripletSparseMatrix(void)
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


TripletSparseMatrix::TripletSparseMatrix( int nn, int mm, int nzz)
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

TripletSparseMatrix::TripletSparseMatrix(double* aArg, int nArg, int mArg, int nzArg, int* RowIndxArg, int* ColIndxArg)
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


}



void TripletSparseMatrix::resize(int nnew, int mnew, int nznew)
{

    double *anew;
    int    *RowIndxNew;
    int    *ColIndxNew;
    int    dmax;

    anew = new double[nznew];
    RowIndxNew = new int[nznew];
    ColIndxNew = new int[nznew];

    dmax = std::min( nznew, nz );

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



TripletSparseMatrix::TripletSparseMatrix( const TripletSparseMatrix& A) // copy constructor
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
TripletSparseMatrix::~TripletSparseMatrix()
{
    if (a!=NULL)
 	delete a;
    if (RowIndx!= NULL)
        delete RowIndx;
    if (ColIndx!= NULL)
        delete ColIndx;

}




void TripletSparseMatrix::InsertNonZero(int i, int j, double val)
{
    int k;

    double* anew;
    int* RowNew;
    int* ColNew;
    int nznew=nz+1;
    int eflag = 0;

    for (k=0; k< nz; k++)
    {
        if (RowIndx[k]==i && ColIndx[k]==j) {
               a[k]= val;
               eflag = 1;
               break;
        }
    }

    if (!eflag) {
    	
      anew =  new double[nznew];
      RowNew= new int[nznew];
      ColNew= new int[nznew];

    	memcpy(anew,   a,       nz*sizeof(double) );
    	memcpy(RowNew, RowIndx, nz*sizeof(int) );
    	memcpy(ColNew, ColIndx, nz*sizeof(int) );

    	delete a;
    	delete RowIndx;
    	delete ColIndx;

    	a = anew;
    	RowIndx = RowNew;
    	ColIndx = ColNew;
    	nz = nznew;
    	asize = nznew;

    	RowIndx[nznew-1] = i;
      ColIndx[nznew-1] = j;
      a[nznew-1]       = val;
      
    }

}


// Operators

TripletSparseMatrix& TripletSparseMatrix::operator += (const TripletSparseMatrix &rval)
{
    int i;

    if ( (this->n != rval.n) || (this->m != rval.m) )
        fprintf(stderr,"\nTripletSparseMatrix::operator += error: matrix dimensions do not agree");

    for (i=0; i<rval.nz; i++)
    {
        if ( (*this)(rval.RowIndx[i],rval.ColIndx[i])!=0.0 )
               (*this)(rval.RowIndx[i],rval.ColIndx[i]) += rval.a[i];
        else
               this->InsertNonZero( rval.RowIndx[i], rval.ColIndx[i], rval.a[i] );
    }
    this->Compress();
    return (*this);
}


TripletSparseMatrix TripletSparseMatrix::operator+ (const TripletSparseMatrix& Other_matrix) const
{
   TripletSparseMatrix Result;

   Result=(*this);

   return ((Result)+=Other_matrix);

}

TripletSparseMatrix& TripletSparseMatrix::operator -= (const TripletSparseMatrix &rval)
{
    int i;

    if ( (this->n != rval.n) || (this->m != rval.m) )
        fprintf(stderr,"\nTripletSparseMatrix::operator -= error: matrix dimensions do not agree");

    for (i=0; i<rval.nz; i++)
    {
        if ( (*this)(rval.RowIndx[i],rval.ColIndx[i])!=0.0 )
             (*this)(rval.RowIndx[i],rval.ColIndx[i]) -= rval.a[i];
        else
               this->InsertNonZero( rval.RowIndx[i], rval.ColIndx[i], -rval.a[i] );
    }
    this->Compress();
    return (*this);
}


TripletSparseMatrix TripletSparseMatrix::operator- (const TripletSparseMatrix& Other_matrix) const
{

   TripletSparseMatrix Result;

   Result=(*this);

   return ((Result)-=Other_matrix);
}




TripletSparseMatrix& TripletSparseMatrix::operator*= (double Arg)
{
   int i;
   for (i=0; i<nz; i++)
   {
      a[i]*= Arg;
   }
   return (*this);
}


TripletSparseMatrix TripletSparseMatrix::operator* (double Arg) const
{
   TripletSparseMatrix Result;

   (Result)=(*this);

   return ( (Result)*=Arg );
}

TripletSparseMatrix operator *(double Arg, const TripletSparseMatrix& A)
{
   return (A*Arg);
}


TripletSparseMatrix TripletSparseMatrix::operator* (const MatrixXd& A) const
{
  int k;
  int i, j, icol, q;
  MatrixXd x;
  MatrixXd r(n,1);
  TripletSparseMatrix sp(n, A.cols(), 0);


  for (icol=0; icol<A.cols(); icol++) {
     x = A.col(icol);
     r.setZero();
     for (k = 0; k < nz; k++) {
         j    = ColIndx[k];
         i    = RowIndx[k];
         r(i) += a[k]*x(j);
     }
     for (q=0; q<length(r); q++) {
        sp.InsertNonZero(q,icol, r(q) );
     }
  }

  sp.Compress();

  return (sp);

}


TripletSparseMatrix elemProduct(const TripletSparseMatrix A, const TripletSparseMatrix& B)
{
          int i, j, k;

          double sum;

          double *a = A.GetPr();

          int *RowIndx = A.GetRowIndx_C_Array();
          int *ColIndx = A.GetColIndx_C_Array();

          TripletSparseMatrix sp( A.rows(), B.cols(), 0 );

          if (A.rows() != B.rows() ||  A.cols() != B.cols())
                sp_error_message("Inconsistent matrix dimensions in TripletSparseMatrix elemProduct()");



          for (k = 0; k < A.GetNonZero(); k++) {
                i = RowIndx[k];
                j = ColIndx[k];
                sum = a[k]*B(i,j);
                if (sum!=0.0) sp.InsertNonZero(i, j, sum);
          }

          return (sp);
}


TripletSparseMatrix& TripletSparseMatrix::operator/= (double Arg)
{
   int i;

   if (Arg==0.0) sp_error_message("Division by zero in TripletSparseMatrix::operator/");

   for (i=0; i<nz; i++)
   {
      a[i]/= Arg;
   }
   return (*this);
}

TripletSparseMatrix TripletSparseMatrix::operator/ (double Arg) const
{
   TripletSparseMatrix Result;

   if (Arg==0.0) sp_error_message("Division by zero in TripletSparseMatrix::operator/");



   Result=(*this);

   return ( (Result)/=Arg );
}

TripletSparseMatrix& TripletSparseMatrix::operator= (const TripletSparseMatrix& Other_matrix)
{
     this->resize(Other_matrix.n, Other_matrix.m, Other_matrix.nz);
     memcpy(a, Other_matrix.a, nz*sizeof(double) );
     memcpy(RowIndx, Other_matrix.RowIndx, nz*sizeof(int) );
     memcpy(ColIndx, Other_matrix.ColIndx, nz*sizeof(int) );
     return (*this);
}

double& TripletSparseMatrix::operator() (int i,int j)
{
   int k;

    int ii;

    int eflag = 0;
    
    if (i>=n || i<0 || j >=m || j<0 ) {

         sp_error_message("Index out of range in operator TripletSparseMatrix::operator(int, int)");    
    
    } 

    for (k=0; k< nz; k++)
    {
       if (RowIndx[k]==i && ColIndx[k]==j) {
              ii = k;
              eflag = 1;
              break;
       }

    }

//    if( !eflag ) {
//
//    	for (k=0; k<nz; k++) {
//
//        	  if (RowIndx[k]==0 || ColIndx[k]==0) {
//           		RowIndx[k]=i;
//           		ColIndx[k]=j;
//                        ii=k;
//                        eflag=1;
//                        break;
//       		  }
//        }
//
//    }

    if (!eflag) {
       this->InsertNonZero(i,j,0.0);
       ii = nz-1;
    }

    return a[ii];

}

double TripletSparseMatrix::operator() (int row,int col) const
{
      int i;
      double retval = 0.0;

      if ( (row>=n) || (row<0) || (col<0) || (col>=m) ) {
          sp_error_message("Out of range index in TripletSparseMatrix::operator()");
      }

      for (i=0; i< nz; i++) {
          if (RowIndx[i]==row && ColIndx[i]==col)
          {
               retval = a[i];
               break;
          }
      }

      return retval;
}



void TripletSparseMatrix::Print(const char* text)
{
   int i;

   fprintf(stderr,"\nSparse matrix %s",text);
   fprintf(stderr,"\nNumber of rows: %li", n);
   fprintf(stderr,"\nNumber of columns: %li", m);
   fprintf(stderr,"\nNumber of non-zero elements: %li", nz);
   if (n*m!=0) fprintf(stderr,"\nDensity: %lf%%", (  ((double) nz*100)/(double) (n*m)  ) );
   if (nz>0) fprintf(stderr,"\n(Row,Col)\tValue");

   for (i=0; i<nz; i++)
   {
       fprintf(stderr,"\n(%li,%li)\t\t%e", RowIndx[i], ColIndx[i], a[i]);

   }
   fprintf(stderr,"\n");

}

void TripletSparseMatrix::SaveSparsityPattern(const char* text) const
{
   int i, j;
   FILE *outfile;

   outfile = fopen(text,"w");

   for (i=0; i<n; i++)  {
       for(j=0; j<m; j++)  {
           if ( (*this)(i,j)!=0.0 )
               fprintf(outfile,"%c ",'*');
           else
               fprintf(outfile,"%c ",'0');
       }
       fprintf(outfile,"\n");
   }
   fclose(outfile);
}


void TripletSparseMatrix::Compress()
{

    int i, j;
    int deleted_count=0;

    double* anew;
    int* RowNew;
    int* ColNew;
    int nznew;
    int nztemp=nz;

    for (i=0; i< nztemp; i++)
    {
        if (fabs(a[i])<=(1.1*PSOPT_extras::GetEPS()) ) {
             for(j=i; j<nztemp-1; j++)
	     {
                    a[j] = a[j+1];
                    RowIndx[j]=RowIndx[j+1];
                    ColIndx[j]=ColIndx[j+1];
	     }
 	     deleted_count++;
             nztemp--;
             i--;
        }
    }

    nznew = nz-deleted_count;
    anew = new double[nznew];
    RowNew= new int[nznew];
    ColNew= new int[nznew];

    memcpy(anew,a, nznew*sizeof(double) );
    memcpy(RowNew, RowIndx, nznew*sizeof(int) );
    memcpy(ColNew, ColIndx, nznew*sizeof(int) );

    delete a;
    delete RowIndx;
    delete ColIndx;

    a = anew;
    RowIndx = RowNew;
    ColIndx = ColNew;

    nz = nznew;
    asize = nznew;

}



TripletSparseMatrix speye(int n)
{
// Returns a sparse identity matrix
     TripletSparseMatrix sp;
     int i;
     sp.resize(n, n, n);

     for (i=0;i<sp.GetNonZero();i++) {
         sp.GetRowIndx_C_Array()[i]=i;
         sp.GetColIndx_C_Array()[i]=i;
         sp.GetPr()[i]=1.0;
     }

     return (sp);
}


void TripletSparseMatrix::Load(const char* fname)
{
  int nrow, ncol, nnz;

  int k;

  FILE *fp;

  if ( (fp = fopen(fname,"r")) == NULL )

  {  sp_error_message( "Error opening file in TripletSparseMatrix::Load()"); }


   fscanf(fp,"%li", &nrow);
   fscanf(fp,"%li", &ncol);
   fscanf(fp,"%li", &nnz);

   this->resize(nrow, ncol, nnz);

   for (k=0;k<nnz;k++) {
            fscanf(fp,"%li",  &RowIndx[k] );
            fscanf(fp,"%li",  &ColIndx[k] );
            fscanf(fp,"%lf", &a[k]);
   }


  fclose(fp);

}

void TripletSparseMatrix::Save(const char* fname) const
{

  int k;

  FILE *fp;

  if ( (fp = fopen(fname,"w")) == NULL )

  {  sp_error_message( "Error opening file in TripletSparseMatrix::Save()"); }


//   fprintf(fp,"%li\t%li\t%li\n", n, m, nz);

   for (k=0;k<nz;k++) {
            fprintf(fp,"%li\t%li\t%e\n",  RowIndx[k]+1, ColIndx[k]+1, a[k] );
   }

   fclose(fp);

}



TripletSparseMatrix tra(const TripletSparseMatrix& A)
{
// Returns a Sparse matrix object with the transposition of input sparse matrix A.

   TripletSparseMatrix sp(A.cols(), A.rows(), A.GetNonZero() );

   memcpy(sp.RowIndx, A.ColIndx, A.nz*sizeof(int) );
   memcpy(sp.ColIndx, A.RowIndx, A.nz*sizeof(int) );
   memcpy(sp.a      , A.a,       A.nz*sizeof(double) );

   return (sp);

}





MatrixXd TripletSparseMatrix::col(int j) const
{
// Returns a dense column vector (MatrixXd object) with the j-th column of the calling
// sparse matrix.
  MatrixXd Temp;

  Temp.resize( this->rows(), 1 );
  Temp.setZero();
  for(int i=0;i<nz;i++) {
       if (ColIndx[i]==j) {
       		Temp( RowIndx[i], 0 ) = a[i];
       }
  }
  return Temp;

}

MatrixXd TripletSparseMatrix::row(int i) const
{
// Returns a dense row vector (MatrixXd object) with the i-th row of the calling
// sparse matrix.
  MatrixXd Temp;

  Temp.resize( 1, this->cols());
  Temp.setZero();
  for(int j=0;j<nz;j++) {
       if (RowIndx[j]==i) {
       		Temp(0, ColIndx[j] ) = a[i];
       }
  }
  return Temp;
}




void TripletSparseMatrix::Transpose()
{
// Transposes a sparse matrix object
   int s;
   int *temp = new int[nz];


   memcpy(temp, RowIndx,    nz*sizeof(int) );
   memcpy(RowIndx, ColIndx, nz*sizeof(int) );
   memcpy(ColIndx, temp,    nz*sizeof(int) );

   s = n;
   n = m;
   m = s;

   return;

}


void sp_error_message(const char *error_text)
{
     error_message( error_text );
}
