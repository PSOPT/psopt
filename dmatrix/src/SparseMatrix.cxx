//
//  SparseMatrix.cpp
//  dmatrix
//
//  Created by Philipp Waxweiler on 21.07.20.
//

#include "SparseMatrix.h"
#include <commonlib.h>
#include "helper.h"

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


void SparseMatrix::InsertNonZero(int i, int j, double val)
{
    int k;

    double* anew;
    int* RowNew;
    int* ColNew;
    int nznew=nz+1;
    int eflag = 0;

    anew =  new double[nznew];
    RowNew= new int[nznew];
    ColNew= new int[nznew];

    for (k=0; k< nz; k++)
    {
        if (RowIndx[k]==i && ColIndx[k]==j) {
               a[k]= val;
               eflag = 1;
               break;
        }
    }

    if (!eflag) {
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

    // Now resize the matrix if necessary
    for (k=0;k<nz;k++) {
        if (RowIndx[k]>n)
           n = RowIndx[k];
        if (ColIndx[k]>m)
           m = ColIndx[k];
    }

}


// Operators

SparseMatrix& SparseMatrix::operator += (const SparseMatrix &rval)
{
    int i;

    if ( (this->n != rval.n) || (this->m != rval.m) )
        fprintf(stderr,"\nSparseMatrix::operator += error: matrix dimensions do not agree");

    for (i=0; i<rval.nz; i++)
    {
        if ( (*this)(rval.RowIndx[i],rval.ColIndx[i])!=0.0 )
               this->a[i] += rval.a[i];
        else
               this->InsertNonZero( rval.RowIndx[i], rval.ColIndx[i], rval.a[i] );
    }
    this->Compress();
    return (*this);
}


SparseMatrix& SparseMatrix::operator+ (const SparseMatrix& Other_matrix) const
{
   SparseMatrix* Result;

   Result = new SparseMatrix;

   (*Result)=(*this);

   return ((*Result)+=Other_matrix);

}

SparseMatrix& SparseMatrix::operator -= (const SparseMatrix &rval)
{
    int i;

    if ( (this->n != rval.n) || (this->m != rval.m) )
        fprintf(stderr,"\nSparseMatrix::operator -= error: matrix dimensions do not agree");

    for (i=0; i<rval.nz; i++)
    {
        if ( (*this)(rval.RowIndx[i],rval.ColIndx[i])!=0.0 )
               this->a[i] -= rval.a[i];
        else
               this->InsertNonZero( rval.RowIndx[i], rval.ColIndx[i], -rval.a[i] );
    }
    this->Compress();
    return (*this);
}


SparseMatrix& SparseMatrix::operator- (const SparseMatrix& Other_matrix) const
{
   SparseMatrix* Result;

   Result = new SparseMatrix;

   (*Result)=(*this);

   return ((*Result)-=Other_matrix);
}


SparseMatrix& SparseMatrix::operator* (DMatrix& A) const
{
  int k;
  int i, j, icol, q;
  DMatrix x;
  DMatrix r(n,1);
  SparseMatrix* sp;

  sp = new SparseMatrix( n, A.GetNoCols(), 0 );


  for (icol=1; icol<=A.m; icol++) {
     x = A(colon(),(long) icol);
     r.FillWithZeros();
     for (k = 0; k < nz; k++) {
         j    = ColIndx[k];
         i    = RowIndx[k];
         r(i) += a[k]*x(j);
     }
     for (q=1; q<=length(r); q++) {
        sp->InsertNonZero(q,icol, r(q) );
     }
  }

  sp->Compress();

  return (*sp);

}

SparseMatrix& SparseMatrix::operator*= (double Arg)
{
   int i;
   for (i=0; i<nz; i++)
   {
      a[i]*= Arg;
   }
   return (*this);
}


SparseMatrix& SparseMatrix::operator* (double Arg) const
{
   SparseMatrix* Result;

   Result = new SparseMatrix;

   (*Result)=(*this);

   return ( (*Result)*=Arg );
}

SparseMatrix& operator *(double Arg, const SparseMatrix& A)
{
   return (A*Arg);
}


cs* convert2triplet(const cs* C)
{
    int m, n, nz, k, j, *Cp, *Ci, *w, *Ti, *Tj ;
    int jlim;
    double *Cx, *Tx ;
    cs *T ;
    if (C->nz>=0) return (NULL) ;        /* check inputs */

    m = C->m ; n = C->n ;
    nz = C->p[n];
    T = cs_spalloc (m, n, nz, 1, 1) ;        /* allocate result */
    Ti = T->i ; Tj = T->p ; Tx = T->x ;
    T->nz = nz;
    w = (int*) cs_calloc (n, sizeof (int)) ;            /* get workspace */
    if (!T || !w) return (cs_done (T, w, NULL, 0)) ;    /* out of memory */
    Cp = C->p ; Ci = C->i ; Cx = C->x ;
    for (k=0; k<C->n; k++)
    {
       if (k==C->n) jlim = nz;
       else jlim = Cp[k+1];
       for (j=Cp[k]; j<jlim; j++)
       if(j>=0)
       {
          Tj[ j ]  =  k;
          Tx[ j ]  =  Cx[ j ];
          Ti[ j ]  =  Ci[ j ];
       }
    }

    return (cs_done (T, w, NULL, 1)) ;        /* success; free w and return T */


}



SparseMatrix& Product(const SparseMatrix& A,const SparseMatrix& B)
{

    /* cs *cs_multiply (const cs *A, const cs *B) */
    int i;
    SparseMatrix *r = new SparseMatrix(A.GetNoRows(), B.GetNoCols(),0);

    int an = A.GetNoRows();
    int am = A.GetNoCols();
    int bn = B.GetNoRows();
    int bm = B.GetNoCols();
    int *ira  = A.GetRowIndx_C_Array();
    int *irb  = B.GetRowIndx_C_Array();
    int *ica  = A.GetColIndx_C_Array();
    int *icb  = B.GetColIndx_C_Array();
    int anz = A.GetNonZero();
    int bnz = B.GetNonZero();
    double *a = A.GetPr();
    double *b = B.GetPr();
    int *irr;
    int *icr;




    cs* Acs = cs_spalloc( an, am, anz, 1, 1 );
    cs* Bcs = cs_spalloc( bn, bm, bnz, 1, 1 );
    cs* CAcs;
    cs* CBcs;
    cs* TABcs;

    Acs->nz=anz;
    Bcs->nz=bnz;

    for(i=0;i<anz;i++)
    {
          Acs->i[i]=ira[i]-1;
          Acs->p[i]=ica[i]-1;
          Acs->x[i]=a[i];
    }

    for(i=0;i<bnz;i++)
    {
          Bcs->i[i]=irb[i]-1;
          Bcs->p[i]=icb[i]-1;
          Bcs->x[i]=b[i];
    }

    cs*  ABcs;

    CAcs = cs_compress (Acs);
    CBcs = cs_compress (Bcs);
    ABcs = cs_multiply( CAcs, CBcs );

    r->Resize(an, bm, ABcs->nzmax);

    TABcs = convert2triplet(ABcs);

    irr = r->GetRowIndx_C_Array();
    icr = r->GetColIndx_C_Array();

    for(i=0;i< TABcs->nz; i++) {
        irr[i]= TABcs->i[i]+1;
        icr[i]= TABcs->p[i]+1;
    }


    memcpy( r->GetPr()      , ABcs->x, TABcs->nzmax*sizeof(double) );

    cs_spfree (ABcs);
    cs_spfree (TABcs);
/*    cs_spfree (Acs);*/
    cs_spfree (Bcs);
    cs_spfree (CAcs);
    cs_spfree (CBcs);


    return (*r);

}


// {
//           int i, j;
//
//           double sum;
//
//           SparseMatrix *sp = new SparseMatrix( A.GetNoRows(), B.GetNoCols(), 0 );
//
//           VecElem *row_p = new VecElem;
//           VecElem *col_p = new VecElem;
//
//           if (A.GetNoCols() != B.GetNoRows())
//                 sp_error_message("Inconsistent matrix dimensions in SparseMatrix::operator*");
//
//
//
//           for (i = 1; i <= A.GetNoRows(); i++) {
//                     for (j = 1; j <= B.GetNoCols(); j++) {
//                        //     row_p = A->rows [i];
//                        //     col_p = B->cols [j];
//                            row_p->index = 0;
//                            col_p->index = 0;
//                            row_p->previous = 0;
//                            col_p->previous = 0;
//                            sum = 0;
//
//                            A.    Find_Next_NonZero_in_Row(i, row_p);
//                            B.    Find_Next_NonZero_in_Col(j, col_p);
//                            while ((row_p->index != -1) && (col_p->index != -1)) {
//                                    if (row_p->index == col_p->index) {
//                                            sum += row_p->value * col_p->value;
//                                            A.    Find_Next_NonZero_in_Row(i, row_p);
//                                            B.    Find_Next_NonZero_in_Col(j, col_p);
//                                    }
//                                    else if (row_p->index < col_p->index) {
//                                            A.    Find_Next_NonZero_in_Row(i, row_p);
//                                    }
//                                    else {
//                                         B.Find_Next_NonZero_in_Col(j, col_p);
//                                    }
//                            }
//                            if (sum!=0.0) sp->InsertNonZero(i, j, sum);
//                    }
//          }
//
//          return (*sp);
//
// }



SparseMatrix& SparseMatrix::operator* (const SparseMatrix& B) const
{

         return Product(*this,B);

}


SparseMatrix& TProduct(SparseMatrix& A,SparseMatrix& B)
{
      SparseMatrix *sp= new SparseMatrix(A.GetNoRows(),B.GetNoCols(),0);

          A.Transpose();

          {*sp=Product(A,B);}

          A.Transpose();

          return (*sp);

}


SparseMatrix& ProductT(SparseMatrix& A,SparseMatrix& B)
{

      SparseMatrix *sp= new SparseMatrix(A.GetNoCols(),B.GetNoCols(),0);

          B.Transpose();

          {*sp=Product(A,B);}

          B.Transpose();

          return (*sp);


}

SparseMatrix& elemProduct(const SparseMatrix A, const SparseMatrix& B)
{
          int i, j, k;

          double sum;

          double *a = A.GetPr();

          int *RowIndx = A.GetRowIndx_C_Array();
          int *ColIndx = A.GetColIndx_C_Array();

          SparseMatrix *sp = new SparseMatrix( A.GetNoRows(), B.GetNoCols(), 0 );

          if (A.GetNoRows() != B.GetNoRows() ||  A.GetNoCols() != B.GetNoCols())
                sp_error_message("Inconsistent matrix dimensions in SparseMatrix elemProduct()");



          for (k = 0; k < A.GetNonZero(); k++) {
                i = RowIndx[k];
                j = ColIndx[k];
                sum = a[k]*B(i,j);
                if (sum!=0.0) sp->InsertNonZero(i, j, sum);
          }

          return (*sp);
}


SparseMatrix& SparseMatrix::operator/= (double Arg)
{
   int i;

   if (Arg==0.0) sp_error_message("Division by zero in SparseMatrix::operator/");

   for (i=0; i<nz; i++)
   {
      a[i]/= Arg;
   }
   return (*this);
}

SparseMatrix& SparseMatrix::operator/ (double Arg) const
{
   SparseMatrix* Result;

   if (Arg==0.0) sp_error_message("Division by zero in SparseMatrix::operator/");

   Result = new SparseMatrix;

   (*Result)=(*this);

   return ( (*Result)/=Arg );
}

SparseMatrix& SparseMatrix::operator= (const SparseMatrix& Other_matrix)
{
     this->Resize(Other_matrix.n, Other_matrix.m, Other_matrix.nz);
     memcpy(a, Other_matrix.a, nz*sizeof(double) );
     memcpy(RowIndx, Other_matrix.RowIndx, nz*sizeof(int) );
     memcpy(ColIndx, Other_matrix.ColIndx, nz*sizeof(int) );
     return (*this);
}

double& SparseMatrix::operator() (int i,int j)
{
   int k;

    int ii;

    int eflag = 0;

    for (k=0; k< nz; k++)
    {
       if (RowIndx[k]==i && ColIndx[k]==j) {
              ii = k;
              eflag = 1;
              break;
       }

    }

    if( !eflag ) {

        for (k=0; k<nz; k++) {

              if (RowIndx[k]==0 || ColIndx[k]==0) {
                   RowIndx[k]=i;
                   ColIndx[k]=j;
                        ii=k;
                        eflag=1;
                        break;
                 }
        }

    }

    if (!eflag) {
       this->InsertNonZero(i,j,0.0);
       ii = nz-1;
    }

    return a[ii];

}

double SparseMatrix::operator() (int row,int col) const
{
      int i;
      double retval = 0.0;

      if ( (row>n) || (row<0) || (col<0) || (col>m) ) {
          error_message("Out of range index in SparseMatrix::operator()");
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

//double SparseMatrix::operator() (int row,int col)
//{
//
//}

void SparseMatrix::Print(const char* text)
{
   int i;

   fprintf(stderr,"\nSparse matrix %s",text);
   fprintf(stderr,"\nNumber of rows: %i", n);
   fprintf(stderr,"\nNumber of columns: %i", m);
   fprintf(stderr,"\nNumber of non-zero elements: %i", nz);
   if (n*m!=0) fprintf(stderr,"\nDensity: %f%%", (  ((double) nz*100)/(double) (n*m)  ) );
   if (nz>0) fprintf(stderr,"\n(Row,Col)\tValue");

   for (i=0; i<nz; i++)
   {
       fprintf(stderr,"\n(%i,%i)\t\t%e", RowIndx[i], ColIndx[i], a[i]);

   }
   fprintf(stderr,"\n");

}

void SparseMatrix::SaveSparsityPattern(const char* text) const
{
   int i, j;
   FILE *outfile;

   outfile = fopen(text,"w");

   for (i=1; i<=n; i++)  {
       for(j=1; j<=m; j++)  {
           if ( (*this)(i,j)!=0.0 )
               fprintf(outfile,"%c ",'*');
           else
               fprintf(outfile,"%c ",'0');
       }
       fprintf(outfile,"\n");
   }
   fclose(outfile);
}


void SparseMatrix::Compress()
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
        if (fabs(a[i])<=(DMatrix::GetEPS()) ) {
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

SparseMatrix& sparse( const DMatrix& A )
{
// Returns a SparseMatrix object created from a DMatrix object
     SparseMatrix *sp = new SparseMatrix(A);
     return (*sp);
}

SparseMatrix& speye(int n)
{
// Returns a sparse identity matrix
     SparseMatrix *sp = new SparseMatrix;
     int i;
     sp->Resize(n, n, n);

     for (i=0;i<sp->nz;i++) {
         sp->RowIndx[i]=i+1;
         sp->ColIndx[i]=i+1;
         sp->a[i]=1.0;
     }

     return (*sp);
}

SparseMatrix& spones(const SparseMatrix& s)
{
// Returns a sparse identity matrix
     SparseMatrix *sp = new SparseMatrix(s);
     int i;

     for (i=0;i<sp->nz;i++) {
         sp->a[i]=1.0;
     }

     return (*sp);
}

SparseMatrix& spconvert(DMatrix& A)
{
   int i;
   if (A.GetNoCols()!=3)
      sp_error_message("Incorrect number of columns of input matrix in spconvert()");
   DMatrix I;
   DMatrix J;
   DMatrix Aij;
   int nzmax = A.GetNoRows();

   I = A(colon(),1);
   J = A(colon(),2);
   Aij = A(colon(),3);

   int maxI = (int) Max(I);
   int maxJ = (int) Max(J);

   SparseMatrix* sp = new SparseMatrix(maxI, maxJ, nzmax);

   for (i=0; i<nzmax; i++)
   {
       sp->RowIndx[i] = (int) I(i+1);
       sp->ColIndx[i] = (int) J(i+1);
       sp->a[i]       = Aij(i+1);
   }

   sp->Compress();

   return (*sp);
}


void SparseMatrix::Load(const char* fname)
{
  int nrow, ncol, nnz;

  int k;

  FILE *fp;

  if ( (fp = fopen(fname,"r")) == NULL )

  {  sp_error_message( "Error opening file in SparseMatrix::Load()"); }


   fscanf(fp,"%i", &nrow);
   fscanf(fp,"%i", &ncol);
   fscanf(fp,"%i", &nnz);

   this->Resize(nrow, ncol, nnz);

   for (k=0;k<nnz;k++) {
            fscanf(fp,"%i",  &RowIndx[k] );
            fscanf(fp,"%i",  &ColIndx[k] );
            fscanf(fp,"%lf", &a[k]);
   }


  fclose(fp);

}

void SparseMatrix::Save(const char* fname) const
{

  int k;

  FILE *fp;

  if ( (fp = fopen(fname,"w")) == NULL )

  {  sp_error_message( "Error opening file in SparseMatrix::Save()"); }


   fprintf(fp,"%i\t%i\t%i\n", n, m, nz);

   for (k=0;k<nz;k++) {
            fprintf(fp,"%i\t%i\t%e\n",  RowIndx[k], ColIndx[k], a[k] );
   }

   fclose(fp);

}



SparseMatrix& tra(const SparseMatrix& A)
{
// Returns a Sparse matrix object with the transposition of input sparse matrix A.

   SparseMatrix *sp = new SparseMatrix(A.GetNoCols(), A.GetNoRows(), A.GetNonZero() );

   memcpy(sp->RowIndx, A.ColIndx, A.nz*sizeof(int) );
   memcpy(sp->ColIndx, A.RowIndx, A.nz*sizeof(int) );
   memcpy(sp->a      , A.a,       A.nz*sizeof(double) );

   return (*sp);

}

SparseMatrix& inv(SparseMatrix& A)
{
// Returns the inverse (if it exists) of a sparse matrix A.
   int i,j;

   if (A.GetNoCols()!=A.GetNoRows())
      sp_error_message("Sparse matrix must be square in function inv()");

   void* LUSOLv;

   DMatrix x(A.GetNoRows(),1);

   DMatrix r(A.GetNoRows(),1);

   DMatrix b(A.GetNoRows(),1);

   SparseMatrix *sp = new SparseMatrix(A.GetNoRows(), A.GetNoCols(), 0);

   LUSOLv = SparseLUFactor(A);

   for(j=1; j<=A.GetNoCols(); j++)
   {
       b.FillWithZeros();
       b(j) = 1.0;
       SparseLUSolveGivenFactor(LUSOLv, A, b, x, r);
       for(i=1;i<=A.GetNoRows();i++)
       {
            if (x(i)!=0.0)
                sp->InsertNonZero(i,j,x(i));
       }
   }

   return (*sp);

}

SparseMatrix& SparseMatrix::operator %(SparseMatrix& B)
{
// Left divition implementation A%B is conceptually equivalent to inv(A)*B
// but it is numerically more efficient
   int i,j,k;

   if (this->GetNoCols()!=this->GetNoRows())
      sp_error_message("Sparse matrix must be square in function SparseMatrix::operator%");

   if (this->GetNoCols()!=B.GetNoRows())
      sp_error_message("Incorrect dimensions in SparseMatrix::operator% ");

   void* LUSOLv;

   DMatrix x(this->GetNoCols(),1);

   DMatrix r(B.GetNoRows(),1);

   DMatrix b(B.GetNoRows(),1);

   SparseMatrix *sp = new SparseMatrix(this->GetNoRows(), B.GetNoCols(), 0);

   LUSOLv = SparseLUFactor(*this);

   for(j=1; j<=B.GetNoCols(); j++)
   {
       for (k=1;k<=B.GetNoRows();k++) b(k) = B(k,j);
       SparseLUSolveGivenFactor(LUSOLv, *this, b, x, r);
       for(i=1;i<=this->GetNoRows();i++)
       {
            if (x(i)!=0.0)
                sp->InsertNonZero(i,j,x(i));
       }
   }

   return (*sp);

}


SparseMatrix& SparseMatrix::operator/ (SparseMatrix& B)
{
// Right divition implementation A/B is conceptually equivalent to A*inv(B)
// but it is numerically more efficient

   if (this->GetNoCols()!=this->GetNoRows())
      sp_error_message("Sparse matrix must be square in function SparseMatrix::operator/");

   if (this->GetNoCols()!=B.GetNoRows())
      sp_error_message("Incorrect dimensions in SparseMatrix::operator/ ");

   SparseMatrix *sp = new SparseMatrix(this->GetNoRows(), B.GetNoCols(), 0);

   B.Transpose();

   this->Transpose();

   {*sp = B%(*this);}

   sp->Transpose();

   B.Transpose();

   this->Transpose();

   return (*sp);

}

SparseMatrix& SparseMatrix::operator||(SparseMatrix& B) const
{
  // Returns the horizontal concatenation of two sparse matrices

   int i;

   if (this->GetNoRows()!=B.GetNoRows())
      sp_error_message("Number of rows must be equal in SparseMatrix::operator||");

   SparseMatrix *sp = new SparseMatrix(n, m+B.m, nz+B.nz);

   memcpy(sp->a, a, nz*sizeof(double) );
   memcpy(sp->RowIndx, RowIndx, nz*sizeof(int) );
   memcpy(sp->ColIndx, ColIndx, nz*sizeof(int) );

   memcpy(sp->a+nz, B.a, B.nz*sizeof(double) );
   memcpy(sp->RowIndx+nz, B.RowIndx, B.nz*sizeof(int) );

   for(i=0; i<B.nz; i++)
   {
       sp->ColIndx[nz+i] = m + B.ColIndx[i];
   }

   return (*sp);

}



SparseMatrix& SparseMatrix::operator&&(SparseMatrix& B) const
{
   // Returns the vertical concatenation of two sparse matrices

   int i;

   if (m!=B.m)
      sp_error_message("Number of cols must be equal in SparseMatrix::operator&&");

   SparseMatrix *sp = new SparseMatrix(n + B.n, m, nz+B.nz);

   memcpy(sp->a, a, nz*sizeof(double) );
   memcpy(sp->RowIndx, RowIndx, nz*sizeof(int) );
   memcpy(sp->ColIndx, ColIndx, nz*sizeof(int) );

   memcpy(sp->a+nz, B.a, B.nz*sizeof(double) );
   memcpy(sp->ColIndx+nz, B.ColIndx, B.nz*sizeof(int) );

   for(i=0; i<B.nz; i++)
   {
       sp->RowIndx[nz+i] = n + B.RowIndx[i];
   }

   return (*sp);

}

SparseMatrix& SparseMatrix::sub_matrix(int istart, int iend, int jstart, int jend) const
{
  // Returns a submatrix of the calling object
   int i,j;
   double value;

   if ((istart<1)||(istart>n)||(iend<1)||(iend>n)||(jstart<1)||(jstart>m)||(istart>iend)||(jstart>jend))
      sp_error_message("Index error in SparseMatrix::sub_matrix");

   int nnew = iend-istart+1;
   int mnew = jend-jstart+1;

   SparseMatrix *sp = new SparseMatrix(nnew, mnew, 0);

   for (i=1; i<=nnew; i++) {
          for(j=1; j<=mnew; j++) {
                 value = (*this)(istart+i-1,jstart+j-1);
                 if ( value!=0.0  ) {
              sp->InsertNonZero(i,j,value);
                 }
      }
   }

   return (*sp);



}


DMatrix& SparseMatrix::Column(int j) const
{
// Returns a dense column vector (DMatrix object) with the j-th column of the calling
// sparse matrix.
  DMatrix* Temp;
  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
  Temp->Resize( this->GetNoRows(), 1 );
  Temp->FillWithZeros();
  for(int i=0;i<nz;i++) {
       if (ColIndx[i]==j) {
               (*Temp)( RowIndx[i], 1 ) = a[i];
       }
  }
  return (*Temp);

}

DMatrix& SparseMatrix::Row(int i) const
{
// Returns a dense row vector (DMatrix object) with the i-th row of the calling
// sparse matrix.
  DMatrix* Temp;
  Temp = DMatrix::GetTempPr(ChkTmpIndx(DMatrix::IncrementAuxIndx()));
  Temp->Resize( 1, this->GetNoCols());
  Temp->FillWithZeros();
  for(int j=0;j<nz;j++) {
       if (RowIndx[j]==i) {
               (*Temp)(1, ColIndx[j] ) = a[i];
       }
  }
  return (*Temp);
}

void SparseMatrix::set_sub_matrix(const SparseMatrix& B, int istart, int iend, int jstart, int jend)
{
   int i, j;
   double value;

   if ((istart<1)||(istart>n)||(iend<1)||(iend>n)||(jstart<1)||(jstart>m)||(istart>iend)||(jstart>jend))
      sp_error_message("Index error in SparseMatrix::set_sub_matrix");

   int nn = iend-istart+1;
   int mm = jend-jstart+1;

   if (B.n != nn || B.m != mm)
      sp_error_message("Range must be consistent with input matrix size in SparseMatrix::set_sub_matrix");

   for (i=1; i<=nn; i++) {
          for(j=1; j<=mm; j++) {
                 value = B(i,j);
                 if ( value!=0.0  ) {
              this->InsertNonZero(istart+i-1,jstart+j-1,value);
                 }
      }
   }

   return;

}



void SparseMatrix::Transpose()
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

SparseMatrix& SparseMatrix::operator^(double x) const
{
  // Element-wise power of a sparse matrix
   int i;
   SparseMatrix *sp = new SparseMatrix(this->GetNoRows(), this->GetNoCols(), this->GetNonZero() );

   memcpy(sp->RowIndx, this->RowIndx, this->nz*sizeof(int) );
   memcpy(sp->ColIndx, this->ColIndx, this->nz*sizeof(int) );
   for (i=0; i<this->nz; i++)
   {
         sp->a[i] = pow( this->a[i], x );
   }

   return (*sp);


}




MYBOOL isNum(char val)
{
  int ord;
  ord = (int) val - 48;
  return( (MYBOOL) ((ord >= 0) && (ord <= 9)) );
}
void SparseLUSolve(SparseMatrix& A, DMatrix& bb, DMatrix& xx, DMatrix& rr,  int argc, char *argv[])
{

/* Output device */
  FILE *outunit = stderr;

  int diagnostics_level=1;

/* Overall dimensions allocated */
  int    maxm = MAXROWS, maxn = MAXCOLS, maxnz = MAXNZ,
         replace = 0, randcol = 0;
  MYBOOL ftran = TRUE;

/* Storage for A, b */
  REAL   *Aij, *b, *xexact;
  int    *iA, *jA;

/* Storage for LUSOL */
  LUSOLrec *LUSOL = NULL;

/* Define local storage variables */
  int  i     , inform, j     , k     , i1   ,
       m     ,  lenb, lenx,
       n     , nelem , nnzero;
  REAL Amax  , test  ,
       bnorm , rnorm , xnorm,
       *rhs  , *r    , *x;
  MYBOOL printsolution = FALSE, success = TRUE;

/* Create the LUSOL object and set user options */
  LUSOL = LUSOL_create(outunit, 0, LUSOL_PIVMOD_TPP, 0);
  LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = 10;
  i = 1;
  n = A.GetNoCols();
  m = A.GetNoRows();
  nnzero = A.GetNonZero();


  while((n == 0) && (i < argc)) {
    if(strcmp("-p", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < 0 || m > LUSOL_PIVMOD_MAX)
          continue;
        LUSOL->luparm[LUSOL_IP_PIVOTTYPE] = m;
      }
    }
    else if(strcmp("-t", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        Amax = atof(argv[i]);
        if(Amax < 1 || Amax > 100)
          continue;
        LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = Amax;
      }
    }
    else if(strcmp("-m", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < LUSOL_MULT_nz_a || m > 100)
          continue;
        LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = atoi(argv[i]);
      }
    }
    else if(strcmp("-s", argv[i]) == 0)
      printsolution = TRUE;
    else if(strcmp("-b", argv[i]) == 0)
      ftran = FALSE;
    else if(strcmp("-r", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < 0 || m > 10)
          continue;
      }
      else
        m = 1;
      srand((unsigned) time( NULL ));
      replace = 2*m;
    }
  } /* end while */



  maxn = A.GetNoCols();
  maxm = A.GetNoRows();
  maxnz= A.GetNonZero();

/* Create the arrays */

  Aij = A.GetPr() - BLAS_BASE;
  iA =  A.GetRowIndx_C_Array() - BLAS_BASE;
  jA =  A.GetColIndx_C_Array() - BLAS_BASE;
  if(ftran)
    lenb = maxm;
  else
    lenb = maxn;
  b   = bb.GetPr()-BLAS_BASE;
  rhs = (REAL *) calloc(lenb+BLAS_BASE, sizeof(REAL));

  if(ftran)
    lenx = maxn;
  else
    lenx = maxm;
  xexact = (REAL *) calloc(lenx+BLAS_BASE, sizeof(REAL));
  r = rr.GetPr() - BLAS_BASE;
  x = xx.GetPr() - BLAS_BASE;


/* -----------------------------------------------------------------
   Load A into (a, indc, indr).
   ----------------------------------------------------------------- */

  if(!LUSOL_assign(LUSOL, iA, jA, Aij, nnzero, TRUE)) {
    fprintf(outunit, "Error: LUSOL failed due to insufficient memory.\n");
    goto x900;
  }

/* ------------------------------------------------------------------
   Factor  A = L U.
   ------------------------------------------------------------------ */
  nelem = nnzero;
  LU1FAC( LUSOL, &inform );
  if (inform > LUSOL_INFORM_SERIOUS) {
    fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    goto x900;
  }
  /* Get the largest element in A; we use it below as an estimate
     of ||A||_inf, even though it isn't a proper norm. */
  Amax = LUSOL->parmlu[LUSOL_RP_MAXELEM_A];

/* ------------------------------------------------------------------
   SOLVE  A x = b.
   Save b first because lu6sol() overwrites rhs when computing x.
   ------------------------------------------------------------------ */
Resolve:
#if 1
  MEMCOPY(x, b, lenb+BLAS_BASE);
  if(ftran)
    inform = LUSOL_ftran(LUSOL, x, NULL, FALSE);
  else
    inform = LUSOL_btran(LUSOL, x, NULL);
#else
  MEMCOPY(rhs, b, lenb+BLAS_BASE);
  if(ftran)
    LU6SOL( LUSOL, LUSOL_SOLVE_Aw_v, rhs, x, NULL, &inform );
  else
    LU6SOL( LUSOL, LUSOL_SOLVE_Atv_w, x, rhs, NULL, &inform );
#endif
  if (inform > LUSOL_INFORM_SERIOUS) {
    fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    goto x900;
  }
  if(printsolution) {
    char* msg;
    sprintf(msg,"\nSolution vector=");
    blockWriteREAL(outunit, msg, x, 1, lenb);
  }

/* ------------------------------------------------------------------
   Set r = b - Ax.
   Find norm of r and x.
   ------------------------------------------------------------------ */
  MEMCOPY(r, b, lenb+BLAS_BASE);
  for (k = 1; k <= nnzero; k++) {
    i    = iA[k];
    j    = jA[k];
    if(ftran)
      r[i] -= Aij[k]*x[j];
    else
      r[j] -= Aij[k]*x[i];
  }
  bnorm  = dnormi( lenb, b );
  rnorm  = dnormi( lenb, r );
  xnorm  = dnormi( lenx, x );

/* ------------------------------------------------------------------
   Report the findings.
   ------------------------------------------------------------------ */
  if(randcol > 0)
    fprintf(outunit, "\n\nColumn %d was %s\n",
                      randcol,
                      (mod(replace,2) == 1 ? "replaced with random data" : "restored"));
  if (diagnostics_level>=2) {
      fprintf(outunit, "\nLU size statistics (%d reallocations):\n",
                   LUSOL->expanded_a);
  }
  test = LUSOL->luparm[LUSOL_IP_NONZEROS_U0]+LUSOL->luparm[LUSOL_IP_NONZEROS_ROW];
  if (diagnostics_level>=2) {
      fprintf(outunit, "L0-size = %d   U0-size = %d   LU-nonzeros = %d   Fill-in = %.1fx\n",
                   LUSOL->luparm[LUSOL_IP_NONZEROS_L0],
                   LUSOL->luparm[LUSOL_IP_NONZEROS_U0],
                   (int) test, test/nnzero);
  }
  test   = rnorm / (Amax*xnorm);
  if (diagnostics_level>=2) {
      fprintf(outunit, "\nAccuracy statistics:\n");
      fprintf(outunit, "%s with a factor tolerance of %g gave a relative error of %g\n",
                   LUSOL_pivotLabel(LUSOL), LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij], test);
      fprintf(outunit, "Amax = %g   bnorm = %g   rnorm = %g   xnorm = %g\n",
                   Amax, bnorm, rnorm, xnorm);
      fprintf(outunit, "\n");
  }
  if (test <= 1.0e-8 ) {}
      /* fprintf(outunit, "The equations were solved with very high accuracy.\n");*/
  else if (test <= 1.0e-6) {}
  //    fprintf(outunit, "The equations were solved with reasonably good accuracy.\n");
  else {
    if (test <= 1.0e-4)
      fprintf(outunit, "SparseLUSolve: Questionable accuracy; the LU factors may not be good enough.\n");
    else
      fprintf(outunit, "SparseLUSolve: Poor accuracy; the LU factorization probably failed.\n");
    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVMOD_TRP)
      fprintf(outunit, "SparseLUSolve: Try a smaller factor tolerance (current is %g).\n",
                       LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]);
    else
      fprintf(outunit, "SparseLUSolve: Try a smaller factor tolerance and/or TRP pivoting.\n");
  }

 /* Check if we should replace a column and resolve */
  if(replace > 0) {
    replace--;
    if(mod(replace, 2) == 1) {
      /* Randomly find a column and replace the data with the data in b */
      rnorm   = rand();
      randcol = (int) (n * rnorm / (RAND_MAX+1)) + 1;
#if 1
      MEMCLEAR(x, m+1);
      for(i = 1; i < m; i++)
        x[i] = Amax * rand() / RAND_MAX;
      inform = LUSOL_replaceColumn(LUSOL, randcol, x);
#else
      inform = LUSOL_replaceColumn(LUSOL, randcol, b);
#endif
    }
    else {
      /* Set the previously replaced column back and resolve */
      for (k = 1; k <= nnzero; k++) {
        i    = iA[k];
        j    = jA[k];
        if(j == randcol)
          x[i] = Aij[k];
      }
      inform = LUSOL_replaceColumn(LUSOL, randcol, x);
    }
    if(inform != LUSOL_INFORM_LUSUCCESS)
      fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    else
      goto Resolve;
  }


/* Free memory */
x900:
  if(!success)
    fprintf(outunit, "Insufficient memory or data file not found.\n");

  free(xexact);

  LUSOL_free(LUSOL);


}


void* SparseLUFactor(SparseMatrix& A, int argc, char *argv[])
{
/* Output device */
  FILE *outunit = stderr;

/* Overall dimensions allocated */
  int    maxm = MAXROWS, maxn = MAXCOLS, maxnz = MAXNZ,
         replace = 0;
  MYBOOL ftran = TRUE;

/* Storage for A, b */
  REAL   *Aij,  *xexact;
  int    *iA, *jA;

/* Storage for LUSOL */
  LUSOLrec *LUSOL = NULL;

/* Define local storage variables */
  int  i     , inform,  i1   ,
       m     , lenx,
       n     , nelem , nnzero;
  REAL Amax;

  MYBOOL printsolution = FALSE;

/* Create the LUSOL object and set user options */
  LUSOL = LUSOL_create(outunit, 0, LUSOL_PIVMOD_TPP, 0);
  LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = 10;
  i = 1;
  n = A.GetNoCols();
  m = A.GetNoRows();
  nnzero = A.GetNonZero();


  while((n == 0) && (i < argc)) {
    if(strcmp("-p", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < 0 || m > LUSOL_PIVMOD_MAX)
          continue;
        LUSOL->luparm[LUSOL_IP_PIVOTTYPE] = m;
      }
    }
    else if(strcmp("-t", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        Amax = atof(argv[i]);
        if(Amax < 1 || Amax > 100)
          continue;
        LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = Amax;
      }
    }
    else if(strcmp("-m", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < LUSOL_MULT_nz_a || m > 100)
          continue;
        LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = atoi(argv[i]);
      }
    }
    else if(strcmp("-s", argv[i]) == 0)
      printsolution = TRUE;
    else if(strcmp("-b", argv[i]) == 0)
      ftran = FALSE;
    else if(strcmp("-r", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < 0 || m > 10)
          continue;
      }
      else
        m = 1;
      srand((unsigned) time( NULL ));
      replace = 2*m;
    }
  } /* end while */



  maxn = A.GetNoCols();
  maxm = A.GetNoRows();
  maxnz= A.GetNonZero();

/* Create the arrays */

  Aij = A.GetPr() - BLAS_BASE;
  iA =  A.GetRowIndx_C_Array() - BLAS_BASE;
  jA =  A.GetColIndx_C_Array() - BLAS_BASE;


  if(ftran)
    lenx = maxn;
  else
    lenx = maxm;
  xexact = (REAL *) calloc(lenx+BLAS_BASE, sizeof(REAL));

/* -----------------------------------------------------------------
   Load A into (a, indc, indr).
   ----------------------------------------------------------------- */

  if(!LUSOL_assign(LUSOL, iA, jA, Aij, nnzero, TRUE)) {
    fprintf(outunit, "Error: LUSOL failed due to insufficient memory.\n");
    goto x900;
  }

/* ------------------------------------------------------------------
   Factor  A = L U.
   ------------------------------------------------------------------ */
  nelem = nnzero;
  LU1FAC( LUSOL, &inform );
  if (inform > LUSOL_INFORM_SERIOUS) {
    fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    goto x900;
  }

x900:
  return ((void*) LUSOL);

}

void SparseLUSolveGivenFactor(void* LUSOLv, SparseMatrix& A, DMatrix& bb, DMatrix& xx, DMatrix& rr,  int argc, char *argv[])
{
/* Output device */
  FILE *outunit = stderr;

  int diagnostics_level=1;

/* Overall dimensions allocated */
  int    maxm = MAXROWS, maxn = MAXCOLS, maxnz = MAXNZ,
         replace = 0, randcol = 0;
  MYBOOL ftran = TRUE;

/* Storage for A, b */
  REAL   *Aij, *b, *xexact;
  int    *iA, *jA;

/* Storage for LUSOL */
  LUSOLrec *LUSOL = (LUSOLrec*) LUSOLv;

/* Define local storage variables */
  int  i     , inform, j     , k     , i1   ,
       m     , lenb, lenx,
       n     , nnzero;
  REAL Amax  , test  ,
       bnorm , rnorm , xnorm,
       *rhs  , *r    , *x;
  MYBOOL printsolution = FALSE, success = TRUE;

/* set user options */

  i = 1;
  n = LUSOL->n;
  m = LUSOL->m;
  nnzero = LUSOL->nelem;


  while((n == 0) && (i < argc)) {
    if(strcmp("-p", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < 0 || m > LUSOL_PIVMOD_MAX)
          continue;
        LUSOL->luparm[LUSOL_IP_PIVOTTYPE] = m;
      }
    }
    else if(strcmp("-t", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        Amax = atof(argv[i]);
        if(Amax < 1 || Amax > 100)
          continue;
        LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij] = Amax;
      }
    }
    else if(strcmp("-m", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < LUSOL_MULT_nz_a || m > 100)
          continue;
        LUSOL->luparm[LUSOL_IP_SCALAR_NZA] = atoi(argv[i]);
      }
    }
    else if(strcmp("-s", argv[i]) == 0)
      printsolution = TRUE;
    else if(strcmp("-b", argv[i]) == 0)
      ftran = FALSE;
    else if(strcmp("-r", argv[i]) == 0) {
      i1 = i+1;
      if((i1 < argc) && isNum(argv[i1][1])) {
        i = i1;
        m = atoi(argv[i]);
        if(m < 0 || m > 10)
          continue;
      }
      else
        m = 1;
      srand((unsigned) time( NULL ));
      replace = 2*m;
    }
  } /* end while */



  maxn = LUSOL->n;
  maxm = LUSOL->m;
  maxnz= LUSOL->nelem;

/* Create the arrays */

  Aij = A.GetPr() - BLAS_BASE;
  iA =  A.GetRowIndx_C_Array() - BLAS_BASE;
  jA =  A.GetColIndx_C_Array() - BLAS_BASE;

  if(ftran)
    lenb = maxm;
  else
    lenb = maxn;
  b   = bb.GetPr()-BLAS_BASE;
  rhs = (REAL *) calloc(lenb+BLAS_BASE, sizeof(REAL));

  if(ftran)
    lenx = maxn;
  else
    lenx = maxm;
  xexact = (REAL *) calloc(lenx+BLAS_BASE, sizeof(REAL));
  r = rr.GetPr() - BLAS_BASE;
  x = xx.GetPr() - BLAS_BASE;


  /* Get the largest element in A; we use it below as an estimate
     of ||A||_inf, even though it isn't a proper norm. */
  Amax = LUSOL->parmlu[LUSOL_RP_MAXELEM_A];

/* ------------------------------------------------------------------
   SOLVE  A x = b.
   Save b first because lu6sol() overwrites rhs when computing x.
   ------------------------------------------------------------------ */
Resolve:
#if 1
  MEMCOPY(x, b, lenb+BLAS_BASE);
  if(ftran)
    inform = LUSOL_ftran(LUSOL, x, NULL, FALSE);
  else
    inform = LUSOL_btran(LUSOL, x, NULL);
#else
  MEMCOPY(rhs, b, lenb+BLAS_BASE);
  if(ftran)
    LU6SOL( LUSOL, LUSOL_SOLVE_Aw_v, rhs, x, NULL, &inform );
  else
    LU6SOL( LUSOL, LUSOL_SOLVE_Atv_w, x, rhs, NULL, &inform );
#endif
  if (inform > LUSOL_INFORM_SERIOUS) {
    fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    goto x900;
  }
  if(printsolution) {
    char* msg;
    sprintf(msg,"\nSolution vector=");
    blockWriteREAL(outunit, msg , x, 1, lenb);
  }

/* ------------------------------------------------------------------
   Set r = b - Ax.
   Find norm of r and x.
   ------------------------------------------------------------------ */
  MEMCOPY(r, b, lenb+BLAS_BASE);
  for (k = 1; k <= nnzero; k++) {
    i    = iA[k];
    j    = jA[k];
    if(ftran)
      r[i] -= Aij[k]*x[j];
    else
      r[j] -= Aij[k]*x[i];
  }
  bnorm  = dnormi( lenb, b );
  rnorm  = dnormi( lenb, r );
  xnorm  = dnormi( lenx, x );

/* ------------------------------------------------------------------
   Report the findings.
   ------------------------------------------------------------------ */
  if(randcol > 0)
    fprintf(outunit, "\n\nColumn %d was %s\n",
                      randcol,
                      (mod(replace,2) == 1 ? "replaced with random data" : "restored"));
  if (diagnostics_level>=2) {
      fprintf(outunit, "\nLU size statistics (%d reallocations):\n",
                   LUSOL->expanded_a);
  }
  test = LUSOL->luparm[LUSOL_IP_NONZEROS_U0]+LUSOL->luparm[LUSOL_IP_NONZEROS_ROW];
  if (diagnostics_level>=2) {
      fprintf(outunit, "L0-size = %d   U0-size = %d   LU-nonzeros = %d   Fill-in = %.1fx\n",
                   LUSOL->luparm[LUSOL_IP_NONZEROS_L0],
                   LUSOL->luparm[LUSOL_IP_NONZEROS_U0],
                   (int) test, test/nnzero);
  }
  test   = rnorm / (Amax*xnorm);
  if (diagnostics_level>=2) {
      fprintf(outunit, "\nAccuracy statistics:\n");
      fprintf(outunit, "%s with a factor tolerance of %g gave a relative error of %g\n",
                   LUSOL_pivotLabel(LUSOL), LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij], test);
      fprintf(outunit, "Amax = %g   bnorm = %g   rnorm = %g   xnorm = %g\n",
                   Amax, bnorm, rnorm, xnorm);
      fprintf(outunit, "\n");
  }
  if (test <= 1.0e-8 ) {}
      /* fprintf(outunit, "The equations were solved with very high accuracy.\n");*/
  else if (test <= 1.0e-6) {}
  //    fprintf(outunit, "The equations were solved with reasonably good accuracy.\n");
  else {
    if (test <= 1.0e-4)
      fprintf(outunit, "SparseLUSolve: Questionable accuracy; the LU factors may not be good enough.\n");
    else
      fprintf(outunit, "SparseLUSolve: Poor accuracy; the LU factorization probably failed.\n");
    if(LUSOL->luparm[LUSOL_IP_PIVOTTYPE] == LUSOL_PIVMOD_TRP)
      fprintf(outunit, "SparseLUSolve: Try a smaller factor tolerance (current is %g).\n",
                       LUSOL->parmlu[LUSOL_RP_FACTORMAX_Lij]);
    else
      fprintf(outunit, "SparseLUSolve: Try a smaller factor tolerance and/or TRP pivoting.\n");
  }

 /* Check if we should replace a column and resolve */
  if(replace > 0) {
    replace--;
    if(mod(replace, 2) == 1) {
      /* Randomly find a column and replace the data with the data in b */
      rnorm   = rand();
      randcol = (int) (n * rnorm / (RAND_MAX+1)) + 1;
#if 1
      MEMCLEAR(x, m+1);
      for(i = 1; i < m; i++)
        x[i] = Amax * rand() / RAND_MAX;
      inform = LUSOL_replaceColumn(LUSOL, randcol, x);
#else
      inform = LUSOL_replaceColumn(LUSOL, randcol, b);
#endif
    }
    else {
      /* Set the previously replaced column back and resolve */
      for (k = 1; k <= nnzero; k++) {
        i    = iA[k];
        j    = jA[k];
        if(j == randcol)
          x[i] = Aij[k];
      }
      inform = LUSOL_replaceColumn(LUSOL, randcol, x);
    }
    if(inform != LUSOL_INFORM_LUSUCCESS)
      fprintf(outunit, "Error:\n%s\n", LUSOL_informstr(LUSOL, inform));
    else
      goto Resolve;
  }


/* Free memory */
x900:
  if(!success)
    fprintf(outunit, "Insufficient memory or data file not found.\n");

  free(xexact);

}

void sp_error_message(const char *error_text)
{
     error_message( error_text );
}


DMatrix& null( const SparseMatrix& A )
{
      DMatrix* Af = new DMatrix;
      DMatrix* r  = new DMatrix;

      *Af = full(A);

      *r = null(*Af);

      delete Af;

      return(*r);
}

DMatrix& SVD(const SparseMatrix& A, DMatrix* U, DMatrix* V)
{
    DMatrix* Af= new DMatrix;
    DMatrix*  s= new DMatrix;

    *Af = full(A);

    *s = SVD(*Af, U, V);

    delete Af;

    return *s;

}

DMatrix& QR( const SparseMatrix& A )
{

    DMatrix* Af= new DMatrix;
    DMatrix*  s= new DMatrix;

    *Af = full(A);

    *s = QR(*Af);

    delete Af;

    return *s;


}

DMatrix& LQ( const SparseMatrix& A, DMatrix* Q )
{
    DMatrix* Af= new DMatrix;
    DMatrix*  s= new DMatrix;

    *Af = full(A);

    *s = LQ(*Af, Q);

    delete Af;

    return *s;


}


DMatrix& orth( const SparseMatrix& A )
{
    DMatrix* Af= new DMatrix;
    DMatrix*  s= new DMatrix;

    *Af = full(A);

    *s = orth(*Af);

    delete Af;

    return *s;

}

DMatrix& schur(const SparseMatrix& A, DMatrix* U )
{
    DMatrix* Af= new DMatrix;
    DMatrix*  s= new DMatrix;

    *Af = full(A);

    *s = schur(*Af, U);

    delete Af;

    return *s;

}

DMatrix& eig(const SparseMatrix& A, DMatrix* V  )
{
    DMatrix* Af= new DMatrix;
    DMatrix*  s= new DMatrix;

    *Af = full(A);

    *s = eig(*Af, V);

    delete Af;

    return *s;


}

double enorm(const SparseMatrix& A)
{
   int i;
   double r = 0;
   for (i=0;i<A.nz;i++)
       r+= pow(A.a[i],2.0);
   r = sqrt(r);

   return r;
}

double norm(const SparseMatrix& A)
{
    double r;

    DMatrix* Af= new DMatrix;

    *Af = full(A);

    r = norm(*Af);

    delete Af;

    return r;

}

SparseMatrix& Abs(const SparseMatrix& A)
{
    int i;

    SparseMatrix*  R= new SparseMatrix(A);

    for (i=0; i<R->nz; i++)
       R->a[i] = fabs(R->a[i]);

    return *R;
}

double cond( const SparseMatrix& A )
{
    double r;

    DMatrix* Af= new DMatrix;

    *Af = full(A);

    r = cond(*Af);

    delete Af;

    return r;

}

double rcond( const SparseMatrix& A )
{
    double r;

    DMatrix* Af= new DMatrix;

    *Af = full(A);

    r = rcond(*Af);

    delete Af;

    return r;

}

int rank_sparse( const SparseMatrix& A )
{
    int r;

    DMatrix* Af= new DMatrix;

    *Af = full(A);

    r = rank_matrix(*Af);

    delete Af;

    return r;

}

SparseMatrix& sprand(const SparseMatrix& S)
{
    int i;

    SparseMatrix*  R= new SparseMatrix(S);

    for (i=0; i<R->nz; i++)
       R->a[i] = DMatrix::random_uniform();

    return *R;
}

SparseMatrix& sprand(int n, int m, double density)
{
    int i, j;

    double r1;

    int nzguess = (int) density*n*m;

    SparseMatrix*  R= new SparseMatrix(n,m,nzguess);

    for (i=1; i<=n; i++)
    {
         for(j=1; j<=m; j++)
       {
               r1 = DMatrix::random_uniform();

               if (r1 < density)
               {
              R->InsertNonZero(i,j, DMatrix::random_uniform());
                }
         }
    }

    return *R;

}

SparseMatrix& sprandn(const SparseMatrix& S)
{
    int i;

    SparseMatrix*  R= new SparseMatrix(S);

    for (i=0; i<R->nz; i++)
       R->a[i] = DMatrix::random_gaussian();

    return *R;
}

SparseMatrix& sprandn(int n, int m, double density)
{
    int i, j;

    double r1;

    int nzguess = (int) density*n*m;

    SparseMatrix*  R= new SparseMatrix(n,m,nzguess);

    for (i=1; i<=n; i++)
    {
         for(j=1; j<=m; j++)
       {
               r1 = DMatrix::random_uniform();

               if (r1 < density)
               {
              R->InsertNonZero(i,j, DMatrix::random_gaussian());
                }
         }
    }

    return *R;
}
