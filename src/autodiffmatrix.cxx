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

// Implementation of AutoDiffMatrix class functions




AutoDiffMatrix::AutoDiffMatrix(int ninput, int minput)
{
   n = ninput;
   m = minput;

   a = new adouble[n*m];
   if (a==NULL) error_message("Allocation error in AutoDiffMatrix::AutoDiffMatrix()");


}

AutoDiffMatrix::AutoDiffMatrix(int ninput)
{
   n = ninput;
   m = 1;

   a = new adouble[n*m];

}

AutoDiffMatrix::AutoDiffMatrix()
{
   n = 0;
   m = 0;

   a = NULL;

}

adouble& AutoDiffMatrix::operator() (int i, int j)
{
  if ( i<0 || j<0 || i>n-1 || j>m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[((j)*n+i) ];
}

adouble& AutoDiffMatrix::elem(int i, int j)
{
  if ( i<0 || j<0 || i>n-1 || j>m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[((j)*n+i) ];
}

adouble& AutoDiffMatrix::operator() (int i)
{
  if ( i<0 || i>n*m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[ i ];
}

adouble AutoDiffMatrix::operator() (int i, int j) const
{
  if ( i<0 || j<0 || i>n-1 || j>m-1 )   // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[((j)*n+i)];
}

adouble AutoDiffMatrix::elem(int i, int j) const
{
  if ( i<0 || j<0 || i>n-1 || j>m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[((j)*n+i)  ];
}

adouble AutoDiffMatrix::operator() (int i) const
{
  if ( i<0 || i>n*m-1 )  // EIGEN_UPDATE
     error_message("\nIndex error in AutoDiffMatrix::operator()");

  return a[ i ];
}

AutoDiffMatrix::~AutoDiffMatrix()
{
   delete[] a;
}


void inverse_ad(const AutoDiffMatrix& minput, AutoDiffMatrix* minv)
{
  int i,p,k,n,m;

  adouble pivot;

  n = minput.rows();
  m = minput.cols();

  if (n!=m) { error_message("\n Matrix must be square in inverse_ad()"); }

  for( i=0;i<n;i++) {
    for(k=0;k<n;k++) {
         (*minv)(i,k) = minput(i,k);
    }
  }

  adouble *a = minv->GetPr();

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

            if (  pivot==0.0) error_message("singular matrix in inverse_ad()");
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


    if (a[p*n+p]==0.0) error_message("singular matrix");

//    Temp->elem(p,p)=1/Temp->elem(p,p);
    a[p*n+p]=1/pivot;


  }
}


void product_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* AB)
{

  // Matrix product optimised for vector columnwise storage

  int  i,j,k,na,nb,mb;


  int  rowA, colB;

  adouble sum;

  na=A.rows();

  mb=B.cols();

  nb=B.rows();

  const adouble* apr = A.GetConstPr();

  const adouble* bpr = B.GetConstPr();

  if (A.cols() != B.rows() )

  { error_message("Incoherent matrix dimensions in product_ad()"); }

  for (i=0; i< na; i++)
  {


    for (j=0; j<mb; j++)
    {

      colB = j*B.rows();

      sum=0.0;

      rowA = i;

      for (k=0; k<nb; k++)
      {

          sum += apr[ rowA ]*bpr[ colB + k ];
          rowA += A.rows();

      }

      (*AB)(i,j) = sum;

    }

  }


}


void subtract_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* AmB)
{

  int i;

  int n = A.rows();

  int m = B.cols();

  int nelem = n*m;

  const adouble *b = B.GetConstPr();

  const adouble *a = A.GetConstPr();


  if (m != B.cols() || n != B.rows())

  { error_message("Incoherent matrix dimensions in subtract_ad()"); }



  for (i=0; i< nelem; i++)
  {

     AmB->GetPr()[i] = a[i] - b[i];

  }

}

void sum_ad(const AutoDiffMatrix& A,const AutoDiffMatrix& B, AutoDiffMatrix* ApB)
{
  int i;

  int n = A.rows();

  int m = B.cols();

  int nelem = n*m;

  const adouble *b = B.GetConstPr();

  const adouble *a = A.GetConstPr();


  if (m != B.cols() || n != B.rows())

  { error_message("Incoherent matrix dimensions in add_ad()"); }



  for (i=0; i< nelem; i++)
  {

     ApB->GetPr()[i] = a[i] + b[i];

  }
}

void AutoDiffMatrix::Print(const char* text) const
{
// Prints an AutoDiffMatrix
// Eg. A.Print(" A = ");

  int i,j;

  fprintf(stderr,"\n%s\n",text);
  for (i=0;i<n;i++)  // EIGEN_UPDATE
  {
    for (j=0;j<m;j++)
    {
      fprintf(stderr," %f ", (this->elem(i,j).value()) );
      fprintf(stderr,"\t");
    }
    fprintf(stderr,"\n");
  }



}

void AutoDiffMatrix::resize(int nn, int mm)
{
    if (a!=NULL) delete[] a;

    n = nn;
    m = mm;

    a = new adouble[n*m];
}

void AutoDiffMatrix::setZero()
{
    for(int i=0;i< n*m; i++) {
       a[i] = 0.0;
    }
}

void AutoDiffMatrix::setOne()
{
    for(int i=0;i< n*m; i++) {
       a[i] = 1.0;
    }
}