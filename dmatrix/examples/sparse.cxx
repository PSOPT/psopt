#include "dmatrixv.h"

// This file is provided as a template for testing the examples given
// in the DMatrix class manual

int main(void)
{

// Paste the example code within the box below
/* _________________________________________________________ */

    int I[6] = {1, 1, 2, 2, 3, 3};
    int J[6] = {1, 2, 3, 1, 2, 3};
    double Aij[6] = { -	0.12264700E+04,
          		0.82867600E+01,
          		0.62150700E+01,
          		0.19468000E+03,
         	       -0.23990100E+04,
         	       -0.32849100E+04,
 		     };

    SparseMatrix A(Aij, 3, 3, 6, I, J);

    SparseMatrix B, C;

    DMatrix b(3,1, -0.187248E+05,
                    0.112479E+05,
        	    0.132356E+05 );

    DMatrix r(3);

    DMatrix x(3);

    void* LUSOLv;


    A.Print("A=");
    
    SparseLUSolve(A, b, x, r);

    x.Print("x");
    r.Print("r");

    LUSOLv = SparseLUFactor(A);

    SparseLUSolveGivenFactor(LUSOLv, A, b, x, r);

    x.Print("x");
    r.Print("r");

    B=1.0*A;

 //   B(2,2)= 10.0;
    B.InsertNonZero(2,2, 10.0);

    B.Print("B=");

    C = A - B;

    C.Print("C=A-B");

    SparseMatrix D = speye(10);

    D.Print("speye(10)");

    SparseMatrix E = spones(A);

    E.Print("spones(A)");

    DMatrix F(4,3);

    F(1,1) = 1; F(1,2)=1; F(1,3)=10.0;
    F(2,1) = 2; F(2,2)=2; F(2,3)=20.0;
    F(3,1) = 3; F(3,2)=3; F(3,3)=30.0;
    F(4,1) = 4; F(4,2)=4; F(4,3)=0.0;

    F.Print("F");
 
    SparseMatrix G = spconvert(F);

    G.Print("spconvert(F)");

    SparseMatrix H;

    H = A*x;

    H.Print("A*x");

    DMatrix fullA;

    fullA = full(A);

    fullA.Print("full(A)");

    A.Save("A.dat");

    SparseMatrix A2;

    A2.Load("A.dat");

    A2.Print("Read A.dat");

    SparseMatrix At = tra(A);

    At.Print("tra(A)");

    SparseMatrix AAt;

    AAt = A*At;

    AAt.Print("A*At");

    (fullA*tra(fullA)).Print("full A*At");

    elemProduct(A,A).Print("elemProduct(A,A)");

    (A^2.0).Print("A^2");

    SparseMatrix invA;

    invA = inv(A);

    (invA).Print("inv(A)");

    (inv(fullA)).Print("inv( full(A))");

    (A%A).Print("A%A");

    (A/A).Print("A/A");

    (A||A).Print("A || A");

    (A&&A).Print("A && A");

    SparseMatrix S(6,6,0);

    S.set_sub_matrix(A, 4, 6, 4, 6);

    S.Print("S");

    SparseMatrix Ssub = S.sub_matrix(4,6,4,6);

    Ssub.Print("S(4:6,4:6)");

    null(A).Print("null(A)");

    null(fullA).Print("null(fullA)");

    SVD(A).Print("SVD(A)");

    SparseMatrix R1 = sprand(10,10,0.1);

    SparseMatrix R2 = sprandn(10,10,0.1);

    R1.Print("sprand(10,10,0.1)");

    R2.Print("sprandn(10,10,0.1)");

    fprintf(stderr, "\nrank(R1)=%i\n", rank(R1));

    SparseMatrix A1000=sprand(1000,1000,0.02);
    SparseMatrix B1000=sprand(1000,1000,0.02);

    SparseMatrix AB1000;
    tic();
    AB1000 = A1000*B1000;
    toc();


/* _________________________________________________________ */

 return 0;
 
}    
