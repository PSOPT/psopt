#include "dmatrixv.h"

// This file is provided as a template for testing the examples given
// in the DMatrix class manual

int main(void)
{

// Paste the example code within the box below
/* _________________________________________________________ */

    int* I[5] = {1, 1, 3, 5, 7};
    int* J[5] = {2, 3, 6, 9, 12};
    double* Aij[5] = {1.0, 2.0, 3.0, 4.0, 5.0}

    SparseMatrix A(5, Aij, I, J, 12, 12, 5);

    DMatrix Adense;

    Adense = A;

    A.Print("A");


/* _________________________________________________________ */

 return 0;
 
}    
