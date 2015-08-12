// DMatrix class
// Example consisting of the iterative solution of the
// discrete Riccati algebraic equation
#include "dmatrixv.h"
int main() {

// Initialize the DMatrix objects

DMatrix A(3,3,      0.8,    1.0,  0.7 ,
                   -1.0,    0.5,  0.2 ,
                    0.1,    0.2,  0.4   );
DMatrix B( 3,2,    0.5, 0.2 ,
                   0.9, 0.4,
                   0.1, -0.5 );
DMatrix Q = identity(3);
DMatrix R  = 0.1*identity(2);
DMatrix S  = zeros(3,3);
DMatrix Sn(3,3);
double epsilon = 0.00001;
int iter_max = 200;
int j;
// Recursive loop
for ( j = 1; j<= iter_max; j++ ) 
{
        printf("\r iter = %d", j);
        // Riccati difference equation
        Sn= tra(A)*(S-S*B*inv(tra(B)*S*B+R)*tra(B)*S)*A+Q;
        // Check convergence
        if ( enorm( Sn - S ) <= epsilon )  break;
        // Update S
        S = Sn;
}
  
S.Print(" Solution S=");
return 0;

}
// Program ends
