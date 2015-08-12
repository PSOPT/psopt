#include "dmatrixv.h"


int main(void)
{

   DMatrix A, B;
 

   A = randu(10,1);
  
   A.Print("randu(10,1)");

   B = randn(10,1);

   B.Print("randn(10,1)");

   mean(B).Print("mean(B)");

   Std(B).Print("std(B)");

   return 0;

}    
