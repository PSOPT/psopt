// colon.cxx

// V.M. Becerra, September 1996

// Program to test some indexing functions of DMatrix


#include "dmatrixv.h"

#include <stdio.h>


int main()

{



    DMatrix A( 5, 5 );



    DMatrix B( 4, 5 );



    DMatrix C( 5, 5 );



    DMatrix D(2,2, 1.0, 2.0,

                   3.0, 4.0 );


	DMatrix E;

	DMatrix F;


    A = eye( 5 );



    B = A( colon(1,4), colon() );



    C( colon(1,4), colon(1,4 ) ) = eye(4);



    C( colon(5,5), colon(5,5) ) = 10.0;



    A.Print("A");



    B.Print("A(1:4,:)");



    C.Print("C");



    printf("\nPress ENTER");getchar();



    D.Print("D");



    D( colon() ).Print("D(:)");



    D( colon(4,-1,1) ).Print("D(4:-1:1)");



    colon(0.0, 1.5, 9.0 ).Print("0:1.5:9.0");



    printf("\nPress ENTER");getchar();



    A(1, colon()) = ones(1, A.GetNoCols() );

    A(colon(), 2 ) = 2.0*ones(A.GetNoRows(),1);


    A.Print(" A(1,:) = ones() ; A(:,2)= 2*ones() ");

	 printf("\nPress ENTER");getchar();

    // Now grow the matrix by one column...

	 A(colon(),6) = ones(5,1);

	 A.Print("A with additional column");


	 printf("\nPress ENTER");getchar();

	 E = colon(1,16);

         E= reshape(E,4,4);

	 E.Print("E");

	 F = colon(1,4);

	 F.Print("1:4=");

	 E(F).Print("E(1:4)");

	 return 0;

}





