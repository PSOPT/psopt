//
//  helper.h
//  dmatrix
//
//  Created by Philipp Waxweiler on 22.07.20.
//

#ifndef helper_h
#define helper_h

#ifndef MAX
#define MAX(a, b) ( (a)>(b)?  (a):(b) )
#endif
#ifndef MIN
#define MIN(a, b) ( (a)<(b)?  (a):(b) )
#endif

#ifdef MATLAB_MEX_FILE

#define ERROR_MESSAGE error_message

#define PRINTF mexPrintf


#else

#define ERROR_MESSAGE error_message

#define PRINTF printf

#endif

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

inline long ChkTmpIndx( long taindx )
{
   if ( taindx >= N_TEMP_OBJECTS )

     ERROR_MESSAGE(" Temporary arrays error: Increase N_TEMP_OBJECTS");

#ifdef DEBUG_TEMPS

   printf("\n Temp Created --> indx: %d", taindx );

#endif

   return  taindx;
}



#endif /* helper_h */
