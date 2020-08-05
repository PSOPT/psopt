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

// Implementation of PSOPT_extras class function


#include "psopt.h"

DEC_THREAD int     PSOPT_extras::errorFlag   = 0;	   

DEC_THREAD int     PSOPT_extras::print_level = 1;        

DEC_THREAD clock_t PSOPT_extras::start_clock = 0;



void PSOPT_extras::SetPrintLevel( int plevel )
{
	 print_level=plevel;
	 return;
}


int PSOPT_extras::PrintLevel() {
	return print_level;
}

#ifndef OUTPUT_STREAM
#define OUTPUT_STREAM stderr
#endif

void error_message(const char *error_text)
{
// Error handling routine

  string m1;
  string m2;

  m1 =  "\n**** Run time error";
  m1 += "\n**** To trace this error, set up your debugger to break at";
  m1 += "\n**** function 'error_message', then do a backtrace.";
  m1 += "\n**** A diagnostic message is given below:";
  m2 = error_text;
  m2 =  "\n**** ====> " + m2 + " <====\n\n";

  if ( PSOPT_extras::PrintLevel() ) {

  fprintf(OUTPUT_STREAM,"%s", m1.c_str());
  fprintf(OUTPUT_STREAM,"%s", m2.c_str());
  FILE* err_file = fopen("error_message.txt","w");
  fprintf(err_file,"%s", m1.c_str() );
  fprintf(err_file,"%s", m2.c_str() );
  fclose(err_file);


 }

  PSOPT_extras::RiseErrorFlag();

  throw ErrorHandler(m1+m2);


}

ErrorHandler::ErrorHandler(const string m)
{
  error_message = m;
}


void PSOPT_extras::tic(void)
{
// Chronometer starts
   PSOPT_extras::SetStartTicks(clock());

}

double PSOPT_extras::toc()
{
// Chronometer stops

   clock_t  stop_ticks     = clock();
   clock_t  elapsed_ticks;
   double elapsed_time;

   elapsed_ticks= stop_ticks-PSOPT_extras::GetStartTicks();

   elapsed_time = ((double) elapsed_ticks/CLOCKS_PER_SEC);

   if (1) {

       if (PSOPT_extras::PrintLevel() ) {
       	fprintf(stderr,"\n Elapsed time is %e seconds\n", elapsed_time );
       }

   }

   return (elapsed_time);

}


