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



#ifndef MAX
#define MAX(a, b) ( (a)>(b)?  (a):(b) )
#endif
#ifndef MIN
#define MIN(a, b) ( (a)<(b)?  (a):(b) )
#endif
//


#undef max
#undef min
#undef abs


#ifndef WIN32
#include <stdlib.h>
#endif

#ifdef WIN32
extern "C" {
_CRTIMP  int * __cdecl errno(void) { static int i=0; return &i; };
}
#endif

#include "psopt.h"
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>


using namespace std;
using namespace Eigen;


void plot(const MatrixXd& xa, const MatrixXd& ya, const string& title, const char* xlabel, const char* ylabel, const char* legend, const char* terminal, const char* output)
{

         MatrixXd x = xa;

	      MatrixXd y = ya;

         FILE *gscript;

         double range_min = 0.001;

	      char legend_i[100];

         MatrixXd XY;

         if ( y.rows() < y.cols() )
         {


                   x = x.transpose().eval();
                   y = y.transpose().eval();
                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY << x, y;
         }
         else {

                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY << x, y;
         }

         int ny = MIN( y.rows(), y.cols() );

         int pos = 0;

         double MinY = Min(y) - 0.05*fabs(Min(y));

         double MaxY = Max(y) + 0.05*fabs(Min(y));

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         if (fabs(MinY-MaxY)< range_min) MaxY=MinY+range_min;

         Save(XY, "XY.dat");

         gscript = fopen("gnuplot.scp","w");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset style data lines");

         fprintf(gscript,"\nset multiplot");

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset xlabel '%s'", xlabel);

         fprintf(gscript,"\nset ylabel '%s'", ylabel);

         fprintf(gscript,"\nset key box");

         fprintf(gscript,"\nset key outside");

         fprintf(gscript,"\nplot [ ] [%f:%f] ",MinY, MaxY);

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,"'XY.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with linespoints");

                if(i<(ny+1) ) fprintf(gscript,", ");

         }

         fprintf(gscript,"\nunset multiplot");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");



}


void multiplot(const MatrixXd& xa, const MatrixXd& ya, const string& title, const char* xlabel, const char* ylabel, const char* legend, int nrows, int ncols, const char* terminal,  const char* output )
{
         MatrixXd x = xa;

	      MatrixXd y = ya;

         FILE *gscript;

         double range_min = 0.001;

	      char legend_i[100];
   	   char ylabel_i[100];

         MatrixXd XY, xcopy, ycopy;

         xcopy = x;

         ycopy = y;

         if ( y.rows() < y.cols() )
         {

                   x = x.transpose().eval();
                   y = y.transpose().eval();
                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY<< x, y;
                   xcopy = x;
                   ycopy = y;
         }
         else {

                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY << x, y;
         }

         int ny = MIN( y.rows(), y.cols() );

         int pos = 0;
         int pos_y = 0;

         int nnrows = nrows;
         int nncols  = ncols;

         if (nrows==0) nnrows = ny;
         if (ncols ==0) nncols = 1;

         if (nnrows*nncols != ny)
             error_message("multiplot() error: number of independent variables to plot must correspond with the grid size nrows*ncols");


         Save(XY, "XY.dat");

         gscript = fopen("gnuplot.scp","w");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset style data lines");

         fprintf(gscript,"\nset multiplot layout %i,%i columnsfirst", nnrows,nncols);

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset xlabel '%s'", xlabel);

         fprintf(gscript,"\nset key box");

         fprintf(gscript,"\nset key outside");

         for(int i=2;i<=ny+1;i++) {


            MatrixXd yi = ycopy.col(i-2);

	         double MinY = Min( yi ) - 0.05*fabs(Min( yi ));

	         double MaxY = Max( yi ) + 0.05*fabs(Max( yi ));

	         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
	         {
	             MinY = -range_min;
	             MaxY =  range_min;
	         }



                 if (ylabel!=NULL) {

                  sscanf(&ylabel[pos_y], "%s", ylabel_i);

                  pos_y+= strlen(ylabel_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(ylabel_i,"%i",i-1);
                   else
                     sprintf(ylabel_i,"");
                }

                fprintf(gscript,"\nset ylabel '%s'", ylabel_i);

                fprintf(gscript,"\nplot [ ] [%f:%f] ",MinY, MaxY);

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,"'XY.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with linespoints");


         }

         fprintf(gscript,"\nunset multiplot");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");



}

void plot(const MatrixXd& x1a, const MatrixXd& y1a, const MatrixXd& x2a, const MatrixXd& y2a, const string& title, const char* xlabel, const char* ylabel, const char* legend, const char* terminal, const char* output)
{

         MatrixXd x1 = x1a;
	 		MatrixXd y1 = y1a;
	 		MatrixXd x2 = x2a;
	 		MatrixXd y2 = y2a;

         FILE *gscript;

         double range_min = 0.001;

	 		char legend_i[100];

         MatrixXd XY;

         if ( y1.rows() < y1.cols() )
         {

                   x1 = x1.transpose().eval();
                   y1 = y1.transpose().eval();
                   XY.resize(x1.rows(), x1.cols()+y1.cols());
                   XY << x1,y1;
         }
         else {

                  XY.resize(x1.rows(), x1.cols()+y1.cols());
                  XY << x1, y1;
         }

         int ny = MIN( y1.rows(), y1.cols() );

         int pos = 0;

         double MinY = MIN ( Min(y1), Min(y2) );

         double MaxY = MAX ( Max(y1), Max(y2) );

         MinY -= 0.05*fabs(MinY);
         MaxY += 0.05*fabs(MaxY);

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY, "XY1.dat");

         gscript = fopen("gnuplot.scp","w");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset style data lines");

         fprintf(gscript,"\nset multiplot");

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset xlabel '%s'", xlabel);

         fprintf(gscript,"\nset ylabel '%s'", ylabel);

         fprintf(gscript,"\nset key box");

         fprintf(gscript,"\nset key outside");

         fprintf(gscript,"\nplot [ ] [%f:%f] ",MinY, MaxY);

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,"'XY1.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with linespoints");

                if(i<=(ny+1) ) fprintf(gscript,", ");

         }

         if ( y2.rows() < y2.cols() )
         {

                   x2 = x2.transpose().eval();
                   y2 = y2.transpose().eval();
                   XY.resize(x2.rows(), x2.cols()+y2.cols());
                   XY<< x2, y2;
         }
         else {

                 XY.resize(x2.rows(), x2.cols()+y2.cols());
                 XY << x2 , y2;

         }

         ny = MIN( y2.rows(), y2.cols());




         Save(XY, "XY2.dat");

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {
                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,"'XY2.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with linespoints");

                if(i<(ny+1) ) fprintf(gscript,", ");

         }



         fprintf(gscript,"\nunset multiplot");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");



}


void spplot(const MatrixXd& x1a, const MatrixXd& y1a, const MatrixXd& x2a, const MatrixXd& y2a, const string& title, const char* xlabel, const char* ylabel, const char* legend, const char* terminal, const char* output)
{

         MatrixXd x1 = x1a;
	     	MatrixXd y1 = y1a;
	     	MatrixXd x2 = x2a;
	     	MatrixXd y2 = y2a;

         FILE *gscript;

         double range_min = 0.001;

	     	char legend_i[100];

         MatrixXd XY;

         if ( y1.rows() < y1.cols() )
         {

                   x1 = x1.transpose().eval();
                   y1 = y1.transpose().eval();
                   XY.resize(x1.rows(), x1.cols()+y1.cols());
                   XY<< x1,y1 ;
         }
         else {

                  XY.resize(x1.rows(), x1.cols()+y1.cols());
                  XY << x1, y1;
         }

         int ny = MIN( y1.rows(), y1.cols() );

         int pos = 0;

         double MinY = MIN ( Min(y1), Min(y2) );

         double MaxY = MAX ( Max(y1), Max(y2) );

         MinY -= 0.05*fabs(MinY);
         MaxY += 0.05*fabs(MaxY);

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY, "XY1.dat");

         gscript = fopen("gnuplot.scp","w");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset style data lines");

         fprintf(gscript,"\nset multiplot");

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset xlabel '%s'", xlabel);

         fprintf(gscript,"\nset ylabel '%s'", ylabel);

         fprintf(gscript,"\nset key box");

         fprintf(gscript,"\nset key outside");

         fprintf(gscript,"\nplot [ ] [%f:%f] ",MinY, MaxY);

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,"'XY1.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with linespoints");

                if(i<=(ny+1) ) fprintf(gscript,", ");

         }

         if ( y2.rows() < y2.cols() )
         {

                   x2 = x2.transpose().eval();
                   y2 = y2.transpose().eval();
                   XY.resize(x2.rows(), x2.cols()+y2.cols());
                   XY << x2, y2;
         }
         else {

                 XY.resize(x2.rows(), x2.cols()+y2.cols());
                 XY << x2, y2;
         }

         ny = MIN( y2.rows(), y2.cols());



         Save(XY, "XY2.dat");

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {
                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,"'XY2.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with points");

                if(i<(ny+1) ) fprintf(gscript,", ");

         }



         fprintf(gscript,"\nunset multiplot");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");



}



void plot(const MatrixXd& x1a, const MatrixXd& y1a, const MatrixXd& x2a, const MatrixXd& y2a, const MatrixXd& x3a, const MatrixXd& y3a,
          const string& title, const char* xlabel, const char* ylabel, const char* legend, const char* terminal, const char* output)
{
         MatrixXd x1 = x1a;
	      MatrixXd y1 = y1a;
	      MatrixXd x2 = x2a;
	      MatrixXd y2 = y2a;
	      MatrixXd x3 = x3a;
	      MatrixXd y3 = y3a;

         FILE *gscript;

         double range_min = 0.001;

	      char legend_i[100];

         MatrixXd XY;

         if ( y1.rows() < y1.cols() )
         {

                   x1 = x1.transpose().eval();
                   y1 = y1.transpose().eval();
                   XY.resize(x1.rows(), x1.cols()+y1.cols());
                   XY<< x1,y1 ;
         }
         else {

                XY.resize(x1.rows(), x1.cols()+y1.cols());
                XY << x1 , y1;
         }

         int ny = MIN( y1.rows(), y1.cols() );

         int pos = 0;

         double MinY = MIN ( Min(y1), Min(y2) );

         double MaxY = MAX ( Max(y1), Max(y2) );

         MinY -= 0.05*fabs(MinY);
         MaxY += 0.05*fabs(MaxY);

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY,"XY1.dat");

         gscript = fopen("gnuplot.scp","w");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset style data lines");

         fprintf(gscript,"\nset multiplot");

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset xlabel '%s'", xlabel);

         fprintf(gscript,"\nset ylabel '%s'", ylabel);

         fprintf(gscript,"\nset key box");

         fprintf(gscript,"\nset key outside");

         fprintf(gscript,"\nplot [ ] [%f:%f] ",MinY, MaxY);

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,"'XY1.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with linespoints");

                if(i<(ny+1) ) fprintf(gscript,", ");

         }

         if ( y2.rows() < y2.cols() )
         {

                   x2 = x2.transpose().eval();
                   y2 = y2.transpose().eval();
                   XY.resize(x2.rows(), x2.cols()+y2.cols());
                   XY<< x2,y2 ;
         }
         else {

                 XY.resize(x2.rows(), x2.cols()+y2.cols());
                 XY << x2 , y2;
         }

         ny = MIN( y2.rows(), y2.cols());



         Save(XY,"XY2.dat");

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {
                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,", 'XY2.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with linespoints");

                if(i<(ny+1) ) fprintf(gscript,", ");

         }


         if ( y3.rows() < y3.cols() )
         {

                   x3 = x3.transpose().eval();
                   y3 = y3.transpose().eval();
                   XY.resize(x3.rows(), x3.cols()+y3.cols());
                   XY<< x3,y3 ;
         }
         else {

               XY.resize(x3.rows(), x3.cols()+y3.cols());
               XY << x3, y3;
         }

         ny = MIN( y3.rows(), y3.cols());




         Save(XY,"XY3.dat");

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {
                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }
         	fprintf(gscript,", 'XY3.dat' using 1:%i title \"%s\" ",  i, legend_i);

                fprintf(gscript,"with linespoints");

                if(i<(ny+1) ) fprintf(gscript,", ");

         }



         fprintf(gscript,"\nunset multiplot");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");



}


void polar(const MatrixXd& theta_a, const MatrixXd& r_a, const string& title,  const char* legend, const char* terminal, const char* output)
{
         MatrixXd theta = theta_a;
	 		MatrixXd r = r_a;

         FILE *gscript;

         double range_min = 0.001;

	 		char legend_i[100];

         MatrixXd & x = theta;

         MatrixXd & y = r;

         MatrixXd XY;

         if ( y.rows() < y.cols() )
         {

                   x = x.transpose().eval();
                   y = y.transpose().eval();
                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY<< x,y ;
         }
         else {

                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY << x,y;
         }

         int ny = MIN( y.rows(), y.cols() );

         int pos = 0;

         double MinY = Min(y);

         double MaxY = Max(y);

         MinY -= 0.05*MinY;
         MaxY += 0.05*MaxY;

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY,"XY.dat");

         gscript = fopen("gnuplot.scp","w");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset style data lines");

         fprintf(gscript,"\nset multiplot");

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset polar");

         fprintf(gscript,"\nplot [ ] ");

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }

         	fprintf(gscript,"'XY.dat' using 1:%i title \"%s\" ",  i, legend_i);


                if(i<(ny+1) ) fprintf(gscript,", ");

         }

         fprintf(gscript,"\nunset multiplot");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");



}


void polar(const MatrixXd& theta_a, const MatrixXd& r_a, const MatrixXd& theta2_a, const MatrixXd& r2_a, const string& title,  const char* legend, const char* terminal, const char* output)
{
         MatrixXd theta = theta_a;
	 		MatrixXd r = r_a;
	 		MatrixXd theta2 = theta2_a;
	 		MatrixXd r2 = r2_a;

         FILE *gscript;

         double range_min = 0.001;

	   	char legend_i[100];

         MatrixXd  x = theta;

         MatrixXd  y = r;

         MatrixXd XY;

         if ( y.rows() < y.cols() )
         {

                   x = x.transpose().eval();
                   y = y.transpose().eval();
                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY<< x,y ;
         }
         else {

                 XY.resize(x.rows(), x.cols()+y.cols());
                 XY << x,y;
         }

         int ny = MIN( y.rows(), y.cols() );

         int pos = 0;

         double MinY = Min(y);

         double MaxY = Max(y);

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY,"XY.dat");

         gscript = fopen("gnuplot.scp","w");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset style data lines");

         fprintf(gscript,"\nset multiplot");

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset polar");

         fprintf(gscript,"\nplot [ ] ");

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }

         	fprintf(gscript,"'XY.dat' using 1:%i title \"%s\" ",  i, legend_i);


                if(i<(ny+1) ) fprintf(gscript,", ");

         }

         x = theta2;

         y = r2;


         if ( y.rows() < y.cols() )
         {

                   x = x.transpose().eval();
                   y = y.transpose().eval();
                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY<< x,y ;
         }
         else {

                 XY.resize(x.rows(), x.cols()+y.cols());
                 XY << x,y;
         }

         ny = MIN( y.rows(), y.cols() );

         pos = 0;

         MinY = Min(y);

         MaxY = Max(y);

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY,"XY2.dat");



         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }

         	fprintf(gscript,", 'XY2.dat' using 1:%i title \"%s\" ",  i, legend_i);

                if(i<(ny+1) ) fprintf(gscript,", ");

         }

         fprintf(gscript,"\nunset multiplot");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");



}


void polar(const MatrixXd& theta_a, const MatrixXd& r_a, const MatrixXd& theta2_a, const MatrixXd& r2_a, const MatrixXd& theta3_a, const MatrixXd& r3_a, const string& title,  const char* legend, const char* terminal, const char* output)
{

         MatrixXd theta = theta_a;
	 		MatrixXd r = r_a;

	 		MatrixXd theta2 = theta2_a;
	 		MatrixXd r2 = r2_a;

	 		MatrixXd theta3 = theta3_a;
	 		MatrixXd r3 = r3_a;

         FILE *gscript;

         double range_min = 0.001;

	 		char legend_i[100];

         MatrixXd  x = theta;

         MatrixXd  y = r;

         MatrixXd XY;

         if ( y.rows() < y.cols() )
         {

                   x = x.transpose().eval();
                   y = y.transpose().eval();
                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY<< x,y ;
         }
         else {

                 XY.resize(x.rows(), x.cols()+y.cols());
                 XY << x,y;
         }

         int ny = MIN( y.rows(), y.cols() );

         int pos = 0;

         double MinY = Min(y);

         double MaxY = Max(y);

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY,"XY.dat");

         gscript = fopen("gnuplot.scp","w");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset style data lines");

         fprintf(gscript,"\nset multiplot");

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset polar");

         fprintf(gscript,"\nplot [ ] ");

         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }

         	fprintf(gscript,"'XY.dat' using 1:%i title \"%s\" ",  i, legend_i);


                if(i<(ny+1) ) fprintf(gscript,", ");

         }

         x = theta2;

         y = r2;


         if ( y.rows() < y.cols() )
         {

                   x = x.transpose().eval();
                   y = y.transpose().eval();
                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY<< x,y ;
         }
         else {

                 XY.resize(x.rows(), x.cols()+y.cols());
                 XY << x,y;
         }

         ny = MIN( y.rows(), y.cols() );

         pos = 0;

         MinY = Min(y);

         MaxY = Max(y);

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY,"XY2.dat");


         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }

         	fprintf(gscript,", 'XY2.dat' using 1:%i title \"%s\" ",  i, legend_i);

                if(i<(ny+1) ) fprintf(gscript,", ");

         }

         x = theta3;

         y = r3;


         if ( y.rows() < y.cols() )
         {

                   x = x.transpose().eval();
                   y = y.transpose().eval();
                   XY.resize(x.rows(), x.cols()+y.cols());
                   XY<< x,y ;
         }
         else {

                 XY.resize(x.rows(), x.cols()+y.cols());
                 XY << x,y;
         }

         ny = MIN( y.rows(), y.cols() );

         pos = 0;

         MinY = Min(y);

         MaxY = Max(y);

         if (fabs(MinY)<range_min && fabs(MaxY)<range_min )
         {
             MinY = -range_min;
             MaxY =  range_min;
         }

         Save(XY,"XY3.dat");



         for(int i=2;i<=ny+1;i++) {

                if (legend!=NULL) {

                  sscanf(&legend[pos], "%s", legend_i);

                  pos+= strlen(legend_i)+1;

                }

                else {
                   if (ny>1)
                     sprintf(legend_i,"%i",i-1);
                   else
                     sprintf(legend_i,"");
                }

         	fprintf(gscript,", 'XY3.dat' using 1:%i title \"%s\" ",  i, legend_i);

                if(i<(ny+1) ) fprintf(gscript,", ");

         }

         fprintf(gscript,"\nunset multiplot");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");



}






void surf(const MatrixXd& xa, const MatrixXd& ya, const MatrixXd& za, const string& title, const char* xlabel, const char* ylabel, const char* zlabel, const char* terminal, const char* output, const char* view)
{
	 // This function creates surface plots given the co-ordinate values (x,y) and the height matrix z.

	 		MatrixXd x = xa;
	 		MatrixXd y = ya;
			MatrixXd z = za;

         FILE *gscript, *datafile;

         double range_min = 0.001;

         int i,j;

         MatrixXd X, Y;

         if (x.rows()>1 && x.cols()>1) error_message("surf(): MatrixXd object x must be a vector");

         if (y.rows()>1 && y.cols()>1) error_message("surf(): MatrixXd object y must be a vector");


           X = stack_columns(x);


           Y = stack_columns(y);

         if ( z.rows()*z.cols() != length(X)*length(y) ) {
				error_message("surf(): input MatrixXd object z has inconsistent dimensions with length(x) and length(y)");
         }

         int nz = MIN( z.rows(), z.cols() );

         int pos = 0;

         double MinZ = Min(z);

         double MaxZ = Max(z);

         if (fabs(MinZ)<range_min && fabs(MaxZ)<range_min )
         {
             MinZ = -range_min;
             MaxZ =  range_min;
         }

         datafile = fopen("XYZ.dat","w");

         if (datafile==NULL) error_message("surf(): error creating gnuplot data file");

         for (i=0; i<length(X); i++) {
				for(j=0;j<length(Y);j++) {
                    fprintf(datafile,"%f  %f  %f\n",X(i),Y(j), z(i,j) );
	         } 
		      fprintf(datafile,"\n");
	      }

         fclose(datafile);

         gscript = fopen("gnuplot.scp","w");

         if (gscript==NULL) error_message("surf(): error creating gnuplot script file");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset xlabel '%s'", xlabel);

         fprintf(gscript,"\nset ylabel '%s'", ylabel);

         fprintf(gscript,"\nset zlabel '%s'", zlabel);

         fprintf(gscript,"\nunset key");

         fprintf(gscript,"\nset hidden3d");

         fprintf(gscript,"\nset grid");

	 if (view != NULL) {
             fprintf(gscript,"\nset view %s", view);
	 }

         fprintf(gscript,"\nset zrange [%f:%f] ",MinZ, MaxZ);

         fprintf(gscript,"\nsplot 'XYZ.dat' with lines");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");

}

void plot3(const MatrixXd& xa, const MatrixXd& ya, const MatrixXd& za, const string& title, const char* xlabel, const char* ylabel, const char* zlabel, const char* terminal, const char* output, const char* view)
{
	 // This function creates 3d plots given the co-ordinate values (x,y) and the height vector z.

	 		MatrixXd x = xa;
	 		MatrixXd y = ya;
	 		MatrixXd z = za;

         FILE *gscript, *datafile;

         double range_min = 0.001;

         int i,j;

         MatrixXd X, Y, Z;

         if (x.rows()>1 && x.cols()>1) error_message("plot3(): MatrixXd object x must be a vector");

         if (y.rows()>1 && y.cols()>1) error_message("plot3(): MatrixXd object y must be a vector");

         if (z.rows()>1 && z.cols()>1) error_message("plot3(): MatrixXd object z must be a vector");


         X = stack_columns(x); 


         Y = stack_columns(y);


         Z = stack_columns(z);

         if ( length(X) != length(Y) || length(x)!=length(Z) ) {
		error_message("plot3(): inconsistent array lengths");
         }

         int nz = MIN( z.rows(), z.cols() );

         int pos = 0;

         double MinZ = Min(z);

         double MaxZ = Max(z);

         if (fabs(MinZ)<range_min && fabs(MaxZ)<range_min )
         {
             MinZ = -range_min;
             MaxZ =  range_min;
         }

         datafile = fopen("XYZ.dat","w");

         if (datafile==NULL) error_message("surf(): error creating gnuplot data file");

         for (i=0; i<length(X); i++) {

               fprintf(datafile,"%f  %f  %f\n",X(i),Y(i), Z(i) );

	      }

         fclose(datafile);

         gscript = fopen("gnuplot.scp","w");

         if (gscript==NULL) error_message("surf(): error creating gnuplot script file");

         if (terminal!=NULL) fprintf(gscript,"\nset terminal %s", terminal);

         if (output!=NULL) fprintf(gscript,"\nset output \"%s\"", output);

         fprintf(gscript,"\nset title '%s'", title.c_str());

         fprintf(gscript,"\nset xlabel '%s'", xlabel);

         fprintf(gscript,"\nset ylabel '%s'", ylabel);

         fprintf(gscript,"\nset zlabel '%s'", zlabel);

         fprintf(gscript,"\nunset key");

         fprintf(gscript,"\nset hidden3d");

      	if (view != NULL) {
             fprintf(gscript,"\nset view %s", view);
	      }

         fprintf(gscript,"\nset grid");

         fprintf(gscript,"\nset zrange [%f:%f] ",MinZ, MaxZ);

         fprintf(gscript,"\nsplot 'XYZ.dat' with lines");

         fclose(gscript);

         system("gnuplot -persist gnuplot.scp ");

}






