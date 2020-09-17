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


#ifndef __snoptProbLocal_H__
#define __snoptProbLocal_H__

#include "snoptProblem.hpp"




/** C++ interfacing with PSOPT.
 */
class snoptProbLocal : public snoptProblemA
{
public:
  /** default constructor */
  snoptProbLocal(Workspace* pr, void* user_data);

  /** default destructor */
  virtual ~snoptProbLocal();
  
  static void snPSOPTusrf_
  ( int    *Status, int *n,    double x[],
    int    *needF,  int *neF,  double F[],
    int    *needG,  int *neG,  double G[],
    char       *cu, int *lencu,
    int    iu[],    int *leniu,
    double ru[],    int *lenru );
  
  static Workspace *workspace;
  static void     *_user_data;  
  

  //@}

private:
  /**@name Methods to block default compiler methods.
   * The compiler automatically generates the following three methods.
   *  Since the default compiler implementation is generally not what
   *  you want (for all but the most simple classes), we usually
   *  put the declarations of these methods in the private section
   *  and never implement them. This prevents the compiler from
   *  implementing an incorrect "default" behavior without us
   *  knowing. (See Scott Meyers book, "Effective C++")
   *
   */
  //@{
  //  IPOPT_PSOPT();
  snoptProbLocal(const IPOPT_PSOPT&);
  snoptProbLocal& operator=(const IPOPT_PSOPT&);
  //@}

};

#endif

