/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2026 Victor M. Becerra

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
           School of Electrical and Mechanical Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    vmbecerra@vmb1.com

**********************************************************************************************/


#include "psopt.h"

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;


using namespace Eigen;


// The rule that turns one variable's bounds into its map onto the scaled problem.
//
// Written once because it was written four times -- controls, states, parameters and
// time -- in four copies that had already begun to differ only in the names of their
// locals. Time keeps its own copy below because it has no shift: the collocation already
// maps [t0,tf] onto a fixed interval, and moving the origin of a variable that the
// transcription rescales for its own purposes buys nothing and would have to be undone
// in convert_to_original_time.
//
// "multiplicative" is PSOPT's historical rule and the default: the factor is one over
// the largest magnitude the variable is allowed, so the scaled variable is mostly inside
// [-1,1] and the origin does not move. "affine" puts a two-sidedly bounded variable onto
// exactly [-1,1], which fixes the case the multiplicative rule cannot: a variable bounded
// between 1000 and 1200 has unit magnitude under the first rule and a variation of 0.17.
//
// The interval is [-1,1] and not the [-1/2,1/2] of CGPOPS, and the difference is not
// cosmetic. A box that is already symmetric about zero has nothing for a shift to fix,
// and on [-a,a] this rule gives exactly the factor 1/a and the shift 0 that the
// multiplicative rule gives -- so "affine" changes only the variables it has something to
// say about. Mapping onto [-1/2,1/2] instead halves the factor on every symmetric box as
// well, which is a rescaling of problems the option was not meant to touch: measured that
// way over the shipped examples it cost a geometric mean of 1.09 in iterations, with the
// regressions concentrated on problems whose boxes are symmetric.
//
// A variable without two finite bounds has no finite centre, so there is nothing to shift
// to and the affine rule falls back to the multiplicative one. So does a variable whose
// bounds coincide, which has zero width.
//
// An absent bound is recognised through PSOPT::no_lower_bound rather than by comparing
// against an IEEE infinity. The two differ: the convention every model in examples/ uses
// for "no bound" is 1.0e19, and a bound written that way reached the comparisons here as
// finite, so a state declared unbounded the documented way was scaled by 1e-19. No
// shipped example writes a variable bound that way -- the convention is used for path
// bounds, which patch 157 dealt with -- so this is a trap rather than a live defect, but
// it is the same trap and it is closed here.
static void variable_map_from_bounds(double zlower, double zupper,
                                     const std::string& mode,
                                     double& scale, double& shift)
{
    scale = 1.0;
    shift = 0.0;

    const bool lo_absent = PSOPT::no_lower_bound(zlower);
    const bool up_absent = PSOPT::no_upper_bound(zupper);

    if ( mode == "affine" && !lo_absent && !up_absent && zupper > zlower ) {
        scale = 2.0/(zupper - zlower);
        shift = 0.5*(zlower + zupper);
        return;
    }

    if ( !lo_absent && !up_absent ) {
        if ( zlower != 0.0 || zupper != 0.0 )
            scale = 1.0/std::max( fabs(zlower), fabs(zupper) );
    }
    else if ( lo_absent && !up_absent && zupper != 0.0 )
        scale = 1.0/fabs(zupper);
    else if ( up_absent && !lo_absent && zlower != 0.0 )
        scale = 1.0/fabs(zlower);
}


void determine_scaling_factors_for_variables(Sol& solution, Prob& problem, Alg& algorithm)
{
     // Scaling factors  for variables computed automatically given the bound information
     // supplied by the user, such that the scaled variables are mostly in the interval [-1,1].
     // If, however, any of the variable bounds is 'inf' in magnitude, then the scaled range
     // will be [-1, inf], [-inf,1],  or [-1, 1]

   int i;

   for(i=0; i<problem.nphases; i++)
   {

	int norder    = problem.phase[i].current_number_of_intervals;
	int ncontrols = problem.phase[i].ncontrols;
	int npath     = problem.phase[i].npath;
	int nstates   = problem.phase[i].nstates;
	int nevents   = problem.phase[i].nevents;
        int nparam    = problem.phase[i].nparameters;

	MatrixXd& control_scaling = problem.phase[i].scale.controls;
	MatrixXd& state_scaling   = problem.phase[i].scale.states;
   MatrixXd& param_scaling   = problem.phase[i].scale.parameters;
   MatrixXd& control_shift   = problem.phase[i].scale.controls_shift;
   MatrixXd& state_shift     = problem.phase[i].scale.states_shift;
   MatrixXd& param_shift     = problem.phase[i].scale.parameters_shift;


	MatrixXd PathJac(npath, nstates+ncontrols);
	MatrixXd h1(npath,norder+1);
	MatrixXd h2(npath,norder+1);
	MatrixXd e1(nevents,1);
	MatrixXd e2(nevents,2);
	MatrixXd EventJac1(nevents, nstates+ncontrols);
	MatrixXd EventJac2(nevents, nstates+ncontrols);

	double zlower;
	double zupper;

	int ii;


	// Control scaling:

	if ( algorithm.scaling=="automatic" || algorithm.scaling!="user" )
	{
	   control_scaling = ones(ncontrols,1); // EIGEN_UPDATE
	   control_shift   = zeros(ncontrols,1);

	   for(ii=0;ii<ncontrols;ii++)
	   {
		zlower = (problem.phase[i].bounds.lower.controls)(ii);
		zupper = (problem.phase[i].bounds.upper.controls)(ii);
		variable_map_from_bounds(zlower, zupper, algorithm.scaling,
		                         control_scaling(ii), control_shift(ii));
	   }
	}




	// State scaling:

	if ( algorithm.scaling=="automatic" || algorithm.scaling!="user" )
	{

		state_scaling = ones(nstates,1);
		state_shift   = zeros(nstates,1);

		for(ii=0;ii<nstates;ii++) // EIGEN_UPDATE
		{
			zlower = (problem.phase[i].bounds.lower.states)(ii);
			zupper = (problem.phase[i].bounds.upper.states)(ii);
			variable_map_from_bounds(zlower, zupper, algorithm.scaling,
			                         state_scaling(ii), state_shift(ii));
		}
	}



	// Parameter scaling

	if ( algorithm.scaling=="automatic" || algorithm.scaling!="user" )
	{
		param_scaling = ones(nparam,1);
		param_shift   = zeros(nparam,1);

		for(ii=0;ii<nparam;ii++)  // EIGEN_UPDATE
		{
			zlower = (problem.phase[i].bounds.lower.parameters)(ii);
			zupper = (problem.phase[i].bounds.upper.parameters)(ii);
			variable_map_from_bounds(zlower, zupper, algorithm.scaling,
			                         param_scaling(ii), param_shift(ii));
		}
	}


	// Time scaling

	if ( algorithm.scaling=="automatic" || algorithm.scaling!="user" )
	{
		problem.phase[i].scale.time = 1.0;

		zlower = (problem.phase[i].bounds.lower.StartTime);
		zupper = (problem.phase[i].bounds.upper.EndTime);
		if ( zlower!=-PSOPT::inf && zupper!= PSOPT::inf ) {
				if (zlower !=0.0 || zupper!=0.0)
				problem.phase[i].scale.time = 1.0/std::max( fabs(zlower), fabs(zupper));
		}
		else if (zlower==-PSOPT::inf && zupper!=PSOPT::inf && zupper!=0.0)
				problem.phase[i].scale.time = 1.0/fabs(zupper);
		else if (zupper==PSOPT::inf && zlower!=-PSOPT::inf && zlower!=0.0)
				problem.phase[i].scale.time = 1.0/fabs(zlower);

	}

  }


}

void determine_objective_scaling(MatrixXd& X,Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace )
{

  // The scaling factor for the objective function is computed such that the
  // scaled gradient at the initial guess has an Euclidean norm of 1.0.

  double nrm_g;
  MatrixXd& GF = *workspace->GFip;
  GF.resize(get_number_nlp_vars(problem, workspace), 1);


  if ( algorithm.scaling=="automatic" || algorithm.scaling!="user" )
  {

	if ( (algorithm.derivatives=="automatic") ) {
	    long      n = length(X);
	    MatrixXd&  GF = *workspace->GFip;
	    problem.scale.objective = -1.0;
	    psopt_ad::ad_record(workspace->ad_f, (int)n, 1, &X(0),
	        [&](const adouble* xin, adouble* yout){ yout[0] = ff_ad(const_cast<adouble*>(xin), workspace); });
	    std::vector<double> gtmp = psopt_ad::ad_gradient(workspace->ad_f, &X(0));
	    for(int t=0;t<(int)n;t++) GF(t) = gtmp[t];

	}

	else {
	  problem.scale.objective = -1.0;
	  ScalarGradient( ff_num, X, &GF , workspace->grw.get(), workspace );
	}

        nrm_g = (GF).norm();

   if ( nrm_g != 0.0 && nrm_g < PSOPT::inf)
	      problem.scale.objective = 1/nrm_g;
	else
	      problem.scale.objective = 1.0;

  }



}



void determine_constraint_scaling_factors(MatrixXd & X, Sol& solution, Prob& problem, Alg& algorithm, Workspace* workspace)
{
// For the differential defect constraints there are two options: The default is to use the
// same scaling factors as those used for the corresponding state. Alternatively, the user
// may specify that the scaling factors for the differential defect constraints be calculated
// based on the Jacobian (see below).
// For all other constraints, the scaling factors are computed such that the corresponding row of
// the scaled Jacobian matrix of the constraints has an Euclidean norm of 1.0.


  if ( algorithm.scaling=="automatic" || algorithm.scaling!="user" )
  {
    int i, j, l, k;

    int nvars = get_number_nlp_vars(problem, workspace);

    int ncons = get_number_nlp_constraints(problem, workspace);

//    MatrixXd& JacCol1 = *workspace->JacCol1;
    MatrixXd& xp      = *workspace->xp;
//    MatrixXd& jac_row_norm = *workspace->JacCol2;
    MatrixXd jtemp;
    
    MatrixXd JacCol1(ncons,1);
    MatrixXd jac_row_norm(ncons,1);

    workspace->use_constraint_scaling = 0;

    jac_row_norm.resize(ncons,1);


     xp = X;
//     clip_vector_given_bounds( xp, xlb, xub);

     if ( useAutomaticDifferentiation(algorithm) && algorithm.constraint_scaling=="automatic") {
        // EXTRA PARAMETER .constraint_scaling ADDED TO ALGORITHM STRUCTURE 27.11.2012.



		   jac_row_norm.setZero();
		
		   double  *x   = &xp(0);
			psopt_ad::ad_record(workspace->ad_gc, nvars, ncons, x,
				[&](const adouble* xin, adouble* yout){ gg_ad(const_cast<adouble*>(xin), yout, workspace); });
			psopt_ad::SparseTriplet J = psopt_ad::ad_sparse_jacobian(workspace->ad_gc, x, /*reuse=*/false);
			for (int t=0;t<J.nnz();t++)
			       jac_row_norm( J.row[t] ) += pow( J.val[t], 2.0);


     }
     else {

	      jac_row_norm.setZero();
	      MatrixXd& xlb = *(workspace->xlb);
		   MatrixXd& xub = *(workspace->xub);
			for(j=0;j<nvars;j++) {  // EIGEN_UPDATE
			    JacobianColumn( gg_num, xp, xlb, xub,j, &JacCol1, workspace->grw.get(), workspace);
			    jac_row_norm+= elemProduct(JacCol1, JacCol1);
			}

     } // end if-else


     jac_row_norm = (jac_row_norm.cwiseSqrt());



     // A row of the Jacobian that is zero at the initial guess carries no information
     // about that constraint's magnitude, and the earlier form of this loop -- which
     // divided by (norm + sqrt(eps)) -- gave it the largest factor the loop can produce,
     // 1/sqrt(eps) = 6.7e+07. That is not a large gradient made comparable to the others;
     // it is an accident of where the guess happens to sit. bryson_max_range asks for
     // u1^2 + u2^2 = 1 and is guessed at u = 0, the one point where that constraint's
     // gradient vanishes, so all fifty of its path rows were multiplied by 6.7e+07 and the
     // problem was handed to the NLP with an initial violation of exactly 6.7108864e+07 --
     // a violation of one, scaled. Such a row is given a factor of 1.0 instead: neutral,
     // which is the only honest reading of no information. Every other factor is clamped
     // to [1.e-7, 1.e7] symmetrically, closing a gap at exactly 1.e7 where the earlier
     // form fell through both branches and left the factor unset.

     const double sqeps    = sqrt(PSOPT_extras::GetEPS());
     const double max_fac  = 1.e7;
     const double min_fac  = 1.e-7;

     for (i=0;i<ncons;i++) // EIGEN_UPDATE
     {
            double fac;

            if ( jac_row_norm(i) <= sqeps ) {
                    // degenerate row: no scale can be inferred from it
                    fac = 1.0;
            }
            else {
                    fac = 1.0/jac_row_norm(i);
                    if ( fac > max_fac ) fac = max_fac;
                    if ( fac < min_fac ) fac = min_fac;
            }

            (*workspace->constraint_scaling)(i) = fac;
     }



     workspace->use_constraint_scaling = 1;


     if ( algorithm.defect_scaling == "jacobian-based" )
          return;


     // By default, use the state scaling factors for the differential defects. See Betts (2001).

     int offset = 0;

 
     for(i=0; i< problem.nphases; i++) {

        MatrixXd& state_scaling = (problem.phase[i].scale.states);

     	  int norder    = problem.phase[i].current_number_of_intervals;

	     int nstates   = problem.phase[i].nstates;

    


        int ncons_phase_i = get_ncons_phase_i(problem,i, workspace);


        for (k=0;k<norder+1;k++) { // EIGEN_UPDATE

  		    for( j=0;j<nstates;j++) {  // EIGEN_UPDATE

                	l = offset + (k)*nstates+j;

        		      (*workspace->constraint_scaling)(l) = state_scaling(j);

        	 }

   	  }

        offset += ncons_phase_i;

     }

  }


}
