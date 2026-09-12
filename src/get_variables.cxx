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


void get_controls(adouble* controls, adouble* xad, int iphase, int k, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
	     MatrixXd& control_scaling = problem.phase[i].scale.controls;

	     int j;

        int iphase_offset= get_iphase_offset(problem,iphase, workspace);

        // get controls

       int ncontrols = problem.phase[i].ncontrols;



        for(j=0;j<ncontrols;j++) {

           controls[j] =  xad[iphase_offset+(k)*ncontrols+j]/control_scaling(j);
        }

}

// The control of element e at its local node p (0..d) under the Nie-Kerrigan basis.
//
// Elements no longer share their end controls. The stored nodal control at a shared node is
// the LEFT element's right-hand value, so every element but the first reads its own left-hand
// value out of the block of duplicates that ir_extra_control_vars sizes. That block sits where
// the Hermite-Simpson midpoint controls would sit, which this basis does not allocate, and
// where the Gauss terminal state would sit, which cannot arise: the integrated-residual
// transcription requires Hermite-Simpson.
//
// Any reader of an element's control has to come through here. Reading get_controls(e*d)
// directly gives the neighbour's value at every interior element boundary.
void get_element_controls(adouble* controls, adouble* xad, int iphase, int e, int p, Workspace* workspace)
{
        const int i = iphase-1;
        Prob& problem = *workspace->problem;
        const int d = workspace->algorithm->ir_local_order;

        if ( p != 0 || e == 0 ) {
            get_controls(controls, xad, iphase, e*d + p, workspace);
            return;
        }

        MatrixXd& control_scaling = problem.phase[i].scale.controls;

        const int ncontrols = problem.phase[i].ncontrols;
        const int nstates   = problem.phase[i].nstates;
        const int norder    = problem.phase[i].current_number_of_intervals;
        const int nparam    = problem.phase[i].nparameters;

        // No duplicates were allocated -- a mesh this basis does not apply to, or no controls
        // at all -- so the shared node is all there is.
        if ( ir_extra_control_vars(norder, ncontrols, *workspace->algorithm) == 0 ) {
            get_controls(controls, xad, iphase, e*d, workspace);
            return;
        }

        const int iphase_offset = get_iphase_offset(problem, iphase, workspace);
        const int base = (nstates+ncontrols)*(norder+1) + nparam;

        for (int j=0; j<ncontrols; j++) {
            controls[j] = xad[iphase_offset + base + (e-1)*ncontrols + j]/control_scaling(j);
        }
}

void get_controls_bar(adouble* controls_bar, adouble* xad, int iphase, int k, Workspace* workspace)
{
   int i = iphase-1;
   Prob& problem = *workspace->problem;
	MatrixXd& control_scaling = problem.phase[i].scale.controls;

	int j;

   int iphase_offset= get_iphase_offset(problem,iphase, workspace);

	int norder    = problem.phase[i].current_number_of_intervals;
	int ncontrols = problem.phase[i].ncontrols;
	int nstates   = problem.phase[i].nstates;
        int nparam    = problem.phase[i].nparameters;

        int offset = (nstates+ncontrols)*(norder+1)+nparam;

        // The midpoint control variables are not part of the decision vector under the
        // Nie-Kerrigan local representation, where this slot belongs to that basis's
        // element-boundary controls instead -- so reading it here would not even fail
        // loudly, it would return a control belonging to a different element. Nothing calls
        // this routine there -- the residual, the cost
        // quadrature, the midpoint path rows, the estimator, integrate() and the control
        // accessor all have their own branch -- and if something ever does, it should say
        // so rather than return whatever is in the next slot.
        if ( !midpoint_control_vars(*workspace->algorithm, workspace) )
            error_message("get_controls_bar called with no midpoint control variables in the "
                          "decision vector; see midpoint_control_vars in util.cxx ");

        for(j=0;j<ncontrols;j++) {
             controls_bar[j] =  xad[iphase_offset+offset+(k)*ncontrols+j]/control_scaling(j);
        }
}


void get_final_controls(adouble* controls, adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        int k = problem.phase[i].current_number_of_intervals;
        get_controls(controls, xad, iphase, k, workspace);
}

void get_initial_controls(adouble* controls, adouble* xad, int iphase, Workspace* workspace)
{
        get_controls(controls, xad, iphase, 0, workspace);
}

void get_states(adouble* states, adouble* xad, int iphase, int k, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem            = *workspace->problem;
	MatrixXd& state_scaling   = problem.phase[i].scale.states;


	int j;

        int iphase_offset= get_iphase_offset(problem, iphase, workspace);


        int nstates = problem.phase[i].nstates;
        int ncontrols=problem.phase[i].ncontrols;
        int norder   =problem.phase[i].current_number_of_intervals;
	int offset1   = ncontrols*(norder+1);
        // get states
        for(j=0;j<nstates;j++) {
           states[j] =  xad[iphase_offset+offset1+(k)*nstates+j]/state_scaling(j);
        }

}

// Gauss: the terminal state x(t_f) is an explicit decision variable appended at the
// end of the phase's variable block (after parameters, before t0,tf), defined by the
// Gauss-quadrature constraint. This accessor reads it; centralising the offset here
// keeps the Gauss layout knowledge in one place.
void get_gauss_terminal_states(adouble* states, adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        MatrixXd& state_scaling = problem.phase[i].scale.states;
        int iphase_offset = get_iphase_offset(problem, iphase, workspace);
        int norder    = problem.phase[i].current_number_of_intervals;
        int ncontrols = problem.phase[i].ncontrols;
        int nstates   = problem.phase[i].nstates;
        int nparam    = problem.phase[i].nparameters;
        int xf_offset = (nstates+ncontrols)*(norder+1) + nparam;
        for (int j=0;j<nstates;j++)
            states[j] = xad[iphase_offset + xf_offset + j]/state_scaling(j);
}

void get_final_states(adouble* states, adouble* xad, int iphase, Workspace* workspace)
{
        Prob& problem = *workspace->problem;
        if ( workspace->algorithm->collocation_method == "Gauss" ) {
            get_gauss_terminal_states(states, xad, iphase, workspace);
            return;
        }
        int k = problem.phase[iphase-1].current_number_of_intervals;
        get_states(states, xad, iphase, k, workspace);
}

void get_initial_states(adouble* states, adouble* xad, int iphase, Workspace* workspace)
{
        get_states(states, xad, iphase, 0, workspace);
}



void get_parameters(adouble* parameters, adouble* xad, int iphase, Workspace* workspace)
{
        Prob& problem = *workspace->problem;

	int iph;

	if ( problem.multi_segment_flag || workspace->auto_linked_flag ) {
	  iph = 1;
	}
	else {
	  iph = iphase;
	}


        int i = iph-1;
        MatrixXd& param_scaling   = problem.phase[i].scale.parameters;


	int j;

	int norder    = problem.phase[i].current_number_of_intervals;
	int ncontrols = problem.phase[i].ncontrols;
	int nstates   = problem.phase[i].nstates;
        int nparam    = problem.phase[i].nparameters;
        int offset2   = (ncontrols+nstates)*(norder+1);



        int iphase_offset = get_iphase_offset(problem,iph, workspace);

        // get parameters
        for(j=0;j<nparam;j++) {
             parameters[j] =  xad[iphase_offset+offset2+j]/param_scaling(j);
        }

}

// The element boundaries in normalised [-1,1] coordinates. See ir_element_boundaries in
// psopt.h for why the caller should hoist this out of an element loop.
//
// The widths carry NO scale factor. Every other variable in the decision vector is scaled
// from its bounds, and these could be too, but they are already order 2/M by construction
// and a factor would be one more thing that has to agree between the bounds, the guess,
// the constraint row and this accessor. Unit scaling is stated here so that nobody has to
// infer it from four places agreeing.
void ir_element_boundaries(adouble* a, adouble* xad, int iphase, Workspace* workspace)
{
        const int i = iphase-1;
        Prob& problem   = *workspace->problem;
        Alg&  algorithm = *workspace->algorithm;

        const int norder = problem.phase[i].current_number_of_intervals;
        const int d      = algorithm.ir_local_order;
        if ( d < 2 || norder < d || (norder % d) != 0 ) return;
        const int M = norder/d;

        const int nflex = ir_flex_mesh_vars(norder, algorithm);

        if ( nflex == 0 ) {
            // Fixed mesh: the boundaries are stored node positions, constants to the tape.
            MatrixXd& sn = workspace->snodes[i];
            for (int e=0; e<=M; e++) a[e] = sn(e*d);
            return;
        }

        const int iphase_offset = get_iphase_offset(problem, iphase, workspace);
        const int nvars_phase_i = get_nvars_phase_i(problem, i, workspace);
        const int base          = iphase_offset + nvars_phase_i - 2 - nflex;

        a[0] = -1.0;
        for (int e=0; e<M; e++) a[e+1] = a[e] + xad[base+e];
}

// Node positions rather than element boundaries; see the declaration in psopt.h.
//
// The local abscissae are the reference element's, not the solved element's: local
// coordinates are affine-invariant, which is the whole reason a moving mesh costs so little
// here. Only the affine map a_e + lgl01(r)*h_e depends on the widths.
bool ir_node_taus(std::vector<adouble>& tau, adouble* xad, int iphase, Workspace* workspace)
{
        const int i = iphase-1;
        Prob& problem   = *workspace->problem;
        Alg&  algorithm = *workspace->algorithm;

        const int norder = problem.phase[i].current_number_of_intervals;
        const int nflex  = ir_flex_mesh_vars(norder, algorithm);

        tau.clear();
        if ( nflex == 0 ) return false;

        const int d = algorithm.ir_local_order;
        const int M = norder/d;

        std::vector<adouble> a(M+1);
        ir_element_boundaries(a.data(), xad, iphase, workspace);

        MatrixXd& lgl01 = workspace->ir_lgl01;       // d+1 reference LGL nodes on [0,1]
        tau.resize(norder+1);
        for (int e=0; e<M; e++) {
            adouble he = a[e+1] - a[e];
            for (int r=0; r<d; r++) tau[e*d + r] = a[e] + lgl01(r)*he;
        }
        tau[norder] = a[M];                          // the phase's right-hand end, exactly +1

        return true;
}

void get_times(adouble *t0, adouble *tf, adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        double   time_scaling    =  problem.phase[i].scale.time;


	     int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);

        int iphase_offset = get_iphase_offset(problem, iphase, workspace);

	     *t0  = xad[iphase_offset + nvars_phase_i-2]/time_scaling;
	     *tf  = xad[iphase_offset + nvars_phase_i-1]/time_scaling;
}

adouble get_initial_time(adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        double   time_scaling    =  problem.phase[i].scale.time;
        adouble t0;



	     int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);

        int iphase_offset = get_iphase_offset(problem,iphase, workspace);

        t0  = xad[iphase_offset + nvars_phase_i-2]/time_scaling;

        return (t0);
}

adouble get_final_time(adouble* xad, int iphase, Workspace* workspace)
{
        int i = iphase-1;
        Prob& problem = *workspace->problem;
        double   time_scaling    =  problem.phase[i].scale.time;
        adouble tf;



	     int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);

        int iphase_offset = get_iphase_offset(problem,iphase, workspace);

		  tf  = xad[iphase_offset + nvars_phase_i-1]/time_scaling;

        return (tf);
}

void get_scaled_decision_variables_and_bounds(MatrixXd& x, MatrixXd& xlb, MatrixXd& xub, Workspace* workspace)
{

    int i;

    Prob & problem = *(workspace->problem);

    adouble* xad = workspace->xad.get();


    int nvar = get_number_nlp_vars(problem, workspace);

    for(i=0;i<nvar;i++){  // EIGEN_UPDATE: index i shifted by -1
      x(i) = xad[i].value();
      xlb(i)= (*workspace->xlb)(i);
      xub(i)= (*workspace->xub)(i);
    }

}


