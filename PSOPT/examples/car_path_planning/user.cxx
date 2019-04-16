

#include "psopt.h"

const double Wheel_Base_Len = 2.67;

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase, Workspace* workspace)
{
   return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters,
                     adouble& time, adouble* xad, int iphase, Workspace* workspace)
{
   const double weight_acc = pow(0.2, 2.0);
   const double weight_omg = pow(0.05, 2.0);
   adouble x = states[ CINDEX(1) ];
   adouble y = states[ CINDEX(2) ];
   adouble v = states[ CINDEX(3) ];
   adouble theta = states[ CINDEX(4) ];
   adouble phi = states[ CINDEX(5) ];

   adouble acc = controls[ CINDEX(1) ];
   adouble omg = controls[ CINDEX(2) ];
   return  0.5 * (weight_acc * acc * acc + weight_omg * omg * omg);
}


//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states,
         adouble* controls, adouble* parameters, adouble& time,
         adouble* xad, int iphase, Workspace* workspace)
{
   adouble x = states[ CINDEX(1) ];
   adouble y = states[ CINDEX(2) ];
   adouble v = states[ CINDEX(3) ];
   adouble theta = states[ CINDEX(4) ];
   adouble phi = states[ CINDEX(5) ];

   adouble acc = controls[ CINDEX(1) ];
   adouble omg = controls[ CINDEX(2) ];

   derivatives[ CINDEX(1) ] = v * cos(theta);
   derivatives[ CINDEX(2) ] = v * sin(theta);
   derivatives[ CINDEX(3) ] = acc;
   derivatives[ CINDEX(4) ] = v * tan(phi) / Wheel_Base_Len;
   derivatives[ CINDEX(5) ] = omg;
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase, Workspace* workspace)

{
   adouble x_0 = initial_states[ CINDEX(1) ];
   adouble y_0 = initial_states[ CINDEX(2) ];
   adouble v_0 = initial_states[ CINDEX(3) ];
   adouble theta_0 = initial_states[ CINDEX(4) ];
   adouble phi_0 = initial_states[ CINDEX(5) ];

   adouble x_f = final_states[ CINDEX(1) ];
   adouble y_f = final_states[ CINDEX(2) ];
   adouble v_f = final_states[ CINDEX(3) ];
   adouble theta_f = final_states[ CINDEX(4) ];
   adouble phi_f = final_states[ CINDEX(5) ];

   e[ CINDEX(1) ] = x_0;
   e[ CINDEX(2) ] = y_0;
   e[ CINDEX(3) ] = v_0;
   e[ CINDEX(4) ] = theta_0;
   e[ CINDEX(5) ] = phi_0;
   e[ CINDEX(6) ] = x_f;
   e[ CINDEX(7) ] = y_f;
   e[ CINDEX(8) ] = v_f;
   e[ CINDEX(9) ] = theta_f;
   e[ CINDEX(10) ] = phi_f;

}



///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages( adouble* linkages, adouble* xad, Workspace* workspace)
{
  // No linkages as this is a single phase problem
}



////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ////CINDEX///////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg  algorithm;
    Sol  solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name        		= "Car-like robot path planning problem";
    problem.outfilename                 = "path_planning.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases   			= 1;
    problem.nlinkages                   = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup  ////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates   		= 5;
    problem.phases(1).ncontrols 		= 2;
    problem.phases(1).nevents   		= 10;
    problem.phases(1).npath     		= 0;
    problem.phases(1).nodes                     = "[30]";

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////
    double states_upper_bound[5] = {20.0, 20.0, 10.0 / 3.6, 2.0 * pi, 25.0 / 180.0 * pi};
    double states_lower_bound[5];
    for (size_t i = 0; i < problem.phases(1).nstates; ++i)
    {
        states_lower_bound[i] = -states_upper_bound[i];
        problem.phases(1).bounds.lower.states(i + 1) = states_lower_bound[i];
        problem.phases(1).bounds.upper.states(i + 1) = states_upper_bound[i];
    }
    double controls_upper_bound[2] = {3.0, 25.0 / 180.0 * pi};
    double controls_lower_bound[2];
    for (size_t i = 0; i < problem.phases(1).ncontrols; ++i)
    {
        controls_lower_bound[i] = -controls_upper_bound[i];
        problem.phases(1).bounds.lower.controls(i + 1) = controls_lower_bound[i];
        problem.phases(1).bounds.upper.controls(i + 1) = controls_upper_bound[i];
    }
    const double x_0 = 0.0;
    const double y_0 = 0.0;
    const double v_0 = 0.0;
    const double theta_0 = 0.0;
    const double phi_0 = 0.0;

    const double x_f = 0.0;
    const double y_f = 0.0;
    const double v_f = 0.0;
    const double theta_f = pi;
    const double phi_f = 0.0;
    double events_upper_bound[10] = {x_0, y_0, v_0, theta_0, phi_0, 
        x_f, y_f, v_f, theta_f, phi_f};
    double events_lower_bound[10] = {x_0, y_0, v_0, theta_0, phi_0,
        x_f, y_f, v_f, theta_f, phi_f};
    for (size_t i = 0; i < problem.phases(1).nevents; ++i)
    {
        problem.phases(1).bounds.lower.events(i + 1) = events_lower_bound[i];
        problem.phases(1).bounds.upper.events(i + 1) = events_upper_bound[i];
    }

    problem.phases(1).bounds.lower.StartTime    = 0.0;
    problem.phases(1).bounds.upper.StartTime    = 0.0;

    problem.phases(1).bounds.lower.EndTime      = 0.0;
    problem.phases(1).bounds.upper.EndTime      = 60.0;



////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////


    problem.integrand_cost 	= &integrand_cost;
    problem.endpoint_cost 	= &endpoint_cost;
    problem.dae             	= &dae;
    problem.events 		= &events;
    problem.linkages		= &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes    			= problem.phases(1).nodes(1);
    int ncontrols                       = problem.phases(1).ncontrols;
    int nstates                         = problem.phases(1).nstates;

    DMatrix x_guess    =  zeros(nstates,nnodes);

    double xf = 5.0;
    double x0 = 0.0;
    for (int index = 1; index <= nnodes; ++index)
    {
        x_guess(1, index) = x0 + (xf - x0) / (nnodes - 1) * (index - 1);
    }
    // x_guess(1,colon()) = x0*ones(1,nnodes);
    // x_guess(2,colon()) = y0*ones(1,nnodes);
    // x_guess(3,colon()) = theta0*ones(1,nnodes);

    problem.phases(1).guess.controls       = zeros(ncontrols,nnodes);
    problem.phases(1).guess.states         = x_guess;
    problem.phases(1).guess.time           = linspace(0.0, xf, nnodes);


////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


    algorithm.nlp_iter_max                = 1000;
    algorithm.nlp_tolerance               = 1.e-3;
    algorithm.nlp_method                  = "IPOPT";
    algorithm.scaling                     = "automatic";
    algorithm.derivatives                 = "automatic";
    // algorithm.mesh_refinement             = "automatic";
    // algorithm.collocation_method = "Legendre";
    // algorithm.collocation_method = "Chebyshev";
    // algorithm.collocation_method = "trapezoidal";
    algorithm.collocation_method = "Hermite-Simpson";
 
//    algorithm.defect_scaling = "jacobian-based";
    algorithm.ode_tolerance               = 1.e-3;



////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

    psopt(solution, problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////


    x      = solution.get_states_in_phase(1);
    u      = solution.get_controls_in_phase(1);
    t      = solution.get_time_in_phase(1);
    lambda = solution.get_dual_costates_in_phase(1);
    H      = solution.get_dual_hamiltonian_in_phase(1);


////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");
    lambda.Save("lambda.dat");
    H.Save("H.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    // plot(t,x,problem.name+": states", "time (s)", "states","x y v");

    // plot(t,u,problem.name+": controls","time (s)", "controls", "u_1 u_2");

    plot(t,x,problem.name+": states", "time (s)", "states","x y v theta phi",
                             "pdf", "car_path_planning_states.pdf");
    plot(x(1, colon()),x(2, colon()),problem.name+": path", "x (m)", "y(m)","path",
                             "pdf", "car_path_planning_path.pdf");

    plot(t,u,problem.name+": controls","time (s)", "controls", "acc omg",
                             "pdf", "car_path_planning_controls.pdf");
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////
