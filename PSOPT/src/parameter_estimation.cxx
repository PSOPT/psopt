/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2015 Victor M. Becerra

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
           University of Reading
           School of Systems Engineering
           P.O. Box 225, Reading RG6 6AY
           United Kingdom
           e-mail: vmbecerra99@gmail.com

**********************************************************************************************/


#include "psopt.h"

adouble endpoint_cost_for_parameter_estimation(adouble* initial_states, adouble* final_states, adouble* parameters,adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace)
{
   // This is the end point cost function for parameter estimation problems.

   adouble retval = 0.0;

   adouble time_k;

   Prob & problem = *(workspace->problem);

   Alg & algorithm = *(workspace->algorithm);

   int nsamples = problem.phases(iphase).nsamples;

   int nobserved = problem.phases(iphase).nobserved;

   int k, j;

   DMatrix& observations = problem.phases(iphase).observations;

   DMatrix& observation_nodes = problem.phases(iphase).observation_nodes;

   DMatrix& covariance = problem.phases(iphase).covariance;

   DMatrix& residual_weights = problem.phases(iphase).residual_weights;

   adouble* interpolated_state = workspace->interp_states_pe[iphase-1];

   adouble* interpolated_control = workspace->interp_controls_pe[iphase-1];

   adouble* observed_variable = workspace->observed_variable[iphase-1];

   adouble* resid = workspace->observed_residual[iphase-1];

   adouble* lam_resid = workspace->lam_resid[iphase-1];


//   DMatrix  Lambda; //  = Sqrt( inv( covariance ) );

   // Note, the covariance matrix needs to be symmetric and positive definite.

//   DMatrix inv_cov= inv(covariance);

//   DMatrix V( inv_cov.GetNoRows(), 2*inv_cov.GetNoCols());

//   DMatrix E = eig( inv_cov, & V );

//   E = E.sub_matrix(1,E.GetNoRows(),1,1);

//   V = V.sub_matrix(1,V.GetNoRows(),1,V.GetNoCols()/2);
//   DMatrix D = diag(E);


//   Lambda = V*Sqrt(D)*tra(V);  // Lambda is the square root of the covariance matrix.


   for(k=1; k<=nsamples;k++) {

        time_k = observation_nodes(k);

       // Interpolate states and controls to find values at measurement instants
        for (j=0; j<problem.phases(iphase).nstates; j++) {

   	            get_interpolated_state(&interpolated_state[j], j+1,  iphase, time_k, xad, workspace);

   	    }

        for (j=0; j<problem.phases(iphase).ncontrols; j++) {

   	            get_interpolated_control(&interpolated_control[j], j+1,  iphase, time_k, xad, workspace);

   	    }

   	    // Evaluate the observations function

   	    if (problem.observation_function != NULL) {

             problem.observation_function(observed_variable, interpolated_state, interpolated_control,  parameters,  time_k,  k, xad,  iphase,workspace);

             // Compute the residual at sample k;

             for (j=0;j<nobserved;j++) {
                  resid[j] = residual_weights(j+1,k)*( observed_variable[j] - observations(j+1, k) );
             }

             if (algorithm.parameter_estimation_norm==2) {

                 retval += dot(resid, resid, nobserved );

   	         }

             else if (algorithm.parameter_estimation_norm==1)  {

                         for (j=0;j<nobserved;j++) {
                            retval += fabs(resid[j]);
                         }

   	         }

             else if ( algorithm.parameter_estimation_norm>2) {
                         adouble dphi = 0.0;
                         double    p    = algorithm.parameter_estimation_norm;
                         for (j=0;j<nobserved;j++) {
                            dphi += pow(resid[j],p);
                         }
                         retval += dphi;

             }


   	    }




   }

   if (algorithm.parameter_estimation_norm>2) {
        double p = algorithm.parameter_estimation_norm;
        retval = pow( retval, 1.0/p);
   }

   // Add regularization term.

   adouble param_norm2 = dot(parameters, parameters, problem.phases(iphase).nparameters );

   retval += problem.phases(iphase).regularization_factor*param_norm2;

   return (retval);

}





void auto_split_observations(Prob& problem, DMatrix& observation_nodes, DMatrix& observations )
{
    int i, kstart, kend;
    int nobs = problem.phases(1).nobserved;

    DMatrix& tobs = observation_nodes;
    DMatrix& xobs = observations;

    kstart = 0;
    kend   = 0;

    for (i=1; i<=problem.nphases; i++ ) {
            if (i>1) kstart = kend; else kstart = 1;
            kend   = kstart + problem.phases(i).nsamples-1;
	    problem.phases(i).observation_nodes      = tobs(1, colon(kstart,kend) );
	    problem.phases(i).observations           = xobs(colon(), colon(kstart,kend) );
	    problem.phases(i).residual_weights       = problem.phases(1).residual_weights;
	    problem.phases(i).covariance             = problem.phases(1).covariance;
	    problem.phases(i).regularization_factor  = problem.phases(1).regularization_factor;
    }


}




void load_parameter_estimation_data(Prob& problem, int iphase, char* filename)
{
    // This function reads data for parameter estimation problems for a given phase "iphase" from a file
    // whose columns contain the following information:
    // TIME   VAR1  WEIGHT1  VAR2   WEIGHT2  VAR3  WEIGHT3
    // where TIME represents the sampling instants, VAR1 is the observed variable # 1, WEIGHT1 is the weight of observed
    // variable 1, and so on.
    // The dimensions of the matrix in the data file is  problem.phases(iphase).nsamples x 2*problem.phases(iphase).nobserved + 1
    //

     DMatrix data, tm, ym, w;

     int nsamples = problem.phases(iphase).nsamples;

     int nobserved = problem.phases(iphase).nobserved;

     data.Resize(nsamples,nobserved*2 + 1);

     data.Load(filename);

     tm = data(colon(),1);

     ym.Resize(nsamples, nobserved);
     w.Resize(nsamples,nobserved);

     for (int i=1;i<=nobserved;i++) {
          ym(colon(), i) = data(colon(), 2+(i-1)*2);
          w (colon(), i) = data(colon(), 3+(i-1)*2);
     }

     problem.phases(iphase).observation_nodes   = tra(tm);
     problem.phases(iphase).observations        = tra(ym);
     problem.phases(iphase).residual_weights    = tra(w);


}



void compute_residual_vector_in_phase(DMatrix& residual_vector, adouble* xad, int iphase, Workspace* workspace)
{


   adouble time_k;

   Prob & problem = *(workspace->problem);

   int nsamples = problem.phases(iphase).nsamples;

   int nobserved = problem.phases(iphase).nobserved;

   int k, j;

   DMatrix& observations = problem.phases(iphase).observations;

   DMatrix& observation_nodes = problem.phases(iphase).observation_nodes;

   DMatrix& residual_weights = problem.phases(iphase).residual_weights;

   adouble* interpolated_state = workspace->interp_states_pe[iphase-1];

   adouble* interpolated_control = workspace->interp_controls_pe[iphase-1];

   adouble* observed_variable = workspace->observed_variable[iphase-1];

   adouble* resid = workspace->observed_residual[iphase-1];

   residual_vector.Resize( nsamples*nobserved, 1);

   adouble* parameters;

   parameters    = workspace->parameters[iphase-1];

   get_parameters(parameters, xad, iphase, workspace);

   for(k=1; k<=nsamples;k++) {

        time_k = observation_nodes(k);

       // Interpolate states and controls to find values at measurement instants
        for (j=0; j<problem.phases(iphase).nstates; j++) {

   	            get_interpolated_state(&interpolated_state[j], j+1,  iphase, time_k, xad, workspace);

   	    }

        for (j=0; j<problem.phases(iphase).ncontrols; j++) {

   	            get_interpolated_control(&interpolated_control[j], j+1,  iphase, time_k, xad, workspace);

   	    }

   	    // Evaluate the observations function

   	    if (problem.observation_function != NULL) {

             problem.observation_function(observed_variable, interpolated_state, interpolated_control,  parameters,  time_k,  k, xad,  iphase, workspace);

             // Compute the residual at sample k;

             for (j=0;j<nobserved;j++) {
                resid[j] = residual_weights(j+1,k)*( observed_variable[j] - observations(j+1, k) );

		residual_vector( (k-1)*nobserved+j+1  ) = resid[j].value();

             }


   	}


   }




}



void rr_num(DMatrix& X, DMatrix* residual_vector, Workspace* workspace)
{



   adouble time_k;

   int index, iphase, dindex, j,  k;

   Prob & problem = *(workspace->problem);

   adouble* parameters;

   DMatrix residual_vector_in_phase;

   adouble* xad = workspace->xad;

   int nvars = get_number_nlp_vars(problem, workspace);

   for(j=0;j<nvars;j++) {
      xad[j] = X(j+1);
   }

   index = 0;

   for(iphase=1; iphase<=problem.nphases;iphase++)
   {
       dindex = problem.phases(iphase).nobserved*problem.phases(iphase).nsamples;
       index = index+ dindex;
   }

   residual_vector->Resize( index, 1);

   index = 0;

   for(iphase=1; iphase<=problem.nphases;iphase++)
   {
       parameters    = workspace->parameters[iphase-1];
       dindex = problem.phases(iphase).nobserved*problem.phases(iphase).nsamples;
       residual_vector_in_phase.Resize( dindex, 1);
       compute_residual_vector_in_phase( residual_vector_in_phase, xad, iphase, workspace);
       (*residual_vector)( colon(index+1, index + dindex), 1) = residual_vector_in_phase;
       index = index+ dindex;
   }



}






void get_scaled_decision_variables_and_bounds(DMatrix& x, DMatrix xlb, DMatrix xub, Workspace* workspace)
{

    int i;

    Prob & problem = *(workspace->problem);

    adouble* xad = workspace->xad;


    int nvar = get_number_nlp_vars(problem, workspace);

    for(i=1;i<=nvar;i++){
      x(i) = xad[i-1].value();
      xlb(i)= (*workspace->xlb)(i);
      xub(i)= (*workspace->xub)(i);
    }

}



void extract_parameter_covariance(DMatrix& Cp, DMatrix& C, Workspace* workspace)
{
     int i, j, ii;
     Prob & problem = *(workspace->problem);
     DMatrix Ip;
     int pcount = 0;

     for(i=0;i< problem.nphases; i++)
     {
         pcount+= problem.phase[i].nparameters;
     }

     Ip.Resize(pcount,1);

     pcount=0;

     for(i=0;i< problem.nphases; i++)
     {
	int norder    = problem.phase[i].current_number_of_intervals;
	int nstates   = problem.phase[i].nstates;
        int ncontrols = problem.phase[i].ncontrols;
        int nparam    = problem.phase[i].nparameters;
	int offset;
        int iphase_offset = get_iphase_offset(problem,i+1, workspace);
	int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);
        int offset1 = (norder+1)*ncontrols;
        int offset2 = (norder+1)*(ncontrols+nstates);

	for (ii=1;ii<=nparam;ii++) {
                        j = iphase_offset+offset2+ii;
			pcount++;
                        Ip(pcount,1)=(double) j;
	}
     }

     Cp = C( Ip, Ip );

}





bool compute_parameter_statistics(DMatrix& Cp, DMatrix& p, DMatrix& plow, DMatrix& phigh, DMatrix& r, Workspace* workspace)
{
      DMatrix Jr, Jc;

      DMatrix X, XL, XU;

      Prob & problem = *(workspace->problem);

      adouble* xad = workspace->xad;

      int iphase, i, j;

      adouble* parameters;

      int pcount = 0;

      double alpha;

      char* msg = new char[60];

      sprintf(msg,"\n>>> Performing statistical analysis of estimated parameters...");

	  psopt_print(workspace,msg);


      for(i=0;i< problem.nphases; i++)
      {
         pcount+= problem.phase[i].nparameters;
      }

      int total_number_of_parameters = pcount;

      DMatrix parameters_full(pcount,1);

      pcount = 0;

      for(i=0;i< problem.nphases; i++)
      {

	 int iphase = i+1;
	 int npar = problem.phase[i].nparameters;

	 parameters    = workspace->parameters[iphase-1];

	 get_parameters( parameters, xad, iphase, workspace );

	 for(j=0; j< npar; j++) {
	   parameters_full(pcount+j+1) = parameters[j].value();
	 }

         pcount+= problem.phase[i].nparameters;
      }

      int nvar = get_number_nlp_vars(problem, workspace);

      X .Resize(nvar,1);
      XL.Resize(nvar,1);
      XU.Resize(nvar,1);


      get_scaled_decision_variables_and_bounds(X, XL, XU, workspace);

      compute_jacobian_of_residual_vector_with_respect_to_variables(Jr, X, XL, XU, workspace);

      compute_jacobian_of_constraints_with_respect_to_variables(Jc, X, XL, XU, workspace);

	  rr_num(X,&r, workspace);

	  p = parameters_full;

      DMatrix JcT(Jc.GetNoCols(),Jc.GetNoRows());

      for(i=1;i<=Jc.GetNoRows();i++) {
          for(j=1;j<=Jc.GetNoCols();j++) {
             JcT(j,i) = Jc(i,j);
          }
      }

      integer M= JcT.GetNoRows();
      integer N= JcT.GetNoCols();
      double* A = JcT.GetPr();
      integer LDA = M;
      double* TAU = new double[MIN(M,N)];
      integer LWORK = 6*MAX(1,N);
      double* WORK= new double[MAX(1,LWORK)];
      integer INFO;

      if ((M-N)<=0) return false;

      dgeqrf_( &M, &N, A, &LDA, TAU, WORK, &LWORK, &INFO );

      DMatrix I2 = ( zeros(N,M-N) && eye(M-N) );

      double *C = I2.GetPr();
      integer LDC = MAX(1,M);

      char SIDE = 'L';
      char TRANS = 'N';

      integer N2 = M-N;
      integer K  = MIN(M,N);



      dormqr_( &SIDE, &TRANS, &M, &N2, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, &INFO, 1, 1);

      // Calculation of covariance matrix w.r.t all decision variables.
      // See the paper:
      // Kostina et al (2003) "Computation of covariance matrices for constrained parameter estimation
      // problems using LSQR".

      DMatrix& Z = I2;

      DMatrix& J1 = Jr;
      DMatrix CC;

      DMatrix F;

      F = (TProductT(Z,J1)*(J1*Z));

      F = pinv(F);
//      F = inv(F);

      CC = Z*F*tra(Z);

      extract_parameter_covariance(Cp, CC, workspace);

      // Compute confidence intervals

      alpha = 0.95; // 95% Confidence level

      pcount=0;

      plow.Resize( total_number_of_parameters, 1);
      phigh.Resize(total_number_of_parameters,1);

      int NN=0;

      for (iphase=1; iphase<=problem.nphases; iphase++)
      {
         NN += problem.phases(iphase).nsamples*problem.phases(iphase).nobserved;
      }

	  int np = total_number_of_parameters;

      for (iphase=1; iphase<=problem.nphases; iphase++)
      {

      double tt;

      tt = inverse_twotailed_t_cdf(alpha, (NN-np) );

      // Compute the limits of the confidence intervals

	  for(i=0;i< problem.phases(iphase).nparameters;i++) {

	    j= pcount+i+1;

	    plow(j,1)   = parameters_full(j) - tt*sqrt( Cp(j,j) );
	    phigh(j,1)  = parameters_full(j) + tt*sqrt( Cp(j,j) );

	  }

	  pcount += problem.phases(iphase).nparameters;


      }



	  return true;


}

