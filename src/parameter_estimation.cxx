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

   MatrixXd& observations = problem.phases(iphase).observations;

   MatrixXd& observation_nodes = problem.phases(iphase).observation_nodes;

   // MatrixXd& covariance = problem.phases(iphase).covariance;

   MatrixXd& residual_weights = problem.phases(iphase).residual_weights;

   adouble* interpolated_state = workspace->interp_states_pe[iphase-1];

   adouble* interpolated_control = workspace->interp_controls_pe[iphase-1];

   adouble* observed_variable = workspace->observed_variable[iphase-1];

   adouble* resid = workspace->observed_residual[iphase-1];



   // Note, the covariance matrix needs to be symmetric and positive definite.



   for(k=0; k<nsamples;k++) {   // EIGEN_UPDATE: k index shifted by -1

        time_k = observation_nodes(k);

       // Interpolate states and controls to find values at measurement instants
        for (j=0; j<problem.phases(iphase).nstates; j++) {

   	            get_interpolated_state(&interpolated_state[j], j,  iphase, time_k, xad, workspace); // EIGEN_UPDATE

   	    }

        for (j=0; j<problem.phases(iphase).ncontrols; j++) {

   	            get_interpolated_control(&interpolated_control[j], j,  iphase, time_k, xad, workspace); // EIGEN_UPDATE

   	    }

   	    // Evaluate the observations function

   	    if (problem.observation_function != NULL) {

             problem.observation_function(observed_variable, interpolated_state, interpolated_control,  parameters,  time_k,  k, xad,  iphase,workspace);

             // Compute the residual at sample k;

             for (j=0;j<nobserved;j++) {
                  resid[j] = residual_weights(j,k)*( observed_variable[j] - observations(j, k) ); // EIGEN_UPDATE
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





void auto_split_observations(Prob& problem, MatrixXd& observation_nodes, MatrixXd& observations )
{
    int i, kstart, kend;
    // int nobs = problem.phases(1).nobserved;

    MatrixXd& tobs = observation_nodes;
    MatrixXd& xobs = observations;

    kstart = 0;
    kend   = 0;

    for (i=1; i<=problem.nphases; i++ ) {
            if (i>1) kstart = kend; else kstart = 0;
            kend   = kstart + problem.phases(i).nsamples-1;

       problem.phases(i).observation_nodes      = tobs.block(0,kstart,1,problem.phases(i).nsamples); 
       problem.phases(i).observations           = xobs.block(0,kstart,xobs.rows(), problem.phases(i).nsamples ); 
	    problem.phases(i).residual_weights       = problem.phases(1).residual_weights;
	    problem.phases(i).covariance             = problem.phases(1).covariance;
	    problem.phases(i).regularization_factor  = problem.phases(1).regularization_factor;
    }


}




void load_parameter_estimation_data(Prob& problem, int iphase, const char* filename)
{
    // This function reads data for parameter estimation problems for a given phase "iphase" from a file
    // whose columns contain the following information:
    // TIME   VAR1  WEIGHT1  VAR2   WEIGHT2  VAR3  WEIGHT3
    // where TIME represents the sampling instants, VAR1 is the observed variable # 1, WEIGHT1 is the weight of observed
    // variable 1, and so on.
    // The dimensions of the matrix in the data file is  problem.phases(iphase).nsamples x 2*problem.phases(iphase).nobserved + 1
    //

     MatrixXd data, tm, ym, w;

     int nsamples = problem.phases(iphase).nsamples;

     int nobserved = problem.phases(iphase).nobserved;

     data.resize(nsamples,nobserved*2 + 1);

     data = load_data(filename, nsamples, 1+nobserved*2);


     tm = data.block(0,0, data.rows(),1);

     ym.resize(nsamples, nobserved);
     w.resize(nsamples,nobserved);

     for (int i=0;i<nobserved;i++) {   // EIGEN_UPDATE index i shifted by -1

            ym.block(0,i,ym.rows(),1) =  data.block(0, 1+i*2, ym.rows(), 1);

            w.block(0,i,w.rows(),1)   =  data.block(0, 2+i*2, w.rows(), 1);
     }

     problem.phases(iphase).observation_nodes   = tm.transpose(); // tra(tm);
     problem.phases(iphase).observations        = ym.transpose(); // tra(ym);
     problem.phases(iphase).residual_weights    = w.transpose();  // tra(w);


}



void compute_residual_vector_in_phase(MatrixXd& residual_vector, adouble* xad, int iphase, Workspace* workspace)
{


   adouble time_k;

   Prob & problem = *(workspace->problem);

   int nsamples = problem.phases(iphase).nsamples;

   int nobserved = problem.phases(iphase).nobserved;

   int k, j;

   MatrixXd& observations = problem.phases(iphase).observations;

   MatrixXd& observation_nodes = problem.phases(iphase).observation_nodes;

   MatrixXd& residual_weights = problem.phases(iphase).residual_weights;

   adouble* interpolated_state = workspace->interp_states_pe[iphase-1];

   adouble* interpolated_control = workspace->interp_controls_pe[iphase-1];

   adouble* observed_variable = workspace->observed_variable[iphase-1];

   adouble* resid = workspace->observed_residual[iphase-1];

   residual_vector.resize( nsamples*nobserved, 1);

   adouble* parameters;

   parameters    = workspace->parameters[iphase-1];

   get_parameters(parameters, xad, iphase, workspace);

   for(k=0; k<nsamples;k++) { // EIGEN_UPDATE: k index shifted by -1

        time_k = observation_nodes(k);

       // Interpolate states and controls to find values at measurement instants
        for (j=0; j<problem.phases(iphase).nstates; j++) {

   	            get_interpolated_state(&interpolated_state[j], j,  iphase, time_k, xad, workspace); // EIGEN_UPDATE

   	    }

        for (j=0; j<problem.phases(iphase).ncontrols; j++) {

   	            get_interpolated_control(&interpolated_control[j], j,  iphase, time_k, xad, workspace); // EIGEN_UPDATE

   	    }

   	    // Evaluate the observations function

   	    if (problem.observation_function != NULL) {

             problem.observation_function(observed_variable, interpolated_state, interpolated_control,  parameters,  time_k,  k, xad,  iphase, workspace);

             // Compute the residual at sample k;

             for (j=0;j<nobserved;j++) {
                resid[j] = residual_weights(j,k)*( observed_variable[j] - observations(j, k) ); // EIGEN_UPDATE

		residual_vector( (k)*nobserved+j  ) = resid[j].value(); // EIGEN_UPDATE

             }


   	}


   }




}



void rr_num(MatrixXd& X, MatrixXd* residual_vector, Workspace* workspace)
{



   adouble time_k;

   int index, iphase, dindex, j;

   Prob & problem = *(workspace->problem);


   MatrixXd residual_vector_in_phase;

   adouble* xad = workspace->xad;

   int nvars = get_number_nlp_vars(problem, workspace);

   for(j=0;j<nvars;j++) {
      xad[j] = X(j); // EIGEN_UPDATE
   }

   index = 0;

   for(iphase=1; iphase<=problem.nphases;iphase++)
   {
       dindex = problem.phases(iphase).nobserved*problem.phases(iphase).nsamples;
       index = index+ dindex;
   }

   residual_vector->resize( index, 1);

   index = 0;

   for(iphase=1; iphase<=problem.nphases;iphase++)
   {

       dindex = problem.phases(iphase).nobserved*problem.phases(iphase).nsamples;
       residual_vector_in_phase.resize( dindex, 1);
       compute_residual_vector_in_phase( residual_vector_in_phase, xad, iphase, workspace);
       (*residual_vector).block(index, 0, dindex, 1) = residual_vector_in_phase;
       index = index+ dindex;
   }



}




void extract_parameter_covariance(MatrixXd& Cp, MatrixXd& C, Workspace* workspace)
{
     int i, j, ii;
     Prob & problem = *(workspace->problem);
     RowVectorXi Ip;
     int pcount = 0;

     for(i=0;i< problem.nphases; i++)
     {
         pcount+= problem.phase[i].nparameters;
     }

     Ip.resize(pcount);

     pcount=0;

     for(i=0;i< problem.nphases; i++)
     {
	     int norder    = problem.phase[i].current_number_of_intervals;
     	  int nstates   = problem.phase[i].nstates;
        int ncontrols = problem.phase[i].ncontrols;
        int nparam    = problem.phase[i].nparameters;

        int iphase_offset = get_iphase_offset(problem,i+1, workspace); // CAREFUL HERE: CHECK THIS FUNCTION'S SECOND PARAM


        int offset2 = (norder+1)*(ncontrols+nstates);

	      for (ii=0;ii<nparam;ii++) { // EIGEN_UPDATE: index ii shifted by -1
                        j = iphase_offset+offset2+ii;

                        Ip(pcount) = j;
                        pcount++;
	      }
     }

     Cp.resize(Ip.size(), Ip.size());
     for (i=0; i< Ip.size(); i++) {
         for (j=0; j< Ip.size(); j++) {
              Cp(i,j) = C(Ip(i), Ip(j));
         }
      }

}





bool compute_parameter_statistics(MatrixXd& Cp, MatrixXd& p, MatrixXd& plow, MatrixXd& phigh, MatrixXd& r, Workspace* workspace)
{
      MatrixXd Jr, Jc;

      MatrixXd X, XL, XU;

      Prob & problem = *(workspace->problem);

      adouble* xad = workspace->xad;

      int iphase, i, j;

      adouble* parameters;

      int pcount = 0;

      double alpha;

      char* msg = new char[100];

      sprintf(msg,"\n>>> Performing statistical analysis of estimated parameters...");

	  psopt_print(workspace,msg);


      for(i=0;i< problem.nphases; i++)
      {
         pcount+= problem.phase[i].nparameters;
      }

      int total_number_of_parameters = pcount;

      MatrixXd parameters_full(pcount,1);

      pcount = 0;

      for(i=0;i< problem.nphases; i++)
      {

	       int iphase = i+1;
	       int npar = problem.phase[i].nparameters;

	       parameters    = workspace->parameters[iphase-1];

	       get_parameters( parameters, xad, iphase, workspace );

	       for(j=0; j< npar; j++) {
	            parameters_full(pcount+j) = parameters[j].value();  // EIGEN_UPDATE
	       }

          pcount+= problem.phase[i].nparameters;
      }

      int nvar = get_number_nlp_vars(problem, workspace);

      X .resize(nvar,1);
      XL.resize(nvar,1);
      XU.resize(nvar,1);


      get_scaled_decision_variables_and_bounds(X, XL, XU, workspace);

      compute_jacobian_of_residual_vector_with_respect_to_variables(Jr, X, XL, XU, workspace);

      compute_jacobian_of_constraints_with_respect_to_variables(Jc, X, XL, XU, workspace);

	  rr_num(X,&r, workspace);

	  p = parameters_full;

     MatrixXd JcT(Jc.cols(),Jc.rows());

      for(i=0;i<Jc.rows();i++) {      // EIGEN_UPDATE
          for(j=0;j<Jc.cols();j++) {
             JcT(j,i) = Jc(i,j);
          }
      }

      long M= JcT.rows();
      long N= JcT.cols();


      if ((M-N)<=0) return false;



       MatrixXd Q = JcT.fullPivHouseholderQr().matrixQ();
      
      MatrixXd I2( N+(M-N) , (M-N) );
      
      I2 << zeros(N,M-N),
            eye(M-N);
            
       MatrixXd Z = Q*I2;            
            

      // Calculation of covariance matrix w.r.t all decision variables.
      // See the paper:
      // Kostina et al (2003) "Computation of covariance matrices for constrained parameter estimation
      // problems using LSQR".



      MatrixXd& J1 = Jr;
      MatrixXd CC;

      MatrixXd F;

      F = Z.transpose()*J1.transpose()*J1*Z;

      F.completeOrthogonalDecomposition().pseudoInverse();

        CC = Z*F*Z.transpose();

      extract_parameter_covariance(Cp, CC, workspace);

      // Compute confidence intervals

      alpha = 0.95; // 95% Confidence level

      pcount=0;

      plow.resize( total_number_of_parameters, 1);
      phigh.resize(total_number_of_parameters,1);

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

       j= pcount+i; // EIGEN_UPDATE

	    plow(j)   = parameters_full(j) - tt*sqrt( Cp(j,j) );
	    phigh(j)  = parameters_full(j) + tt*sqrt( Cp(j,j) );

	  }

	  pcount += problem.phases(iphase).nparameters;


      }



	  return true;


}

