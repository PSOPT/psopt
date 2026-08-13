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
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/


#include "psopt.h"

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;

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

   adouble* interpolated_state = workspace->interp_states_pe[iphase-1].get();

   adouble* interpolated_control = workspace->interp_controls_pe[iphase-1].get();

   adouble* observed_variable = workspace->observed_variable[iphase-1].get();

   adouble* resid = workspace->observed_residual[iphase-1].get();



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

   adouble* interpolated_state = workspace->interp_states_pe[iphase-1].get();

   adouble* interpolated_control = workspace->interp_controls_pe[iphase-1].get();

   adouble* observed_variable = workspace->observed_variable[iphase-1].get();

   adouble* resid = workspace->observed_residual[iphase-1].get();

   residual_vector.resize( nsamples*nobserved, 1);

   adouble* parameters;

   parameters    = workspace->parameters[iphase-1].get();

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

   adouble* xad = workspace->xad.get();

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





// Scale factors of the static parameters, in the order in which
// extract_parameter_covariance lays out the rows and columns of Cp: phase by phase, and
// within a phase in declaration order. The scaled decision variable is s*p (see
// get_parameters), so a covariance computed in scaled variables is brought to physical
// units by dividing entry (i,j) by s_i*s_j.
static void get_parameter_scale_factors(MatrixXd& s, Workspace* workspace)
{
     Prob& problem = *(workspace->problem);
     int n = 0;
     for (int i = 0; i < problem.nphases; i++) n += problem.phase[i].nparameters;
     s = ones( n>0 ? n : 1, 1);
     int c = 0;
     for (int i = 0; i < problem.nphases; i++) {
         MatrixXd& ps = problem.phase[i].scale.parameters;
         for (int j = 0; j < problem.phase[i].nparameters; j++) {
             double f = ( ps.size() > j ) ? ps(j) : 1.0;
             s(c++) = ( f != 0.0 ) ? f : 1.0;
         }
     }
}


void store_parameter_statistics(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
     if ( problem.observation_function == NULL )          return;
     if ( algorithm.parameter_statistics != "yes" )       return;

     int nparam = get_total_number_of_parameters(problem);

     MatrixXd Cp(nparam,nparam), plow(nparam,1), phigh(nparam,1), p(nparam,1), r;

     bool ok = compute_parameter_statistics(Cp, p, plow, phigh, r, workspace);

     solution.parameter_statistics_ok = ok;
     solution.observation_residuals   = r;

     if ( ok ) {
         solution.parameter_covariance      = Cp;
         solution.parameter_confidence_low  = plow;
         solution.parameter_confidence_high = phigh;
         solution.sigma_hat                 = workspace->pe_sigma_hat;
         solution.parameter_dof             = workspace->pe_dof;
         solution.n_observations            = workspace->pe_nobs;
     }
}


bool compute_parameter_statistics(MatrixXd& Cp, MatrixXd& p, MatrixXd& plow, MatrixXd& phigh, MatrixXd& r, Workspace* workspace)
{
      MatrixXd Jr, Jc;

      MatrixXd X, XL, XU;

      Prob & problem = *(workspace->problem);

      adouble* xad = workspace->xad.get();

      int iphase, i, j;

      adouble* parameters;

      int pcount = 0;

      double alpha;

      snprintf(workspace->text,sizeof(workspace->text),"\n>>> Performing statistical analysis of estimated parameters...");

	  psopt_print(workspace,workspace->text);


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

	       parameters    = workspace->parameters[iphase-1].get();

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

      long M= JcT.rows();          // number of decision variables, n_z

      // Null space of the active constraint set, from a rank-revealing QR of Jc^T. With
      // column pivoting Jc^T P = Q R, so the leading rank(Jc) columns of Q span the range
      // of Jc^T and the remaining ones span the null space of Jc. The rank must be
      // computed rather than assumed: the constraint block carries structurally empty rows
      // (the defect array is sized for norder+1 nodes while a local discretization fills
      // norder of them), and the active-bound rows appended by
      // compute_jacobian_of_constraints_with_respect_to_variables may repeat a condition
      // that an event constraint already imposes.

      Eigen::ColPivHouseholderQR<MatrixXd> qr_of_JcT(JcT);
      qr_of_JcT.setThreshold(1.0e-10);

      long rank_Jc = qr_of_JcT.rank();
      long ndof    = M - rank_Jc;    // degrees of freedom of the NLP, n_f = n_z - n_c

      if ( ndof <= 0 ) return false;

      MatrixXd Q = qr_of_JcT.householderQ();
      MatrixXd Z = Q.rightCols(ndof);

      // Calculation of covariance matrix w.r.t all decision variables.
      // See the paper:
      // Kostina et al (2003) "Computation of covariance matrices for constrained parameter estimation
      // problems using LSQR".
      //
      //     C_z = sigma_hat^2 * Z ( Z^T Jr^T Jr Z )^{-1} Z^T
      //
      // Note the inverse: the reduced Gauss-Newton matrix Z^T Jr^T Jr Z is an
      // *information* matrix, and reporting it in place of its inverse gives, roughly, the
      // reciprocal of the covariance.

      MatrixXd& J1 = Jr;
      MatrixXd CC;

      MatrixXd F;

      F = Z.transpose()*J1.transpose()*J1*Z;

      MatrixXd Finv = F.completeOrthogonalDecomposition().pseudoInverse();

      // Residual variance factor. Z (Z^T Jr^T Jr Z)^{-1} Z^T is the covariance only when
      // the residual weights have been set to 1/sigma, so that the weighted residuals have
      // unit variance. When the noise level is not known in advance, which is the usual
      // case, it is estimated from the residuals themselves,
      //
      //     sigma_hat^2 = J* / (N_s - n_f),
      //
      // with J* the sum of squared weighted residuals at the solution, N_s the number of
      // scalar observations and n_f the degrees of freedom. When the weights really are
      // 1/sigma this factor comes out near unity of its own accord, so estimating it is
      // also the right thing to do in that case.

      long n_obs = 0;
      for (int iph=1; iph<=problem.nphases; iph++)
          n_obs += problem.phases(iph).nsamples*problem.phases(iph).nobserved;

      double Jstar  = r.squaredNorm();
      double sigma2 = 1.0;

      if ( n_obs - ndof > 0 ) {
          sigma2 = Jstar/((double)(n_obs - ndof));
      }
      else {
          snprintf(workspace->text,sizeof(workspace->text),
                   "\n>>> Warning: %ld scalar observations against %ld degrees of freedom, so the"
                   "\n    residual variance cannot be estimated; sigma_hat^2 = 1 has been assumed.",
                   n_obs, ndof);
          psopt_print(workspace,workspace->text);
      }

      workspace->pe_sigma_hat = sqrt(sigma2);
      workspace->pe_dof       = ndof;
      workspace->pe_nobs      = n_obs;
      workspace->pe_objective = Jstar;

      CC = sigma2*(Z*Finv*Z.transpose());

      extract_parameter_covariance(Cp, CC, workspace);

      // The Jacobians above are taken with respect to the *scaled* decision variables, so
      // CC and Cp are covariances in scaled units, while the parameter values reported
      // beside them are physical. Convert: with a scaled variable s*p, entry (i,j) is
      // divided by s_i*s_j.

      MatrixXd pscale;
      get_parameter_scale_factors(pscale, workspace);

      for (i=0; i<Cp.rows(); i++) {
          for (j=0; j<Cp.cols(); j++) {
              Cp(i,j) /= ( pscale(i)*pscale(j) );
          }
      }

      // Compute confidence intervals

      alpha = 0.95; // 95% Confidence level

      pcount=0;

      plow.resize( total_number_of_parameters, 1);
      phigh.resize(total_number_of_parameters,1);

      int NN = (int) workspace->pe_nobs;

      // The t distribution carries N_s - n_f degrees of freedom, with n_f the degrees of
      // freedom of the NLP computed above. Counting parameters in place of n_f happens to
      // agree whenever the parameters are the only free quantities, but not in general: a
      // problem with a free initial state, or a free final time, has more.

      long ndof_t = workspace->pe_dof;

      for (iphase=1; iphase<=problem.nphases; iphase++)
      {

      double tt;

      tt = inverse_twotailed_t_cdf(alpha, (int)(NN-ndof_t) );

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

