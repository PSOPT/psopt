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

// Bring std names into this translation unit (formerly leaked via psopt.h).
using namespace std;
using namespace Ipopt;   // for Ipopt return-status enums (Solve_Succeeded)

using namespace Eigen;

void psopt_print(Workspace* workspace, const char* msg)
{
    if (workspace->algorithm->print_level) {
         fprintf(stderr,"%s", msg);
    }
    snprintf(workspace->text,sizeof(workspace->text),"%s"," ");
}


void print_iterations_summary(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace)
{
    // Total output utility: skip-or-stop (per algorithm.on_error) on a failed solve
    // rather than reaching into the empty per-phase arrays.
    if (psopt_solution_failed("print_iterations_summary", solution)) return;

	FILE* outfile  = workspace->psopt_solution_summary_file;
        FILE* outfile2 = workspace->mesh_statistics;
	int jj;

	int sum_n_jacobian_evals   = 0;
	int sum_n_hessian_evals    = 0;
	int sum_n_obj_evals        = 0;
	int sum_n_con_evals        = 0;
	int sum_n_ode_rhs_evals    = 0;
	double sum_CPU_time        = 0.0;


   fprintf(outfile,"\n\n*****************************************************************************************************************");
   fprintf(outfile,"\n************************************* Mesh Refinement Statistics ************************************************");
	fprintf(outfile,"\n*****************************************************************************************************************");

	fprintf(outfile,"\n\nIter\tMethod\tNodes\tNV\tNC\tOEval\tCEval\tJEval\tHEval\tODE RHS\tODE Error\tNLP CPU(sec)");

	for (jj=0;jj< workspace->current_mesh_refinement_iteration;jj++) {

	   fprintf(outfile,"\n%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%e\t%e", jj+1, solution.mesh_stats[jj].method.c_str(), solution.mesh_stats[jj].nnodes, solution.mesh_stats[jj].nvars,
		solution.mesh_stats[jj].ncons, solution.mesh_stats[jj].n_obj_evals, solution.mesh_stats[jj].n_con_evals,
		solution.mesh_stats[jj].n_jacobian_evals, solution.mesh_stats[jj].n_hessian_evals,
		solution.mesh_stats[jj].n_ode_rhs_evals, solution.mesh_stats[jj].epsilon_max,
		solution.mesh_stats[jj].CPU_time);

		sum_n_jacobian_evals   += solution.mesh_stats[jj].n_jacobian_evals;
	        sum_n_hessian_evals    += solution.mesh_stats[jj].n_hessian_evals;
	        sum_n_obj_evals        += solution.mesh_stats[jj].n_obj_evals;
		sum_n_con_evals        += solution.mesh_stats[jj].n_con_evals;
	        sum_n_ode_rhs_evals    += solution.mesh_stats[jj].n_ode_rhs_evals;
	        sum_CPU_time           += solution.mesh_stats[jj].CPU_time;
	}

        fprintf(outfile,"\n__________________________________________________________________________________________________________________\n\n");

	   double diff_CPU_time = solution.cpu_time - sum_CPU_time;

	   fprintf(outfile,"\nAdditional CPU time (sec)\t\t\t\t\t\t\t\t\t%e",
		diff_CPU_time);

		fprintf(outfile,"\nTotals\t-\t-\t-\t-\t%i\t%i\t%i\t%i\t%i\t-\t\t%e",
		sum_n_obj_evals, sum_n_con_evals,
		sum_n_jacobian_evals, sum_n_hessian_evals,
		sum_n_ode_rhs_evals,
		solution.cpu_time);

      fprintf(outfile,"\n__________________________________________________________________________________________________________________\n\n");


	   fclose(outfile);


	fprintf(outfile2,"\n\n*****************************************************************************************************************");
        fprintf(outfile2,"\n************************************* Mesh Refinement Statistics ************************************************");
	fprintf(outfile2,"\n*****************************************************************************************************************");

	fprintf(outfile2,"\n\nIter\tMethod\tNodes\tNV\tNC\tOEval\tCEval\tJEval\tHEval\tODE RHS\tODE Error\tNLP CPU(sec)");

	for (jj=0;jj< workspace->current_mesh_refinement_iteration;jj++) {

	  fprintf(outfile2,"\n%i\t%s\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%e\t%e", jj+1, solution.mesh_stats[jj].method.c_str(), solution.mesh_stats[jj].nnodes, solution.mesh_stats[jj].nvars,
		solution.mesh_stats[jj].ncons, solution.mesh_stats[jj].n_obj_evals, solution.mesh_stats[jj].n_con_evals,
		solution.mesh_stats[jj].n_jacobian_evals, solution.mesh_stats[jj].n_hessian_evals,
		solution.mesh_stats[jj].n_ode_rhs_evals, solution.mesh_stats[jj].epsilon_max,
		solution.mesh_stats[jj].CPU_time);



	}

        fprintf(outfile2,"\n__________________________________________________________________________________________________________________\n\n");


	        fprintf(outfile2,"\nAdditional CPU time (sec)\t\t\t\t\t\t\t\t\t%e",
		diff_CPU_time);

		fprintf(outfile2,"\nTotals\t-\t-\t-\t-\t%i\t%i\t%i\t%i\t%i\t-\t\t%e",
		sum_n_obj_evals, sum_n_con_evals,
		sum_n_jacobian_evals, sum_n_hessian_evals,
		sum_n_ode_rhs_evals,
		solution.cpu_time);

        fprintf(outfile2,"\n__________________________________________________________________________________________________________________\n\n");


	fclose(outfile2);


}


void print_iterations_summary_tex(Prob& problem,Alg& algorithm,Sol& solution, Workspace* workspace)
{
    // Total output utility: skip-or-stop (per algorithm.on_error) on a failed solve
    // rather than reaching into the empty per-phase arrays.
    if (psopt_solution_failed("print_iterations_summary_tex", solution)) return;


	FILE* outfile = workspace->mesh_statistics_tex;
	int jj;

	string caption, label, table_key;

	int sum_n_jacobian_evals   = 0;
	int sum_n_hessian_evals    = 0;
	int sum_n_obj_evals        = 0;
	int sum_n_con_evals        = 0;
	int sum_n_ode_rhs_evals    = 0;
	double sum_CPU_time        = 0.0;



	if (algorithm.print_level==0) return;

	int dotindex = problem.outfilename.find_first_of(".");

	string msname = "mesh_stats_" + problem.outfilename.substr(0,dotindex);

	caption = "\n\\caption{Mesh refinement statistics: " + problem.name + "}";
	label   = "\n\\label{" + msname + "}";
	table_key     = "\n\\newline \\\\ \\emph{Key}: Iter=iteration number, DM= discretization method, M=number of nodes, NV=number of variables, NC=number of constraints, OE=objective evaluations,  \
	              CE = constraint evaluations, JE = Jacobian evaluations, HE = Hessian evaluations, RHS = ODE right hand side \
		      evaluations, $\\epsilon_{\\max}$ = maximum relative ODE error, CPU$_\\mathrm{a}$ = CPU time in seconds spent by NLP algorithm, \
		      CPU$_\\mathrm{b}$ = additional CPU time in seconds spent by PSOPT";



	fprintf(outfile,"\n\\begin{table}");

	fprintf(outfile,"%s",caption.c_str());

        fprintf(outfile,"%s",label.c_str());

	fprintf(outfile,"%s","\n\\renewcommand{\\tabcolsep}{0.15cm}");

        fprintf(outfile,"\n\\small");

        fprintf(outfile,"\n\\begin{tabular}{llllllllllll}");

	fprintf(outfile,"\nIter&DM&M&NV&NC&OE&CE&JE&HE&RHS&$\\epsilon_{\\max}$&CPU$_\\mathrm{a}$ \\\\ \\hline \\\\");

	for (jj=0;jj< workspace->current_mesh_refinement_iteration;jj++) {

	  fprintf(outfile,"\n%i&%s&%i&%i&%i&%i&%i&%i&%i&%i&%.3e&%.3e\\\\", jj+1, solution.mesh_stats[jj].method.c_str(), solution.mesh_stats[jj].nnodes, solution.mesh_stats[jj].nvars,
		solution.mesh_stats[jj].ncons, solution.mesh_stats[jj].n_obj_evals, solution.mesh_stats[jj].n_con_evals,
		solution.mesh_stats[jj].n_jacobian_evals, solution.mesh_stats[jj].n_hessian_evals,
		solution.mesh_stats[jj].n_ode_rhs_evals, solution.mesh_stats[jj].epsilon_max,
		solution.mesh_stats[jj].CPU_time);

		sum_n_jacobian_evals   += solution.mesh_stats[jj].n_jacobian_evals;
	        sum_n_hessian_evals    += solution.mesh_stats[jj].n_hessian_evals;
	        sum_n_obj_evals        += solution.mesh_stats[jj].n_obj_evals;
		sum_n_con_evals        += solution.mesh_stats[jj].n_con_evals;
	        sum_n_ode_rhs_evals    += solution.mesh_stats[jj].n_ode_rhs_evals;
	        sum_CPU_time           += solution.mesh_stats[jj].CPU_time;
	}

	double diff_CPU_time = solution.cpu_time - sum_CPU_time;

      fprintf(outfile,"\n\\hline");
		fprintf(outfile,"\nCPU$_\\mathrm{b}$ &-&-&-&-&-&-&-&-&-&-&%.3e\\\\",
		diff_CPU_time);

		fprintf(outfile,"\n-&-&-&-&-&%i&%i&%i&%i&%i&-&%.3e\\\\",
		sum_n_obj_evals, sum_n_con_evals,
		sum_n_jacobian_evals, sum_n_hessian_evals,
		sum_n_ode_rhs_evals,
		solution.cpu_time);

      fprintf(outfile,"\n\\end{tabular}");

	   fprintf(outfile,"%s",table_key.c_str());

      fprintf(outfile,"\n\\normalsize");

      fprintf(outfile,"\n\\end{table}");

	fclose(outfile);

}



void print_psopt_summary(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
    // Total output utility: skip-or-stop (per algorithm.on_error) on a failed solve
    // rather than reaching into the empty per-phase arrays.
    if (psopt_solution_failed("print_psopt_summary", solution)) return;
    FILE* outfile;
    string auxstr;
    string filename;
    string mesh_stats_file = "mesh_statistics.txt";
    int i;


    if ( !algorithm.print_level ) return;

    if ( problem.outfilename == "" )
       filename = "psopt.txt";
    else
       filename = problem.outfilename;
    outfile = fopen( filename.c_str(), "w");
    if (outfile == NULL) {
        // Previously unchecked: a failed open here fprintf'd to NULL (a segfault).
        // Warn and skip the summary instead; the solve itself is unaffected.
        if (PSOPT_extras::PrintLevel())
            fprintf(stderr,
                "\n**** ====> PSOPT: print_psopt_summary could not open \"%s\" for "
                "writing; skipping the summary. <====\n\n", filename.c_str());
        return;
    }
    fprintf(outfile,"\nPSOPT results summary");
    fprintf(outfile,"\n=====================\n");

    auxstr = "Problem: \t\t\t\t\t" + problem.name ;
    fprintf(outfile,"\n%s", auxstr.c_str());
    fprintf(outfile,"\nCPU time (seconds):\t\t\t\t%e", solution.cpu_time);
    fprintf(outfile,"\nNLP solver used: \t\t\t\t%s", algorithm.nlp_method.c_str());
    fprintf(outfile,"\nPSOPT release number: \t\t\t\t%s", PSOPT_RELEASE_STRING );
    fprintf(outfile,"\nPSOPT build date: \t\t\t\t%s", PSOPT_BUILD_DATE );

    fprintf(outfile,"\nDate and time of this run: \t\t\t%s", solution.end_date_and_time.c_str() );

    if ( algorithm.nlp_method == "IPOPT") {
        if (solution.nlp_return_code  == (int) Solve_Succeeded) {
            fprintf(outfile,"\nOptimal (unscaled) cost function value: \t%e", solution.cost);
            for (i=0;i < problem.nphases; i++) {

            	fprintf(outfile,"\nPhase %i endpoint cost function value:\t\t%e",i+1, solution.endpoint_cost[i]);
            	fprintf(outfile,"\nPhase %i integrated part of the cost: \t\t%e", i+1, solution.integrated_cost[i]);
            	fprintf(outfile,"\nPhase %i initial time:\t\t\t\t%e",i+1, (solution.nodes[i])(0));  // EIGEN_UPDATE
                fprintf(outfile,"\nPhase %i final time:\t\t\t\t%e",i+1, (solution.nodes[i])(length(solution.nodes[i])-1));
		fprintf(outfile,"\nPhase %i maximum relative local error:\t\t%e", i+1, Max(solution.relative_errors[i]) );

	    }
            auxstr = "The problem has been solved!";
        }
        else {
           auxstr = "*** The problem FAILED! - see screen output";
           fprintf(outfile,"\nReturned (unscaled) cost function value: \t%e", solution.cost);
            for (i=0;i < problem.nphases; i++) {

            	fprintf(outfile,"\nPhase %i endpoint cost function value: \t\t%e",i+1, solution.endpoint_cost[i]);
            	fprintf(outfile,"\nPhase %i integrated part of the cost: \t\t%e", i+1, solution.integrated_cost[i]);
            	fprintf(outfile,"\nPhase %i initial time:\t\t\t\t%e",i+1, (solution.nodes[i])(0));
                fprintf(outfile,"\nPhase %i final time:\t\t\t\t%e",i+1, (solution.nodes[i])(length(solution.nodes[i])-1));
		fprintf(outfile,"\nPhase %i maximum relative local error\t\t%e", i+1, Max(solution.relative_errors[i]) );

	    }
        }
    }

    if ( algorithm.nlp_method == "SNOPT") {
        if (solution.nlp_return_code  == 0) {
            fprintf(outfile,"\nOptimal (unscaled) cost function value: \t%e", solution.cost);
            for (i=0;i < problem.nphases; i++) {

            	fprintf(outfile,"\nPhase %i endpoint cost function value:\t\t%e",i+1, solution.endpoint_cost[i]);
            	fprintf(outfile,"\nPhase %i integrated part of the cost: \t\t%e", i+1, solution.integrated_cost[i]);
            	fprintf(outfile,"\nPhase %i initial time:\t\t\t\t%e",i+1, (solution.nodes[i])(0));
                fprintf(outfile,"\nPhase %i final time:\t\t\t\t%e",i+1, (solution.nodes[i])(length(solution.nodes[i])-1));
		fprintf(outfile,"\nPhase %i maximum relative local error\t\t%e", i+1, Max(solution.relative_errors[i]) );

	    }
            auxstr = "Finished successfully";
        }
        else {
           auxstr = "*** The problem FAILED! - see screen output";
           fprintf(outfile,"\nReturned (unscaled) cost function value: \t%e", solution.cost);
            for (i=0;i < problem.nphases; i++) {

            	fprintf(outfile,"\nPhase %i endpoint cost function value: \t\t%e",i+1, solution.endpoint_cost[i]);
            	fprintf(outfile,"\nPhase %i integrated part of the cost: \t\t%e", i+1, solution.integrated_cost[i]);
            	fprintf(outfile,"\nPhase %i initial time:\t\t\t\t%e",i+1, (solution.nodes[i])(0));
                fprintf(outfile,"\nPhase %i final time:\t\t\t\t%e",i+1, (solution.nodes[i])(length(solution.nodes[i])-1));
		fprintf(outfile,"\nPhase %i maximum relative local error\t\t%e", i+1, Max(solution.relative_errors[i]) );

	    }
        }
    }

    fprintf(outfile,"\nNLP solver reports: \t\t\t\t%s\n", auxstr.c_str() );

    fclose(outfile);


#ifndef WIN32
    auxstr = "cat " + mesh_stats_file;
#else
	auxstr = "type " + mesh_stats_file;
#endif

    if ( system(auxstr.c_str()) != 0 )
        fprintf(stderr, "PSOPT warning: command failed: %s\n", auxstr.c_str());


#ifndef WIN32
    auxstr = "cat " + filename;
#else
	auxstr = "type " + filename;
#endif

    if ( system(auxstr.c_str()) != 0 )
        fprintf(stderr, "PSOPT warning: command failed: %s\n", auxstr.c_str());



}


// ============================================================================
//  Save_to_json_file
// ----------------------------------------------------------------------------
//  Writes the full PSOPT solution to a single, language-neutral JSON file so
//  the result can be loaded for post-processing in Python, MATLAB, R, etc.
//  with one read.
//
//  Schema:
//    {
//      "summary": {
//        "problem_name", "nlp_solver", "psopt_release", "psopt_build_date",
//        "run_date_and_time", "cpu_time_seconds", "nlp_return_code",
//        "solved" (bool), "solver_message", "optimal_cost", "n_phases"
//      },
//      "phases": [
//        {
//          "phase", "n_nodes", "n_states", "n_controls", "n_parameters",
//          "endpoint_cost", "integrated_cost", "initial_time", "final_time",
//          "max_relative_local_error",
//          "time":            [ n_nodes ],
//          "states":          [ n_states ][ n_nodes ],
//          "controls":        [ n_controls ][ n_nodes ],
//          "costates":        [ n_states ][ n_nodes ],
//          "relative_errors": [ rows ][ cols ],
//          "parameters":      [ n_parameters ]
//        }, ...
//      ]
//    }
//
//  Orientation: every trajectory matrix is written as an array of rows, where
//  row i is the trajectory of variable i+1 across the nodes -- shape
//  [n_variables][n_nodes], the same orientation PSOPT uses internally and that
//  Save(const MatrixXd&, ...) writes. In NumPy, json.load(...) then
//  np.array(phase["states"]) yields an (n_states, n_nodes) array directly.
//
//  Notes:
//   - Numbers are written with full double precision (%.17g) so values
//     round-trip exactly, unlike the %g used by the tab-delimited Save().
//   - JSON cannot represent NaN/Inf; any non-finite value is written as null
//     so the file stays valid for standard parsers.
// ============================================================================

static void psopt_json_string(FILE* fp, const string& s)
{
    fputc('"', fp);
    for (string::size_type k = 0; k < s.size(); k++) {
        unsigned char c = (unsigned char) s[k];
        switch (c) {
            case '"':  fputs("\\\"", fp); break;
            case '\\': fputs("\\\\", fp); break;
            case '\b': fputs("\\b",  fp); break;
            case '\f': fputs("\\f",  fp); break;
            case '\n': fputs("\\n",  fp); break;
            case '\r': fputs("\\r",  fp); break;
            case '\t': fputs("\\t",  fp); break;
            default:
                if (c < 0x20) fprintf(fp, "\\u%04x", (unsigned int) c);
                else          fputc((int) c, fp);
        }
    }
    fputc('"', fp);
}

static void psopt_json_number(FILE* fp, double x)
{
    // JSON has no NaN/Inf representation; emit null so the file stays valid.
    if (std::isfinite(x)) fprintf(fp, "%.17g", x);
    else                  fputs("null", fp);
}

// Flatten any MatrixXd to a 1-D JSON array (used for the time and parameter vectors).
static void psopt_json_vector(FILE* fp, const MatrixXd& v)
{
    long n = (long) v.size();
    fputc('[', fp);
    for (long k = 0; k < n; k++) {
        if (k) fputc(',', fp);
        psopt_json_number(fp, v(k));
    }
    fputc(']', fp);
}

// Write a MatrixXd as a JSON array of rows: [ [row 0 ...], [row 1 ...], ... ].
static void psopt_json_matrix(FILE* fp, const MatrixXd& m)
{
    fputc('[', fp);
    for (long i = 0; i < (long) m.rows(); i++) {
        if (i) fputc(',', fp);
        fputs("\n        [", fp);
        for (long j = 0; j < (long) m.cols(); j++) {
            if (j) fputc(',', fp);
            psopt_json_number(fp, m(i, j));
        }
        fputc(']', fp);
    }
    fputc(']', fp);
}

void Save_to_json_file(const string& filename, Prob& problem, Sol& solution, Alg& algorithm)
{
    // Total output utility: skip-or-stop (per algorithm.on_error) on a failed solve
    // rather than reaching into the empty per-phase arrays.
    if (psopt_solution_failed("Save_to_json_file", solution)) return;
    FILE* fp = fopen(filename.c_str(), "w");
    if (fp == NULL) {
        // Non-throwing I/O failure: warn and skip, rather than throwing ErrorHandler
        // (which escaped to std::terminate from user code). The write is skipped.
        if (PSOPT_extras::PrintLevel())
            fprintf(stderr,
                "\n**** ====> PSOPT: Save_to_json_file could not open \"%s\" for "
                "writing; skipping (file not written). <====\n\n", filename.c_str());
        return;
    }

    // Determine success exactly as print_psopt_summary() does.
    bool solved = false;
    const char* message = "*** The problem FAILED! - see screen output";
    if (algorithm.nlp_method == "IPOPT") {
        solved = (solution.nlp_return_code == (int) Solve_Succeeded);
        if (solved) message = "The problem has been solved!";
    }
    else if (algorithm.nlp_method == "SNOPT") {
        solved = (solution.nlp_return_code == 0);
        if (solved) message = "Finished successfully";
    }

    fputs("{\n", fp);

    // ---- run summary (mirrors print_psopt_summary) ----
    fputs("  \"summary\": {\n", fp);
    fputs("    \"problem_name\": ",      fp); psopt_json_string(fp, problem.name);                  fputs(",\n", fp);
    fputs("    \"nlp_solver\": ",        fp); psopt_json_string(fp, algorithm.nlp_method);          fputs(",\n", fp);
    fputs("    \"psopt_release\": ",     fp); psopt_json_string(fp, string(PSOPT_RELEASE_STRING));  fputs(",\n", fp);
    fputs("    \"psopt_build_date\": ",  fp); psopt_json_string(fp, string(PSOPT_BUILD_DATE));      fputs(",\n", fp);
    fputs("    \"run_date_and_time\": ", fp);
    {
        // solution.end_date_and_time comes from ctime() and carries a trailing
        // newline; strip it so the JSON value is clean.
        string run_dt = solution.end_date_and_time;
        while (!run_dt.empty() && (run_dt[run_dt.size()-1] == '\n' || run_dt[run_dt.size()-1] == '\r'))
            run_dt.erase(run_dt.size()-1);
        psopt_json_string(fp, run_dt);
    }
    fputs(",\n", fp);
    fputs("    \"cpu_time_seconds\": ",  fp); psopt_json_number(fp, solution.cpu_time);             fputs(",\n", fp);
    fprintf(fp, "    \"nlp_return_code\": %d,\n", solution.nlp_return_code);
    fprintf(fp, "    \"solved\": %s,\n", solved ? "true" : "false");
    fputs("    \"solver_message\": ",    fp); psopt_json_string(fp, string(message));               fputs(",\n", fp);
    fputs("    \"optimal_cost\": ",      fp); psopt_json_number(fp, solution.cost);                 fputs(",\n", fp);
    fprintf(fp, "    \"n_phases\": %d\n", problem.nphases);
    fputs("  },\n", fp);

    // ---- per-phase trajectories and costs ----
    fputs("  \"phases\": [\n", fp);
    for (int i = 1; i <= problem.nphases; i++) {
        MatrixXd& time     = solution.get_time_in_phase(i);
        MatrixXd& states   = solution.get_states_in_phase(i);
        MatrixXd& controls = solution.get_controls_in_phase(i);
        MatrixXd& costates = solution.get_dual_costates_in_phase(i);
        MatrixXd& relerr   = solution.get_relative_local_error_in_phase(i);
        MatrixXd& params   = solution.get_parameters_in_phase(i);

        long nnodes = length(time);

        fputs("    {\n", fp);
        fprintf(fp, "      \"phase\": %d,\n", i);
        fprintf(fp, "      \"n_nodes\": %ld,\n", nnodes);
        fprintf(fp, "      \"n_states\": %ld,\n", (long) states.rows());
        fprintf(fp, "      \"n_controls\": %ld,\n", (long) controls.rows());
        fprintf(fp, "      \"n_parameters\": %ld,\n", (long) params.size());
        fputs("      \"endpoint_cost\": ",   fp); psopt_json_number(fp, solution.endpoint_cost[i-1]);   fputs(",\n", fp);
        fputs("      \"integrated_cost\": ", fp); psopt_json_number(fp, solution.integrated_cost[i-1]); fputs(",\n", fp);
        fputs("      \"initial_time\": ",    fp); psopt_json_number(fp, (nnodes > 0) ? time(0) : 0.0);  fputs(",\n", fp);
        fputs("      \"final_time\": ",      fp); psopt_json_number(fp, (nnodes > 0) ? time(nnodes-1) : 0.0); fputs(",\n", fp);
        fputs("      \"max_relative_local_error\": ", fp);
              psopt_json_number(fp, (relerr.size() > 0) ? Max(relerr) : 0.0);                          fputs(",\n", fp);

        fputs("      \"time\": ",            fp); psopt_json_vector(fp, time);     fputs(",\n", fp);
        fputs("      \"states\": ",          fp); psopt_json_matrix(fp, states);   fputs(",\n", fp);
        fputs("      \"controls\": ",        fp); psopt_json_matrix(fp, controls); fputs(",\n", fp);
        fputs("      \"costates\": ",        fp); psopt_json_matrix(fp, costates); fputs(",\n", fp);
        fputs("      \"relative_errors\": ", fp); psopt_json_matrix(fp, relerr);   fputs(",\n", fp);
        if (solution.smoothness != NULL) {   // populated only when diagnostic_level > 0
            fputs("      \"smoothness\": ", fp); psopt_json_matrix(fp, solution.smoothness[i-1]); fputs(",\n", fp);
        }
        if (solution.stationarity_residual != NULL) {   // populated only when diagnostic_level >= 2
            fputs("      \"stationarity_residual\": ", fp);
            psopt_json_matrix(fp, solution.stationarity_residual[i-1]); fputs(",\n", fp);
        }
        fputs("      \"parameters\": ",      fp); psopt_json_vector(fp, params);   fputs("\n", fp);

        fputs((i < problem.nphases) ? "    },\n" : "    }\n", fp);
    }
    fputs("  ]\n", fp);
    fputs("}\n", fp);

    fclose(fp);
}





void print_algorithm_summary(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
    // Total output utility: skip-or-stop (per algorithm.on_error) on a failed solve
    // rather than reaching into the empty per-phase arrays.
    if (psopt_solution_failed("print_algorithm_summary", solution)) return;
    FILE * outfile = workspace->psopt_solution_summary_file;
    string amrtype;

    fprintf(outfile,"\n*********************************************************************************************************");
    fprintf(outfile,"\n***************************************** PSOPT RUN SUMMARY *********************************************");
    fprintf(outfile,"\n*********************************************************************************************************");
    fprintf(outfile,"\n\n PROBLEM NAME:                   %s\n", problem.name.c_str() );
    fprintf(outfile,"\n PSOPT RELEASE NUMBER:           %s\n", PSOPT_RELEASE_STRING );
    fprintf(outfile,"\n PSOPT BUILD DATE:               %s\n", PSOPT_BUILD_DATE );
    fprintf(outfile,"\n DATE AND TIME OF RUN START:     %s", solution.start_date_and_time.c_str() );
    fprintf(outfile,"\n DATE AND TIME OF RUN END:       %s", solution.end_date_and_time.c_str()   );


    fprintf(outfile,"\n*********************************************************************************************************");
    fprintf(outfile,"\n***************************************** ALGORITHM OPTIONS *********************************************");
    fprintf(outfile,"\n*********************************************************************************************************");
    fprintf(outfile,"\nCOLLOCATION METHOD:             %s", algorithm.collocation_method.c_str()   );
    if (use_global_collocation(algorithm)) {
    fprintf(outfile,"\nDIFFERENTIATION MATRIX:         %s", algorithm.diff_matrix.c_str()   );
    }
    fprintf(outfile,"\nNLP METHOD:                     %s", algorithm.nlp_method.c_str()   );
    if (algorithm.nlp_method == "IPOPT") {
    fprintf(outfile,"\nHESSIAN OPTION:                 %s", algorithm.hessian.c_str()   );
    }
    fprintf(outfile,"\nNLP TOLERANCE:                  %e", algorithm.nlp_tolerance   );
    fprintf(outfile,"\nNLP MAX ITERATIONS:             %i", algorithm.nlp_iter_max   );
    fprintf(outfile,"\nSCALING METHOD:                 %s", algorithm.scaling.c_str()   );
    fprintf(outfile,"\nDEFECT SCALING:                 %s", algorithm.defect_scaling.c_str()   );
    fprintf(outfile,"\nDIFFERENTIATION:                %s", algorithm.derivatives.c_str()   );

    if (algorithm.mesh_refinement == "manual") {
        amrtype = "";
    }
    else if (use_global_collocation(algorithm) ) {
        amrtype = " (global variant)";
    }
    else {
        amrtype = " (local variant)";
    }

    amrtype = algorithm.mesh_refinement + amrtype;

    fprintf(outfile,"\nMESH REFINEMENT:                %s", amrtype.c_str()  );
    if (algorithm.mesh_refinement == "automatic" )  {
    fprintf(outfile,"\nMESH REF. ODE TOLERANCE:        %e", algorithm.ode_tolerance   );
    fprintf(outfile,"\nMESH REF. MAX ITERATIONS:       %i", algorithm.mr_max_iterations   );
    fprintf(outfile,"\nMESH REF. MAX GROWTH FACTOR:    %e", algorithm.mr_max_growth_factor   );
      if (use_global_collocation(algorithm)) {
    fprintf(outfile,"\nMESH REF. MIN ORDER (N_min):    %i", algorithm.mr_min_order   );
    fprintf(outfile,"\nMESH REF. MAX ORDER (N_max):    %i", algorithm.mr_max_order   );
      }
      else {
    fprintf(outfile,"\nMESH REF. KAPPA:                %e", algorithm.mr_kappa   );
    fprintf(outfile,"\nMESH REF. M1:                   %i", algorithm.mr_M1   );
    fprintf(outfile,"\nMESH REF. SWITCH ORDER AT ITER: %i", algorithm.switch_order   );
      }

    }




}

void print_solution_summary(Prob& problem, Alg& algorithm, Sol& solution, Workspace* workspace)
{
    // Total output utility: skip-or-stop (per algorithm.on_error) on a failed solve
    // rather than reaching into the empty per-phase arrays.
    if (psopt_solution_failed("print_solution_summary", solution)) return;
    FILE * outfile = workspace->psopt_solution_summary_file;

    int i,j,k;

    string auxstr;


    fprintf(outfile,"\n*********************************************************************************************************");
    fprintf(outfile,"\n*****************************************  SOLUTION SUMMARY *********************************************");
    fprintf(outfile,"\n*********************************************************************************************************\n");





    fprintf(outfile,"\nTotal CPU time (seconds):\t\t\t%e", solution.cpu_time);

    if ( algorithm.nlp_method == "IPOPT") {
        if (solution.nlp_return_code  == (int) Solve_Succeeded) {
            fprintf(outfile,"\nOptimal (unscaled) cost function value: \t%e", solution.cost);
            for (i=0;i < problem.nphases; i++) {
            	fprintf(outfile,"\nPhase %i endpoint cost function value:\t\t%e",i+1, solution.endpoint_cost[i]);
            	fprintf(outfile,"\nPhase %i integrated part of the cost: \t\t%e", i+1, solution.integrated_cost[i]);
            	fprintf(outfile,"\nPhase %i initial time:\t\t\t\t%e",i+1, (solution.nodes[i])(0));
                fprintf(outfile,"\nPhase %i final time:\t\t\t\t%e",i+1, (solution.nodes[i])(length(solution.nodes[i])-1));
		fprintf(outfile,"\nPhase %i maximum relative local error:\t\t%e", i+1, Max(solution.relative_errors[i]) );

	    }
            auxstr = "The problem has been solved!";
        }
        else {
           auxstr = "*** The problem FAILED! - see screen output";
           fprintf(outfile,"\nReturned (unscaled) cost function value: \t%e", solution.cost);
            for (i=0;i < problem.nphases; i++) {
            	fprintf(outfile,"\nPhase %i endpoint cost function value: \t\t%e",i+1, solution.endpoint_cost[i]);
            	fprintf(outfile,"\nPhase %i integrated part of the cost: \t\t%e", i+1, solution.integrated_cost[i]);
            	fprintf(outfile,"\nPhase %i initial time:\t\t\t\t%e",i+1, (solution.nodes[i])(0));
                fprintf(outfile,"\nPhase %i final time:\t\t\t\t%e",i+1, (solution.nodes[i])(length(solution.nodes[i])-1));
		fprintf(outfile,"\nPhase %i maximum relative local error\t\t%e", i+1, Max(solution.relative_errors[i]) );

	    }
        }
    }

    if ( algorithm.nlp_method == "SNOPT") {
        if (solution.nlp_return_code  == 0) {
            fprintf(outfile,"\nOptimal (unscaled) cost function value: \t%e", solution.cost);
            for (i=0;i < problem.nphases; i++) {
            	fprintf(outfile,"\nPhase %i endpoint cost function value:\t\t%e",i+1, solution.endpoint_cost[i]);
            	fprintf(outfile,"\nPhase %i integrated part of the cost: \t\t%e", i+1, solution.integrated_cost[i]);
            	fprintf(outfile,"\nPhase %i initial time:\t\t\t\t%e",i+1, (solution.nodes[i])(0));
                fprintf(outfile,"\nPhase %i final time:\t\t\t\t%e",i+1, (solution.nodes[i])(length(solution.nodes[i])-1));
		fprintf(outfile,"\nPhase %i maximum relative local error\t\t%e", i+1, Max(solution.relative_errors[i]) );

	    }
            auxstr = "Finished successfully";
        }
        else {
           auxstr = "*** The problem FAILED! - see screen output";
           fprintf(outfile,"\nReturned (unscaled) cost function value: \t%e", solution.cost);
            for (i=0;i < problem.nphases; i++) {
            	fprintf(outfile,"\nPhase %i endpoint cost function value: \t\t%e",i+1, solution.endpoint_cost[i]);
            	fprintf(outfile,"\nPhase %i integrated part of the cost: \t\t%e", i+1, solution.integrated_cost[i]);
            	fprintf(outfile,"\nPhase %i initial time:\t\t\t\t%e",i+1, (solution.nodes[i])(0));
                fprintf(outfile,"\nPhase %i final time:\t\t\t\t%e",i+1, (solution.nodes[i])(length(solution.nodes[i])-1));
		fprintf(outfile,"\nPhase %i maximum relative local error\t\t%e", i+1, Max(solution.relative_errors[i]) );

	    }
        }
    }

    fprintf(outfile,"\nNLP solver reports: \t\t\t\t%s", auxstr.c_str() );
    fprintf(outfile,"\nNLP return code: \t\t\t\t%i", solution.nlp_return_code );
    fprintf(outfile,"\nInternal scaling for objective function \t%e", problem.scale.objective );


    for(i=1;i<=problem.nphases;i++) {


    fprintf(outfile,"\n*****************************************  PHASE %i GRID POINTS (TIME) ***********************************", i);

      for (k=0; k< problem.phases(i).current_number_of_intervals+1;k++) { // EIGEN_UPDATE

	    fprintf(outfile,"\n%e", solution.nodes[i-1](k));


      }


    if (problem.phases(i).ncontrols>0) {


    fprintf(outfile,"\n*****************************************  PHASE %i CONTROLS *******************************************", i);

      for (k=0; k< problem.phases(i).current_number_of_intervals+1;k++) { // EIGEN_UPDATE
	  fprintf(outfile,"\n");
	  for(j=0;j<problem.phases(i).ncontrols;j++) {
	    fprintf(outfile,"%e\t", solution.controls[i-1](j,k));
	  }

      }

    }


    if (problem.phases(i).nstates>0) {

    fprintf(outfile,"\n*****************************************  PHASE %i STATES *******************************************", i);

      for (k=0; k< problem.phases(i).current_number_of_intervals+1;k++) {  // EIGEN_UPDATE
	  fprintf(outfile,"\n");
	  for(j=0;j<problem.phases(i).nstates;j++) {
	    fprintf(outfile,"%e\t", solution.states[i-1](j,k));
	  }

      }
    }

    if (problem.phases(i).nparameters>0) {
    fprintf(outfile,"\n*****************************************  PHASE %i PARAMETERS *******************************************", i);


    fprintf(outfile,"\n");
    for(j=0;j<problem.phases(i).nparameters;j++) {  // EIGEN_UPDATE
	    fprintf(outfile,"\n%e", solution.parameters[i-1](j));
    }


    }


    }

    if ( problem.observation_function!=NULL && algorithm.parameter_statistics == "yes") {

    fprintf(outfile,"\n****************************************** PARAMETER STATISTICS *****************************************");

    int nparam = get_total_number_of_parameters(problem);

    MatrixXd Cp(nparam,nparam), plow(nparam,1), phigh(nparam,1), p(nparam,1), r;

    bool peout = compute_parameter_statistics(Cp, p, plow, phigh,r, workspace);

    int pcount = 0;

    if (peout) {



      fprintf(outfile,"\n>>>>> Parameter covariance matrix (all phases) ");

      for(i=0; i<Cp.rows();i++) { // EIGEN_UPDATE
        fprintf(outfile,"\n");
        for(j=0; j<Cp.cols();j++) {
             fprintf(outfile,"%e\t\t",Cp(i,j));
        }
      }

     FullPivLU<MatrixXd> lu_decomp(Cp);  
        
     fprintf(outfile,"\n\n>>>>> Rank of parameter covariance matrix: %li ", lu_decomp.rank() );

     fprintf(outfile,"\n\n>>> 95 percent statistical confidence limits on estimated parameters ");
     fprintf(outfile,"\nPhase\tParameter\t(Low Confidence Limit) \t(Value) \t\t(High Confidence Limit)");

      pcount = 0;
      
      for (int iphase=1; iphase<=problem.nphases; iphase++)
      {


          for(i=0;i< problem.phases(iphase).nparameters;i++) {

            j= pcount+i; // EIGEN_UPDATE

            fprintf(outfile,"\n%i\t%i\t\t%e\t\t%e\t\t%e", iphase, i+1, plow(j), p(j), phigh(j));

          }

          pcount += problem.phases(iphase).nparameters;


      }

    }

    else {

       fprintf(outfile,"\n *** Unable to compute confidence limits for estimated parameters");

    }

      pcount = 0;

      fprintf(outfile,"\n\n>>> Residual of observed variables");

      fprintf(outfile,"\nPhase\tTime\t\t (Residuals, each column corresponds to an observed variable)");

      for (int iphase=1; iphase<=problem.nphases; iphase++)
      {

          for(i=0;i< problem.phases(iphase).nsamples;i++) {

                  fprintf(outfile,"\n%i\t%e", iphase, problem.phases(iphase).observation_nodes(i) ); // EIGEN_UPDATE

              for(j=0;j< problem.phases(iphase).nobserved; j++ ) {

                  k= pcount+i*problem.phases(iphase).nobserved+j;

                  fprintf(outfile,"\t%e", r(k));

              }


          }

          pcount += problem.phases(iphase).nsamples*problem.phases(iphase).nobserved;


      }


    }


}


void print_constraint_summary(Prob& problem, Sol& solution, Workspace* workspace)
{
    // Total output utility: skip-or-stop (per algorithm.on_error) on a failed solve
    // rather than reaching into the empty per-phase arrays.
    if (psopt_solution_failed("print_constraint_summary", solution)) return;

    int i, j, k, l;

    FILE* outfile = workspace->psopt_solution_summary_file;


    MatrixXd& X = *workspace->Xip;

    MatrixXd& G  = *workspace->Gip;

    MatrixXd& xlb    = *workspace->xlb;

    MatrixXd& xub    = *workspace->xub;

    MatrixXd& cs = *workspace->constraint_scaling;

    int nvars = get_number_nlp_vars(problem, workspace);

    int ncons = get_number_nlp_constraints(problem, workspace);

    double* g_l = new double[ncons];
    double* g_u = new double[ncons];

    double tol    = (workspace->algorithm)->nlp_tolerance;


    for(i=0;i<nvars;i++) {  // EIGEN_UPDATE
       X(i) = solution.xad[i].value();
    }

    gg_num(X, &G, workspace);

    get_constraint_bounds(g_l, g_u, workspace);

    int phase_offset  = 0;

       fprintf(outfile,"\n*********************************************************************************************************");
       fprintf(outfile,"\n*************************************** NLP CONSTRAINT SUMMARY ******************************************");
       fprintf(outfile,"\n*********************************************************************************************************\n");


    for(i=0;i< problem.nphases; i++)
    {
	int norder    = problem.phase[i].current_number_of_intervals;
	int nstates   = problem.phase[i].nstates;
	int nevents   = problem.phase[i].nevents;
	int npath     = problem.phase[i].npath;
	int offset;
   int path_offset = phase_offset+nstates*(norder+1)+nevents;
	int ncons_phase_i = get_ncons_phase_i(problem,i, workspace);

       fprintf(outfile,"\n******************************* PHASE %i SCALED DIFFERENTIAL DEFECT CONSTRAINTS **************************", i+1);
       fprintf(outfile,"\nLOWER\t\tG(X)\t\tUPPER\t\tSCALING\t\tDESCRIPTION");

       for(k=0; k<norder+1; k++)  // EIGEN_UPDATE
       {
                for (j=0; j<problem.phase[i].nstates; j++) {
	   	     l = phase_offset+(k)*nstates+j;
		     fprintf(outfile,"\n%e\t%e\t%e\t%e\tDEFECT FOR STATE DERIV. %i NODE %i", g_l[l], G(l), g_u[l], cs(l), j, k);
		     if ( g_l[l]-tol>G(l) || G(l)>g_u[l]+tol ) fprintf(outfile," **BOUND VIOLATED");
       		     if ( (g_l[l]-tol<G(l) && G(l)<g_l[l]+tol) || (g_u[l]-tol<G(l) && G(l)<g_u[l]+tol) ) fprintf(outfile," **AT BOUND");

		}
       }
       fprintf(outfile,"\n************************************* PHASE %i SCALED PATH CONSTRAINTS ************************************", i+1);
       fprintf(outfile,"\nLOWER\t\tG(X)\t\tUPPER\t\tSCALING\t\tDESCRIPTION");
       for(k=0; k<norder+1; k++) // EIGEN_UPDATE
       {
	        for (j=0; j<npath; j++)
	        {
		     l = path_offset + (k)*npath + j;
		     fprintf(outfile,"\n%e\t%e\t%e\t%e\tPATH CONSTRAINT %i NODE %i", g_l[l], G(l), g_u[l], cs(l),j,k);
		     if ( g_l[l]-tol>G(l) || G(l)>g_u[l]+tol ) fprintf(outfile," **BOUND VIOLATED");
     		     if ( (g_l[l]-tol<G(l) && G(l)<g_l[l]+tol) || (g_u[l]-tol<G(l) && G(l)<g_u[l]+tol) ) fprintf(outfile," **AT BOUND");

                }
       }
       offset = phase_offset+nstates*(norder+1);
       fprintf(outfile,"\n************************************ PHASE %i SCALED EVENT CONSTRAINTS  ***********************************",i+1);
       fprintf(outfile,"\nLOWER\t\tG(X)\t\tUPPER\t\tSCALING\t\tDESCRIPTION");

	for (k=0; k<nevents;k++) {
		j = offset + k;
		fprintf(outfile,"\n%e\t%e\t%e\t%e\tEVENT CONSTRAINT %i", g_l[j], G(j+1), g_u[j], cs(j+1),k+1);
	        if ( g_l[j]-tol>G(j) || G(j)>g_u[j]+tol ) fprintf(outfile," **BOUND VIOLATED");
		if ( (g_l[j]-tol<G(j) && G(j)<g_l[j]+tol) || (g_u[j]-tol<G(j) && G(j)<g_u[j]+tol)) fprintf(outfile," **AT BOUND");

  	}
	// Add tf >= t0 constraint [ t0MIN-tfMAX <= t0-tf <= 0 ]

       fprintf(outfile,"\n*********************************** PHASE %i SCALED tf >= t0 CONSTRAINT ***********************************",i+1);
       fprintf(outfile,"\nLOWER\t\tG(X)\t\tUPPER\t\tSCALING");
       j = phase_offset + ncons_phase_i - 1;
       fprintf(outfile,"\n%e\t%e\t%e\t%e", g_l[j], G(j), g_u[j], cs(j));

       if ( g_l[j]-tol>G(j) || G(j)>g_u[j]+tol ) fprintf(outfile," **BOUND VIOLATED");
       if ( (g_l[j]-tol<G(j) && G(j)<g_l[j]+tol) || (g_u[j]-tol<G(j) && G(j)<g_u[j]+tol)) fprintf(outfile," **AT BOUND");

       phase_offset += ncons_phase_i;
    }

       fprintf(outfile,"\n************************************* SCALED PHASE LINKAGE CONSTRAINTS ***********************************");
       fprintf(outfile,"\nLOWER\t\tG(X)\t\tUPPER\t\tSCALING\t\tDESCRIPTION");

    for(j=0;j<problem.nlinkages;j++)
    {
     	int l = phase_offset+j;
	   fprintf(outfile,"\n%e\t%e\t%e\t%e\tLINKAGE CONSTRAINT %i", g_l[l], G(l), g_u[l], cs(l),j);
	   if ( g_l[l]-tol>G(l) || G(l)>g_u[l]+tol ) fprintf(outfile," **BOUND VIOLATED");
        if ( (g_l[l]-tol<G(l) && G(l)<g_l[l]+tol) || (g_u[l]-tol<G(l) && G(l)<g_u[l]+tol) ) fprintf(outfile," **AT BOUND");

    }




      fprintf(outfile,"\n**************************************UNSCALED NLP DECISION VARIABLES *************************************");
      fprintf(outfile,"\nJ\tLOWER\t\tX(J)\t\tUPPER\t\tSCALING\t\tDESCRIPTION");


    int ii;

    for(i=0;i< problem.nphases; i++)
    {

        	MatrixXd& control_scaling = problem.phase[i].scale.controls;
			MatrixXd& state_scaling   = problem.phase[i].scale.states;
			MatrixXd& param_scaling   = problem.phase[i].scale.parameters;
			double time_scaling      =  problem.phase[i].scale.time;
			int norder    = problem.phase[i].current_number_of_intervals;
			int nstates   = problem.phase[i].nstates;
        	int ncontrols = problem.phase[i].ncontrols;
        	int nparam    = problem.phase[i].nparameters;
        	int iphase_offset = get_iphase_offset(problem,i+1, workspace);
			int nvars_phase_i = get_nvars_phase_i(problem,i, workspace);
        	int offset1 = (norder+1)*ncontrols;
        	int offset2 = (norder+1)*(ncontrols+nstates);


	for (k=0; k<norder+1; k++) { // EIGEN_UPDATE

                for (ii=0;ii<ncontrols;ii++) {  //EIGEN_UPDATE
                        j = iphase_offset+(k)*ncontrols+ii;
			fprintf(outfile,"\n%i\t%e\t%e\t%e\t%e\tCONTROL %i NODE %i PHASE %i", j, xlb(j)/control_scaling(ii),
				  X(j)/control_scaling(ii), xub(j)/control_scaling(ii),  control_scaling(ii), ii, k, i+1);
			if ( xlb(j)-tol>X(j) || X(j)>xub(j)+tol ) fprintf(outfile," **BOUND VIOLATED");
   		        if ( (xlb(j)-tol<X(j) && X(j)<xlb(j)+tol) || (xub(j)-tol<X(j) && X(j)<xub(j)+tol) ) fprintf(outfile," **AT BOUND");

		}

	}

	for (k=0; k<norder+1; k++) {  // EIGEN_UPDATE


                for (ii=0;ii<nstates;ii++) {
                        j = iphase_offset+(k)*nstates+offset1+ii;
			fprintf(outfile,"\n%i\t%e\t%e\t%e\t%e\tSTATE %i NODE %i PHASE %i", j, xlb(j)/state_scaling(ii),
				X(j)/state_scaling(ii), xub(j)/state_scaling(ii), state_scaling(ii), ii, k, i+1);
			if ( xlb(j)-tol>X(j) || X(j)>xub(j)+tol ) fprintf(outfile," **BOUND VIOLATED");
   		        if ( (xlb(j)-tol<X(j) && X(j)<xlb(j)+tol) || (xub(j)-tol<X(j) && X(j)<xub(j)+tol) ) fprintf(outfile," **AT BOUND");

		}


	}

        for (ii=0;ii<nparam;ii++) {  // EIGEN_UPDATE
                        j = iphase_offset+offset2+ii;
			fprintf(outfile,"\n%i\t%e\t%e\t%e\t%e\tPARAMETER %i PHASE %i", j, xlb(j)/param_scaling(ii),
				X(j)/param_scaling(ii), xub(j)/param_scaling(ii), param_scaling(ii), ii, i+1);
			if ( xlb(j)-tol>X(j) || X(j)>xub(j)+tol ) fprintf(outfile," **BOUND VIOLATED");
   		        if ( (xlb(j)-tol<X(j) && X(j)<xlb(j)+tol) || (xub(j)-tol<X(j) && X(j)<xub(j)+tol) ) fprintf(outfile," **AT BOUND");

	}

        if (need_midpoint_controls(*workspace->algorithm, workspace)) {

		for (k=0; k<norder; k++) { // EIGEN_UPDATE

			for (ii=0;ii<ncontrols;ii++) {
				j = iphase_offset+offset2+nparam+(k)*ncontrols+ii;
				fprintf(outfile,"\n%i\t%e\t%e\t%e\t%e\tMIDPOINT CONTROL %i NODE %i PHASE %i", j,
					 xlb(j)/control_scaling(ii), X(j)/control_scaling(ii), xub(j)/control_scaling(ii),
					control_scaling(ii), ii, k, i+1);
				if ( xlb(j)-tol>X(j) || X(j)>xub(j)+tol ) fprintf(outfile," **BOUND VIOLATED");
   		                if ( (xlb(j)-tol<X(j) && X(j)<xlb(j)+tol) || (xub(j)-tol<X(j) && X(j)<xub(j)+tol) ) fprintf(outfile," **AT BOUND");

			}

		}
	}

        j = iphase_offset + nvars_phase_i-2;

	fprintf(outfile,"\n%i\t%e\t%e\t%e\t%e\tt0 PHASE %i", j, xlb(j)/time_scaling, X(j)/time_scaling,
		                                                xub(j)/time_scaling, time_scaling, i+1);
	if ( xlb(j)-tol>X(j) || X(j)>xub(j)+tol ) fprintf(outfile," **BOUND VIOLATED");
        if ( (xlb(j)-tol<X(j) && X(j)<xlb(j)+tol) || (xub(j)-tol<X(j) && X(j)<xub(j)+tol) ) fprintf(outfile," **AT BOUND");


        j = iphase_offset + nvars_phase_i-1;

	fprintf(outfile,"\n%i\t%e\t%e\t%e\t%e\ttf PHASE %i", j, xlb(j)/time_scaling, X(j)/time_scaling,
		                                                xub(j)/time_scaling, time_scaling, i+1);
	if ( xlb(j)-tol>X(j) || X(j)>xub(j)+tol ) fprintf(outfile," **BOUND VIOLATED");
        if ( (xlb(j)-tol<X(j) && X(j)<xlb(j)+tol) || (xub(j)-tol<X(j) && X(j)<xub(j)+tol) ) fprintf(outfile," **AT BOUND");



    }
    
    delete [] g_l;
    delete [] g_u;

}



