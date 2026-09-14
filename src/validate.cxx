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


#ifdef USE_SQP
// Implemented in qp_plugin_loader.cxx, and declared here as SQP_interface.cxx declares
// it: every QP backend lives in a plugin opened at run time, so whether this build can
// use one is a question for the loader and not for the compiler.
bool psopt_qp_plugin_available(const std::string& backend, std::string& message);
#endif

// The QP backends PSOPT knows how to name. One list, used both to check the option and
// to say what a build can actually load. The names were written out twice before -- once
// in the test and once in the message that reports it -- which is the arrangement that
// lets the two drift apart, and a list written twice is a list that will.
static const char* const qp_backend_names[] =
    { "GALAHAD", "ProxQP", "QPALM", "OSQP", "PIQP", "Clarabel", NULL };

static bool is_known_qp_backend(const string& name)
{
    for (const char* const* b = qp_backend_names; *b != NULL; ++b)
        if (name == *b) return true;
    return false;
}

static string qp_backend_name_list()
{
    string s;
    for (const char* const* b = qp_backend_names; *b != NULL; ++b) {
        if (!s.empty()) s += (b[1] == NULL) ? " and " : ", ";
        s += string("\"") + *b + "\"";
    }
    return s;
}


void validate_user_input(Prob& problem, Alg& algorithm, Workspace* workspace)
{
    int i;

    if (algorithm.nlp_method != "IPOPT" && algorithm.nlp_method != "SQP" )
       error_message("Incorrect NLP method specified. The valid values are \"IPOPT\" and \"SQP\". "
                     "SNOPT was supported until 2026 and has been removed; it is commercial, and "
                     "PSOPT's own SQP now fills the same place with no licence to obtain.");
    if (algorithm.collocation_method != "Legendre" && algorithm.collocation_method!="Chebyshev" && algorithm.collocation_method!="trapezoidal" && algorithm.collocation_method!="Hermite-Simpson" && algorithm.collocation_method!="Radau" && algorithm.collocation_method!="Gauss")
       error_message("Incorrect collocation method specified. Valid options are \"Legendre\" , \"Chebyshev\", \"trapezoidal\", \"Hermite-Simpson\", \"Radau\", and \"Gauss\" ");
    // The flexible mesh is a property of the integrated-residual transcription: the widths
    // are the widths OF something, and no other transcription here has elements to move.
    // Asking for it otherwise is a misunderstanding worth naming rather than ignoring,
    // because ir_flex_mesh_vars would simply return zero and the option would appear to work.
    //
    // Both local representations have them. The Nie-Kerrigan element of degree d spans d
    // intervals; the legacy cubic-Hermite element is one interval carrying one cubic, so a
    // flexible mesh there makes every interval width a variable. The second is a larger
    // decision vector for the same node count and a lower order per element, but it is the
    // same facility and it moves a node onto a switch the same way.
    if ( algorithm.ir_flexible_mesh && algorithm.transcription_method != "integrated-residual" )
       error_message("algorithm.ir_flexible_mesh requires algorithm.transcription_method = "
                     "\"integrated-residual\": the flexible mesh moves the element boundaries "
                     "of that transcription, and no other has element boundaries to move ");
    if ( algorithm.ir_flexible_mesh && algorithm.ir_local_order != 0 && algorithm.ir_local_order < 2 )
       error_message("algorithm.ir_local_order must be 0 (cubic Hermite) or at least 2 "
                     "(Nie-Kerrigan): no other local order is defined ");

    if ( flexible_partition_active(algorithm) &&
         ( min_partition_fraction(algorithm) <= 0.0 || min_partition_fraction(algorithm) >= 1.0 ) )
       error_message("the minimum element or segment fraction must lie strictly between 0 and 1: "
                     "it is the floor on a width as a fraction of the uniform width, and a piece "
                     "free to collapse to zero width has its ends coincident. See "
                     "algorithm.ir_min_element_fraction and algorithm.ms_min_segment_fraction ");

    // Automatic refinement of an element basis, and of a flexible mesh, was refused here
    // until the refinement could be written in the same currency as the transcription. It
    // now is: ir_refine_driver splits elements rather than inserting nodes, so the node
    // count moves in multiples of the stride by construction and the partition the flexible
    // mesh solved on is what the next mesh is built from. See ir_element_refinement_active.
    //
    // What has to be checked instead is that the mesh has somewhere to grow into, because
    // the workspace is sized once, from get_max_nodes, and the Betts ceiling it computes is
    // driven by mr_max_growth_factor and mr_max_iterations.
    if ( ir_element_refinement_active(algorithm) && algorithm.mr_max_growth_factor <= 0.0 )
       error_message("algorithm.mesh_refinement = \"automatic\" with the integrated-residual "
                     "transcription needs algorithm.mr_max_growth_factor > 0: it is the budget "
                     "the element refinement is allowed to spend ");

    if (algorithm.scaling != "automatic" && algorithm.scaling!="user")
       error_message("Incorrect scaling option specified. Valid options are \"automatic\" and \"user\" ");
    if (algorithm.transcription_method != "collocation"
        && algorithm.transcription_method != "integrated-residual"
        && algorithm.transcription_method != "multiple-shooting")
       error_message("Incorrect transcription_method specified. Valid options are \"collocation\", "
                     "\"integrated-residual\" and \"multiple-shooting\" ");

    if ( is_multiple_shooting(algorithm) ) {
       if ( algorithm.ms_steps_per_segment < 1 )
          error_message("algorithm.ms_steps_per_segment must be at least 1: it is the number of "
                        "fixed integrator steps taken across one shooting segment ");
       // Every one of these belongs to a transcription that builds a trajectory out of
       // decision variables, which this one does not: between the segment boundaries there is
       // an integrator and nothing to refine, regularise or bound the residual of. Saying so
       // is better than accepting the option and ignoring it.
       if ( algorithm.ms_control_parameterisation != "constant"
            && algorithm.ms_control_parameterisation != "linear"
            && algorithm.ms_control_parameterisation != "quadratic" )
          error_message("algorithm.ms_control_parameterisation must be \"constant\", "
                        "\"linear\" or \"quadratic\" ");
       if ( algorithm.ms_integrator != "RK4" && algorithm.ms_integrator != "RK8" )
          error_message("algorithm.ms_integrator must be \"RK4\" or \"RK8\": those are the "
                        "explicit schemes the segment integrator provides ");
       if ( algorithm.ms_path_samples < 0 )
          error_message("algorithm.ms_path_samples must be zero or positive: it is the number "
                        "of interior points per segment at which the path constraints are also "
                        "enforced ");
       // A new step count is of no use without another solve to use it in, and the outer
       // loop that provides one is the mesh-refinement loop. Refusing is better than
       // accepting the option and silently doing nothing with it, which is the failure mode
       // that made ps_method worth deleting.
       if ( algorithm.ms_adaptive_steps && algorithm.mesh_refinement != "automatic" )
          error_message("algorithm.ms_adaptive_steps needs algorithm.mesh_refinement = "
                        "\"automatic\": the step count is chosen BETWEEN solves -- it has to "
                        "be, since a step count that varied with the decision variables would "
                        "make the constraints non-smooth in them -- so there has to be another "
                        "solve for a new one to be used in ");
       if ( algorithm.ms_adaptive_steps && algorithm.ms_max_steps_per_segment < 1 )
          error_message("algorithm.ms_max_steps_per_segment must be at least 1: it is the "
                        "ceiling on any one segment's step count under ms_adaptive_steps ");
       if ( algorithm.ms_adaptive_steps && algorithm.ode_tolerance <= 0.0 )
          error_message("algorithm.ms_adaptive_steps needs algorithm.ode_tolerance > 0: it is "
                        "the quantity the step count is chosen to reach ");
       if ( algorithm.ms_algebraic_iterations < 1 )
          error_message("algorithm.ms_algebraic_iterations must be at least 1: it is the fixed, "
                        "unrolled number of Broyden iterations the half-explicit scheme spends "
                        "on a phase's algebraic equations at each stage ");

       // The algebraic declaration. Every check here is a sizing statement the user can get
       // wrong silently, and a half-explicit scheme that solves the wrong rows or reads the
       // wrong controls does not fail -- it returns a plausible answer to a different problem.
       for (i = 0; i < problem.nphases; i++) {
          const int nalg = problem.phase[i].nalgebraic;
          if ( nalg == 0 ) continue;
          if ( nalg < 0 )
             error_message("problem.phases(i).nalgebraic must not be negative: it is how many of "
                           "the phase's controls are the algebraic variables of a semi-explicit "
                           "index-1 DAE ");
          if ( nalg > problem.phase[i].ncontrols )
             error_message("problem.phases(i).nalgebraic exceeds the phase's ncontrols: the "
                           "algebraic variables are carried as the LAST nalgebraic controls, so "
                           "there have to be that many ");
          if ( nalg > problem.phase[i].npath )
             error_message("problem.phases(i).nalgebraic exceeds the phase's npath: the algebraic "
                           "equations are the FIRST nalgebraic path constraints, so there have to "
                           "be that many ");
          for (int j = 0; j < nalg; j++) {
             if ( (problem.phase[i].bounds.lower.path)(j) != (problem.phase[i].bounds.upper.path)(j) )
                error_message("the first nalgebraic path constraints of a phase that declares "
                              "nalgebraic must be EQUALITIES -- equal lower and upper bounds. "
                              "They are the algebraic equations of the DAE, solved for the "
                              "algebraic variables at every stage of the segment integrator; an "
                              "inequality has no solution to be solved for ");
          }
          // Half-explicit is explicit in the differential part and inherits its stability
          // restriction exactly. Saying so here is the difference between a user reaching for
          // this because their system is a DAE, which it serves, and reaching for it because
          // their system is stiff, which it does not.
          snprintf(workspace->text, sizeof(workspace->text),
             "\n>>> Note: phase %d declares %d algebraic variable(s); the segment integrator is "
             "\n>>> half-explicit and solves the first %d path constraint(s) for the last %d "
             "\n>>> control(s) at every stage, with %d unrolled Broyden iterations. The algebraic"
             "\n>>> relation then holds wherever those variables are defined, rather than only at"
             "\n>>> the segment boundaries, and the scheme keeps the order of ms_integrator."
             "\n>>> The differential part stays EXPLICIT, so this does nothing for a STIFF"
             "\n>>> system: it adds a class of problem, not a stability region.\n",
             i+1, nalg, nalg, nalg, algorithm.ms_algebraic_iterations);
          psopt_print(workspace, workspace->text);
       }
       // Equality path components are not sampled inside the segments, and a user who set
       // ms_path_samples expecting them to be should hear it from PSOPT rather than from a
       // return code. See ms_samplable_path_indices for why they cannot be.
       if ( algorithm.ms_path_samples > 0 ) {
          for (i = 0; i < problem.nphases; i++) {
             const int npath = problem.phase[i].npath;
             const int nsp   = ms_samplable_path_components(problem, i);
             if ( npath > 0 && nsp < npath ) {
                snprintf(workspace->text, sizeof(workspace->text),
                   "\n>>> Note: phase %d has %d path constraint(s) of which %d are equalities."
                   "\n>>> Equality components are imposed at the segment boundaries and are NOT"
                   "\n>>> sampled inside the segments: an equality demanded where the state is"
                   "\n>>> not a decision variable is an equation with nothing to answer it, and"
                   "\n>>> the NLP would be over-determined. Between the boundaries such a"
                   "\n>>> constraint holds only as well as the control parameterisation makes"
                   "\n>>> it hold -- exactly, if it involves the controls alone and the"
                   "\n>>> parameterisation is \"constant\"; otherwise to O(h).\n",
                   i+1, npath, npath - nsp);
                psopt_print(workspace, workspace->text);
             }
          }
       }
       if ( algorithm.ms_path_samples > algorithm.ms_steps_per_segment - 1 )
          error_message("algorithm.ms_path_samples must be at most "
                        "algorithm.ms_steps_per_segment - 1: the samples are placed at "
                        "integrator step boundaries, so there have to be step boundaries "
                        "inside the segment to place them at ");
       if ( algorithm.ir_regularization > 0.0 )
          error_message("algorithm.ir_regularization has no meaning with "
                        "transcription_method = \"multiple-shooting\": there is no discretised "
                        "residual to penalise, the segment integrator satisfies the dynamics "
                        "exactly for the scheme it uses ");
       if ( algorithm.ir_flexible_mesh )
          error_message("algorithm.ir_flexible_mesh belongs to the integrated-residual "
                        "transcription; multiple shooting has segments rather than elements, and "
                        "the option that moves them is algorithm.ms_flexible_segments ");
       if ( algorithm.ir_local_order != 0 )
          error_message("algorithm.ir_local_order belongs to the integrated-residual "
                        "transcription and has no meaning with multiple shooting ");
       if ( algorithm.ms_flexible_segments && algorithm.ms_path_samples == 0 )
          psopt_print(workspace,
                 "\n>>> Note: algorithm.ms_flexible_segments is on and ms_path_samples is zero, so "
                 "\n>>> the path constraints are imposed only at boundaries that are now free to "
                 "\n>>> move. Consider setting ms_path_samples.\n");
       // The two features are individually sound and interact badly exactly where both are
       // doing their job. A parabola through three values inside the control's box need not
       // stay inside it -- it overshoots by a quarter of the second difference -- and the
       // sharpest second difference a solution can present is a jump, which is precisely what
       // a moving boundary is there to put a node on. Measured on the minimum-time bang-bang
       // problem with bounds [-1,2]: with the boundary on the switch the parabola reaches
       // 2.375, nineteen per cent above its own upper bound, and that is the control the
       // segment integrator is handed. ms_path_samples does not help; it samples the PATH
       // constraints, and this is a variable bound.
       if ( algorithm.ms_flexible_segments && ms_quadratic_controls(algorithm) )
          psopt_print(workspace,
                 "\n>>> Note: ms_control_parameterisation = \"quadratic\" with "
                 "ms_flexible_segments on."
                 "\n>>> A parabola through three in-bounds control values can leave the bounds "
                 "between"
                 "\n>>> them, and it does so most where a moving boundary is most useful -- at a "
                 "corner"
                 "\n>>> of the optimal control. Use \"constant\" where the control rides its "
                 "bounds.\n");
       // Automatic segment refinement was refused here until the refinement could be
       // written in the currency this transcription is stated in. It now is: ms_refine_driver
       // refines on an indicator of the CONTROL PARAMETERISATION's error and of the PATH
       // CONSTRAINTS' coverage, which are the two things the segment count controls, and not
       // on the reported ODE error, which it does not. See ms_refine_driver.
       if ( ms_refinement_active(algorithm) && algorithm.ms_refine_tolerance <= 0.0 )
          error_message("algorithm.ms_refine_tolerance must be positive: it is the tolerance "
                        "automatic segment refinement compares its indicator against, and the "
                        "indicator is dimensionless ");
       if ( ms_refinement_active(algorithm) && algorithm.mr_max_growth_factor <= 0.0 )
          error_message("algorithm.mesh_refinement = \"automatic\" with multiple shooting "
                        "needs algorithm.mr_max_growth_factor > 0: it is the budget the "
                        "segment refinement is allowed to spend ");
       // diagnostic_level was refused here when this transcription had no costates. It has
       // them now (recover_costates_adjoint), and the part of the report that matters most --
       // the rank and conditioning of the constraint Jacobian -- never depended on the
       // transcription at all: it re-tapes the constraints at the final iterate and factorises
       // them, which is the same question whatever wrote the rows. solution_diagnostics says
       // which pieces do not apply here rather than being refused wholesale.
    }
    // nalgebraic is read by the segment integrator and by nothing else. Every other
    // transcription builds the trajectory out of decision variables, so an algebraic relation
    // there is an ordinary equality path constraint imposed at every node and needs no
    // declaration -- and accepting a number that changes nothing would be the worse outcome,
    // since the user would believe a DAE was being solved.
    if ( !is_multiple_shooting(algorithm) ) {
       for (i = 0; i < problem.nphases; i++)
          if ( problem.phase[i].nalgebraic > 0 )
             error_message("problem.phases(i).nalgebraic is read only by "
                           "transcription_method = \"multiple-shooting\", whose segment "
                           "integrator solves the algebraic equations at each of its own "
                           "stages. Under collocation or the integrated residual the "
                           "trajectory is made of decision variables and an algebraic relation "
                           "is an ordinary equality path constraint at every node, so leave "
                           "nalgebraic at zero ");
    }
    if (algorithm.transcription_method == "integrated-residual") {
       if (algorithm.collocation_method != "Hermite-Simpson")
          error_message("integrated-residual transcription currently requires collocation_method = \"Hermite-Simpson\" ");
       if (algorithm.ir_residual_nodes < 2)
          error_message("algorithm.ir_residual_nodes must be >= 2 for integrated-residual transcription ");
    }
    if (algorithm.ir_objective != "residual" && algorithm.ir_objective != "cost")
       error_message("Incorrect ir_objective specified. Valid options are \"residual\" and \"cost\" ");
    if (algorithm.ir_objective == "cost") {
       if (algorithm.transcription_method != "integrated-residual")
          error_message("ir_objective=\"cost\" requires transcription_method=\"integrated-residual\" (DAIR optimality step) ");
       if (algorithm.ir_regularization <= 0.0 && algorithm.ir_residual_bound < 0.0)
          error_message("ir_objective=\"cost\" needs the dynamics enforced: set ir_regularization>0 (penalty form) or ir_residual_bound>=0 (robust constraint form) ");
    }
    if (algorithm.ir_residual_bound >= 0.0) {
       if (algorithm.transcription_method != "integrated-residual" || algorithm.ir_objective != "cost")
          error_message("ir_residual_bound>=0 (robust-DAIR constraint form) requires transcription_method=\"integrated-residual\" and ir_objective=\"cost\" ");
    }
    if (algorithm.ir_dair) {
       if (algorithm.transcription_method != "integrated-residual")
          error_message("ir_dair=true requires transcription_method=\"integrated-residual\" ");
       if (algorithm.collocation_method != "Hermite-Simpson")
          error_message("ir_dair=true requires collocation_method=\"Hermite-Simpson\" ");
       if (algorithm.ir_dair_delta_factor <= 0.0)
          error_message("algorithm.ir_dair_delta_factor must be > 0 ");
    }
    if ( algorithm.objective_form != "as-posed" && algorithm.objective_form != "mayer"
         && algorithm.objective_form != "auto" )
       error_message("algorithm.objective_form must be \"as-posed\", \"mayer\" or \"auto\" ");

    if (algorithm.ir_local_order != 0) {
       if (algorithm.ir_local_order < 2)
          error_message("algorithm.ir_local_order must be 0 (legacy cubic-Hermite IR) or >= 2 (Nie-Kerrigan) ");
       if (algorithm.transcription_method != "integrated-residual")
          error_message("ir_local_order requires transcription_method=\"integrated-residual\" ");
       // The residual box constrains the residual only where it samples it. A degree-d element
       // state has d+1 coefficients per component and its residual xdot-f is of higher degree
       // still, so a rule with too few points leaves the residual free to oscillate between
       // them: the box is then satisfied to its stated tolerance while the true error is orders
       // of magnitude larger. Measured on examples/dae_i3 at d=4, delta=1e-6, and with the
       // holonomic constraint left as a pointwise path constraint rather than folded into the
       // residual -- that is, `./dae_i3 4 1e-6 out` -- the maximum relative local error is
       // 9.7e-2 with m=4, 7.9e-4 with m=5 and 1.5e-8, the box itself, from m=6 onward. Hence
       // m >= d+2, which the former m >= d did not give.
       //
       // The configuration is named because it matters, and because leaving it out has already
       // sent one reader looking for a regression that was not there. Folded in, the same d and
       // delta do not converge at m = 6, 7, 8, 10 or 12: the box cannot be met by a pendulum of
       // admissible length, which is a property of that problem rather than of m and is recorded
       // in the example. And the two figures below the rule can no longer be obtained through this
       // interface, since the rule refuses them; they are the measurement the rule exists
       // because of, not one a reader can repeat.
       if (algorithm.ir_residual_nodes < algorithm.ir_local_order + 2)
          error_message("ir_residual_nodes must be >= ir_local_order+2, so that the residual box "
                        "samples the element residual densely enough to constrain it between the "
                        "sample points ");
    }
    if (algorithm.ir_regularization < 0.0)
       error_message("algorithm.ir_regularization must be >= 0 ");
    if (algorithm.ir_regularization > 0.0) {
       if (algorithm.collocation_method != "Hermite-Simpson")
          error_message("integrated-residual regularization (ir_regularization>0) currently requires collocation_method = \"Hermite-Simpson\" ");
       if (algorithm.ir_residual_nodes < 2)
          error_message("algorithm.ir_residual_nodes must be >= 2 when ir_regularization>0 ");
    }
    if (algorithm.defect_scaling != "state-based" && algorithm.defect_scaling!="jacobian-based")
       error_message("Incorrect differential defect scaling option specified. Valid options are \"state-based\" and \"jacobian-based\" ");
    if (algorithm.derivatives != "automatic" && algorithm.derivatives!="numerical")
       error_message("Incorrect derivatives option specified. Valid options are \"automatic\" and \"numerical\" ");
    if (algorithm.hessian != "exact" && algorithm.hessian!="limited-memory" && algorithm.hessian!="numerical")
       error_message("Incorrect algorithm.hessian option specified. Valid options are \"limited-memory\", \"exact\" and \"numerical\" ");
    if (!is_known_qp_backend(algorithm.qp_solver))
       error_message(("Incorrect algorithm.qp_solver option specified. Valid options are "
                      + qp_backend_name_list() + " ").c_str());
#ifdef USE_SQP
    // A name on that list is not a backend this build can load. Every backend is a
    // plugin opened at run time, so the two questions are independent and only the
    // second one matters to a solve. Asked here, the absence is a set-up failure with
    // the missing name in it and a list of what this build does carry. Left to the
    // first subproblem -- which is where it was discovered until now -- it is a
    // subproblem the SQP declines: the solver stops at iteration zero and reports
    // status 2 through solution.nlp_return_code, while solution.error_flag stays at
    // zero, because nothing was thrown and the set-up really was sound. A caller
    // reading error_flag alone is then told the run succeeded and is handed back the
    // cost of its own initial guess. The comment on psopt_qp_plugin_available has said
    // since it was written that validate() asks this question; until now it did not.
    if (algorithm.nlp_method == "SQP") {
       string why;
       if (!psopt_qp_plugin_available(algorithm.qp_solver, why)) {
          string carried;
          for (const char* const* b = qp_backend_names; *b != NULL; ++b) {
             string ignored;
             if (psopt_qp_plugin_available(*b, ignored)) {
                if (!carried.empty()) carried += ", ";
                carried += *b;
             }
          }
          error_message(("algorithm.qp_solver = \"" + algorithm.qp_solver + "\" is a valid "
                         "name, but this build cannot load that backend: " + why
                       + " This build can load: "
                       + (carried.empty() ? string("no QP backend at all") : carried) + ". ").c_str());
       }
    }
#endif
    // Zero asks for the automatic budget, which scales with the subproblem; anything
    // else is used as given, and a handful of iterations is not a budget.
    if (algorithm.qp_iter_max != 0 && algorithm.qp_iter_max < 10)
       error_message("algorithm.qp_iter_max is too small; it must be at least 10, or 0 for the automatic budget ");
    if (algorithm.trust_region_radius < 0.0)
       error_message("algorithm.trust_region_radius must be positive, or 0 for the default ");
    if (algorithm.trust_region != "box" && algorithm.trust_region != "l2")
       error_message("Incorrect algorithm.trust_region option specified. Valid options are \"box\" and \"l2\" ");
    // The trust region exists to make an indefinite model usable. The quasi-Newton model
    // is built positive definite and is given no region at all, so asking for a Euclidean
    // one there asks for nothing; said here rather than left to be inferred from a run
    // that behaves exactly as it did before.
    if (algorithm.trust_region == "l2" && algorithm.nlp_method == "SQP"
                                       && algorithm.hessian != "exact") {
       snprintf(workspace->text,sizeof(workspace->text),
                "\n*** Warning: algorithm.trust_region = \"l2\" applies only to hessian = \"exact\"; "
                "the quasi-Newton model is given no trust region");
       psopt_print(workspace,workspace->text);
    }
    if (algorithm.sqp_strategy != "M" && algorithm.sqp_strategy != "FM"
                                     && algorithm.sqp_strategy != "F")
       error_message("Incorrect algorithm.sqp_strategy option specified. Valid options are \"M\", \"FM\" and \"F\" ");
    if (algorithm.qp_restoration != "elastic" && algorithm.qp_restoration != "relaxation")
       error_message("Incorrect algorithm.qp_restoration option specified. Valid options are \"elastic\" and \"relaxation\" ");
    if (algorithm.elastic_penalty != "weights" && algorithm.elastic_penalty != "multipliers")
       error_message("Incorrect algorithm.elastic_penalty option specified. Valid options are \"weights\" and \"multipliers\" ");
    if (algorithm.qp_solver != "GALAHAD" && algorithm.nlp_method != "SQP") {
       snprintf(workspace->text,sizeof(workspace->text),"\n*** Warning: algorithm.qp_solver applies only to nlp_method = \"SQP\"");
       psopt_print(workspace,workspace->text);
    }
    if (algorithm.on_error != "fail-fast" && algorithm.on_error != "fail-soft")
       error_message("Incorrect algorithm.on_error option specified. Valid options are \"fail-fast\" and \"fail-soft\" ");
    if (algorithm.hessian == "numerical" && algorithm.nlp_method !="IPOPT") {
       snprintf(workspace->text,sizeof(workspace->text),"\n*** Warning: the 'numerical' algorithm.hessian option is only available with the IPOPT solver");
       psopt_print(workspace,workspace->text);
    }
    if (algorithm.hessian == "exact" && algorithm.nlp_method !="IPOPT" && algorithm.nlp_method !="SQP") {
       snprintf(workspace->text,sizeof(workspace->text),"\n*** Warning: the 'exact' algorithm.hessian option is only available with the IPOPT and SQP solvers");
       psopt_print(workspace,workspace->text);
    }
    if (algorithm.diff_matrix != "standard" && algorithm.diff_matrix != "reduced-roundoff" )
       error_message("Incorrect algorithm.diff_matrix option specified. Valid options are \"standard\" and \"reduced-roundoff\" ");


    if (algorithm.hessian == "exact" && algorithm.derivatives !="automatic") {
       snprintf(workspace->text,sizeof(workspace->text),"\n*** Warning: the 'exact' algorithm.hessian option is only available with automatic derivatives");
       psopt_print(workspace,workspace->text);
    }
    if (algorithm.nlp_tolerance <= 0)
       error_message("algorithm.nlp_tolerance must be positive");
    if (algorithm.nlp_iter_max <= 0)
       error_message("algorithm.iter_max must be positive");

    if (algorithm.nsteps_error_integration <= 0)
       error_message("algorithm.nsteps_error_integration must be positive");

    if (algorithm.ode_tolerance <= 0)
       error_message("algorithm.ode_tolerance must be positive");

    if (algorithm.mr_max_growth_factor <= 0 || algorithm.mr_max_growth_factor > 1 )
       error_message("algorithm.mr_max_growth_factor must be in the interval (0,1]");
    if (algorithm.mr_max_iterations <= 0)
       error_message("algorithm.mr_max_iterations must be positive");
    if (algorithm.mr_min_order < 2)
       error_message("algorithm.mr_min_order must be >= 2");
    if (algorithm.mr_max_order < algorithm.mr_min_order)
       error_message("algorithm.mr_max_order must be >= algorithm.mr_min_order");

    if (algorithm.mr_kappa <= 0 || algorithm.mr_kappa>1.0 )
       error_message("algorithm.mr_kappa must be in the interval (0,1]");

    if (algorithm.mr_M1 <= 0  )
       error_message("algorithm.mr_M1 must be positive");

    if (algorithm.mr_switch_detection != 0 && algorithm.mr_switch_detection != 1 )
       error_message("algorithm.mr_switch_detection must be 0 or 1");

    if (algorithm.switch_order < 0  )
       error_message("algorithm.switch_order must be >= 0");

    if (algorithm.mesh_refinement != "automatic" &&  algorithm.mesh_refinement != "manual" )
       error_message("algorithm.mesh_refinement must either \"manual\" or \"automatic\" ");



    for (i=0;i<problem.nphases;i++)
    {

         if (problem.phase[i].ncontrols>1) {
		      (problem.phase[i].bounds.lower.controls).resize(problem.phase[i].ncontrols,1);
         	(problem.phase[i].bounds.upper.controls).resize(problem.phase[i].ncontrols,1);
         }

         if (problem.phase[i].nstates>1) {
		      (problem.phase[i].bounds.lower.states).resize(problem.phase[i].nstates,1);
         	(problem.phase[i].bounds.upper.states).resize(problem.phase[i].nstates,1);
         }

         if (problem.phase[i].nevents>1) {
	       	(problem.phase[i].bounds.lower.events).resize(problem.phase[i].nevents,1); 
         	(problem.phase[i].bounds.upper.events).resize(problem.phase[i].nevents,1);
         }

         if (problem.phase[i].nparameters>1) {
		      (problem.phase[i].bounds.lower.parameters).resize(problem.phase[i].nparameters,1);
         	(problem.phase[i].bounds.upper.parameters).resize(problem.phase[i].nparameters,1);
         }



         if (problem.phase[i].ncontrols>0 && (problem.phase[i].bounds.lower.controls.array() >  problem.phase[i].bounds.upper.controls.array() ).any() )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible control variable bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if (problem.phase[i].nstates>0 && ( problem.phase[i].bounds.lower.states.array() >  problem.phase[i].bounds.upper.states.array() ).any() )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible state variable bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if (problem.phase[i].nevents>0 && ( problem.phase[i].bounds.lower.events.array() >  problem.phase[i].bounds.upper.events.array() ).any() )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible event bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if (problem.phase[i].nparameters>0 && ( problem.phase[i].bounds.lower.parameters.array() >  problem.phase[i].bounds.upper.parameters.array() ).any() )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible static parameter bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if ( problem.phase[i].bounds.lower.StartTime >  problem.phase[i].bounds.upper.StartTime  )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible start time bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         if ( problem.phase[i].bounds.lower.EndTime >  problem.phase[i].bounds.upper.EndTime  )
         {
                snprintf(workspace->text,sizeof(workspace->text),"Infeasible end time bounds supplied by the user in phase %i",i);
 		error_message(workspace->text);
         }

         // hp-adaptive fixed mesh (Route B, increment 1). When a phase carries an
         // explicit multi-interval mesh via hp_orders / hp_breakpoints, validate it:
         // Radau-only for now; K orders and K-1 breaks; breaks strictly increasing in
         // (0,1); each interval order >= 2 (a single-point interval is degenerate).
         if ( hp_mesh_active(problem.phase[i]) )
         {
            if ( algorithm.collocation_method != "Radau" && algorithm.collocation_method != "Gauss" && algorithm.collocation_method != "Legendre" && algorithm.collocation_method != "Chebyshev" )
            {
               snprintf(workspace->text,sizeof(workspace->text),"hp-adaptive mesh (hp_orders) in phase %i currently requires collocation_method = \"Radau\", \"Gauss\", \"Legendre\" or \"Chebyshev\"",i+1);
               error_message(workspace->text);
            }
            int Khp = (int) problem.phase[i].hp_orders.size();
            if ( (int) problem.phase[i].hp_breakpoints.size() != Khp-1 )
            {
               snprintf(workspace->text,sizeof(workspace->text),"In phase %i, length(hp_breakpoints) must equal length(hp_orders)-1",i+1);
               error_message(workspace->text);
            }
            for (int kk=0; kk<Khp; kk++)
            {
               if ( problem.phase[i].hp_orders(kk) < 2 )
               {
                  snprintf(workspace->text,sizeof(workspace->text),"In phase %i, every hp_orders entry must be >= 2",i+1);
                  error_message(workspace->text);
               }
            }
            for (int kk=0; kk<Khp-1; kk++)
            {
               double bk = problem.phase[i].hp_breakpoints(kk);
               if ( bk <= 0.0 || bk >= 1.0 )
               {
                  snprintf(workspace->text,sizeof(workspace->text),"In phase %i, hp_breakpoints must lie strictly in the open interval (0,1)",i+1);
                  error_message(workspace->text);
               }
               if ( kk>0 && bk <= problem.phase[i].hp_breakpoints(kk-1) )
               {
                  snprintf(workspace->text,sizeof(workspace->text),"In phase %i, hp_breakpoints must be strictly increasing",i+1);
                  error_message(workspace->text);
               }
            }
         }

	 if ( problem.phase[i].nobserved >  0  )
         {

//	        if ( problem.phase[i].observation_nodes(1) != problem.phase[i].bounds.lower.StartTime && problem.phase[i].observation_nodes(1) != problem.phase[i].bounds.upper.StartTime )
//		{
//                 snprintf(workspace->text,sizeof(workspace->text),"Initial observation time must be equal to start time in phase %i",i);
// 		  error_message(workspace->text);
//		}

//	        if ( fabs( problem.phase[i].observation_nodes("end") - problem.phase[i].bounds.lower.EndTime ) > 0.0001 && fabs(problem.phase[i].observation_nodes("end") - problem.phase[i].bounds.upper.EndTime )>0.0001 )
//		{
//		  snprintf(workspace->text,sizeof(workspace->text),"Final observation time must be equal to end time in phase %i",i+1);
// 		  error_message(workspace->text);
//		}

	        if ( problem.phase[i].nsamples != length( problem.phase[i].observation_nodes)  )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"Length of observation nodes vector in phase %i must be equal to problem.phases(%i).nsamples",i+1, i+1);
 		  error_message(workspace->text);
		}

		if (  isEmpty( problem.phase[i].residual_weights ) ) {
                    problem.phase[i].residual_weights = ones( problem.phase[i].nobserved, problem.phase[i].nsamples );
		}

		if ( problem.phase[i].nsamples !=  problem.phase[i].residual_weights.cols() )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"The number of columns of the residual weight vector in phase %i must be equal to problem.phases(%i).nsamples",i+1, i+1);
 		  error_message(workspace->text);
		}

		if ( problem.phase[i].nobserved !=  problem.phase[i].residual_weights.rows() )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"The number of rows of the residual weight vector in phase %i must be equal to the number of observed variables", i+1 );
 		  error_message(workspace->text);
		}

		if (  isEmpty(problem.phase[i].covariance)  ) {
                    problem.phase[i].covariance = eye( problem.phase[i].nobserved );
		}

		if ( problem.phase[i].nobserved !=  problem.phase[i].covariance.rows() && problem.phase[i].nobserved !=  problem.phase[i].covariance.cols()  )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"The number of rows and columns of matrix problem.phases(%i).covariance must be equal to problem.phases(%i).nobserved",i+1, i+1);
 		  error_message(workspace->text);
		}

		if ( !isSymmetric(problem.phase[i].covariance)  )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"Matrix problem.phases(%i).covariance must be symmetric",i+1);
 		  error_message(workspace->text);
		}

		if ( problem.phase[i].regularization_factor< 0  )
		{
		  snprintf(workspace->text,sizeof(workspace->text),"problem.phases(%i).regularization_factor must be positive",i+1);
 		  error_message(workspace->text);
		}

		problem.integrand_cost 	= NULL;
                problem.endpoint_cost 	= &endpoint_cost_for_parameter_estimation;
                problem.phase[i].zero_cost_integrand = true;



         }

    }

   if (problem.nlinkages>0 && ( problem.bounds.lower.linkage.array() >  problem.bounds.upper.linkage.array() ).any() )
   {
         snprintf(workspace->text,sizeof(workspace->text),"Infeasible phase linkage bounds supplied by the user");
         error_message(workspace->text);
   }

   if ( length(problem.bounds.lower.times) !=  length(problem.bounds.upper.times) || (!isEmpty(problem.bounds.lower.times) && length(problem.bounds.lower.times)!=problem.nphases+1) )
   {
         snprintf(workspace->text,sizeof(workspace->text),"Incorrect length of problem.bounds.lower.times or problem.bounds.upper.times");
         error_message(workspace->text);
   }

   if ( !isEmpty(problem.bounds.lower.times) ) {
     for (i=0;i<problem.nphases;i++) { //EIGEN_UPDATE
	    problem.phase[i].bounds.lower.StartTime = problem.bounds.lower.times(i);
	    problem.phase[i].bounds.upper.StartTime = problem.bounds.upper.times(i);
	    problem.phase[i].bounds.lower.EndTime   = problem.bounds.lower.times(i+1);
	    problem.phase[i].bounds.upper.EndTime   = problem.bounds.upper.times(i+1);
     }
   }

}
