/*********************************************************************************************

This file is part of the PSOPT library, a software tool for computational optimal control

Copyright (C) 2009-2020 Victor M. Becerra

This library is free software; you can redistribute it and/or modify it under the terms of the
GNU Lesser General Public License as published by the Free Software Foundation; either version
2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this
library; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
Boston, MA 02111-1307 USA

Author:    Professor Victor M. Becerra
Address:   University of Portsmouth
           School of Energy and Electronic Engineering
           Portsmouth PO1 3DJ
           United Kingdom
e-mail:    v.m.becerra@ieee.org

**********************************************************************************************/

// A DEFECT ROW A TRANSCRIPTION CANNOT FILL IS NOT AN EQUALITY CONSTRAINT
//
// Every transcription reuses the collocation layout, whose defect block holds
// nstates*(norder+1) rows. Only the differentiation-matrix schemes fill all of it: multiple
// shooting has norder matching conditions, trapezoidal and Hermite-Simpson have norder
// intervals, Radau does not collocate its terminal node, Gauss collocates no breakpoint. Each
// of those writes 0.0 into the rows it cannot fill.
//
// A row of zeros with bounds [0,0] is an equality constraint on nothing, and IPOPT counts the
// equality constraints against the variables and refuses outright when there are more of the
// former: "Too few degrees of freedom (n_x, n_c)". The count is nstates per phase, so the
// refusal is invisible on an ordinary problem and decisive on one whose state has been
// REPLICATED -- which is what a scenario-augmented robust design is, M copies of the state
// against one copy of the control.
//
// The tests below replicate a plant K times, driven by one shared control from one shared
// initial condition, so that every copy is the same trajectory and the objective cannot depend
// on K. They are sized so that the phantom count refuses them: at nine nodes and five copies
// of a two-state plant the count is 9 - 10 = -1, while the problem has nine degrees of freedom.
// Before NLP_bounds freed the padded rows these returned -10 and nothing else.

#include "gtest/gtest.h"
#include <psopt.h>

#include <cmath>
#include <string>
#include <vector>

namespace pad {

static int K_COPIES = 1;          // how many copies of the plant

// Minimum effort to bring the first copy's position to one, with the end state free so that
// the only equality events are the shared initial conditions. Copies 2..K appear in no cost
// and in no event: they are driven by the same control from the same state, so they are the
// same trajectory, and a transcription that treats the replication correctly must return an
// objective that does not depend on K at all.
adouble endpoint_cost(adouble* i, adouble* f, adouble*, adouble&, adouble&, adouble*,
                      int, Workspace*)
{
    return 10.0*(f[0] - 1.0)*(f[0] - 1.0);
}

adouble integrand_cost(adouble*, adouble* c, adouble*, adouble&, adouble*, int, Workspace*)
{
    return c[0]*c[0];
}

void dae(adouble* d, adouble*, adouble* s, adouble* c, adouble*, adouble&,
         adouble*, int, Workspace*)
{
    for (int k = 0; k < K_COPIES; k++) {
        d[2*k]     = s[2*k + 1];
        d[2*k + 1] = c[0];
    }
}

void events(adouble* e, adouble* i, adouble*, adouble*, adouble&, adouble&,
            adouble*, int, Workspace*)
{
    for (int k = 0; k < 2*K_COPIES; k++) e[k] = i[k];
}

void linkages(adouble*, adouble*, Workspace*) {}

struct Run { int flag; int rc; double J; MatrixXd x; };

static Run solve(int copies, int nodes, const std::string& scheme,
                 const std::string& transcription, bool free_padded = true)
{
    K_COPIES = copies;
    const int nst = 2*copies;

    Alg algorithm; Sol solution; Prob problem;
    Run out; out.flag = -1; out.rc = -99; out.J = 0.0;

    problem.name        = "padded defect rows, replicated state";
    problem.outfilename = "test_defect_padding.txt";
    problem.nphases     = 1;
    problem.nlinkages   = 0;
    psopt_level1_setup(problem);

    problem.phases(1).nstates   = nst;
    problem.phases(1).ncontrols = 1;
    problem.phases(1).nevents   = nst;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     << nodes;
    psopt_level2_setup(problem, algorithm);

    for (int j = 0; j < nst; j++) {
        problem.phases(1).bounds.lower.states(j) = -10.0;
        problem.phases(1).bounds.upper.states(j) =  10.0;
        problem.phases(1).bounds.lower.events(j) =  0.0;
        problem.phases(1).bounds.upper.events(j) =  0.0;
    }
    problem.phases(1).bounds.lower.controls << -10.0;
    problem.phases(1).bounds.upper.controls <<  10.0;
    problem.phases(1).bounds.lower.StartTime = 0.0; problem.phases(1).bounds.upper.StartTime = 0.0;
    problem.phases(1).bounds.lower.EndTime   = 1.0; problem.phases(1).bounds.upper.EndTime   = 1.0;

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost  = &endpoint_cost;
    problem.dae            = &dae;
    problem.events         = &events;
    problem.linkages       = &linkages;

    problem.phases(1).guess.states   = zeros(nst, nodes);
    problem.phases(1).guess.controls = zeros(1, nodes);
    problem.phases(1).guess.time     = linspace(0.0, 1.0, nodes);

    algorithm.nlp_method            = "IPOPT";
    algorithm.scaling               = "automatic";
    algorithm.derivatives           = "automatic";
    algorithm.nlp_iter_max          = 2000;
    algorithm.nlp_tolerance         = 1.0e-10;
    algorithm.print_level           = 0;
    algorithm.mesh_refinement       = "manual";
    algorithm.collocation_method    = scheme;
    algorithm.free_padded_defect_rows = free_padded;
    if (!transcription.empty()) {
        algorithm.transcription_method        = transcription;
        algorithm.ms_integrator               = "RK4";
        algorithm.ms_steps_per_segment        = 8;
        algorithm.ms_control_parameterisation = "linear";
    }

    out.flag = psopt(solution, problem, algorithm);
    out.rc   = solution.nlp_return_code;
    if (out.flag == 0) {
        out.J = solution.cost;
        out.x = solution.get_states_in_phase(1);
    }
    return out;
}

// The largest difference between any copy's trajectory and the first copy's. They are driven
// by the same control from the same state, so this is round-off or the replication is wrong.
static double copy_spread(const Run& r, int copies)
{
    double worst = 0.0;
    for (int k = 1; k < copies; k++)
        for (int j = 0; j < 2; j++)
            for (int c = 0; c < r.x.cols(); c++)
                worst = std::max(worst, std::fabs(r.x(2*k + j, c) - r.x(j, c)));
    return worst;
}

} // namespace pad


// Ten stored nodes (nine intervals) and eight copies of a two-state plant, one control.
// Variables 8*2*10 + 10 = 170. Genuine equalities: 8*2*9 trapezoid defects + 8*2 initial
// conditions = 160, so the problem has ten degrees of freedom. Counted equalities: the
// defect block holds 8*2*10 rows whatever the scheme fills, so 176 -- six more than there
// are variables, and IPOPT refuses. Both halves are asserted, because the option's whole
// content is the difference between them.
TEST(DefectPadding, AReplicatedStateSolvesOnlyWhenThePaddedRowsAreFreed)
{
    const pad::Run one  = pad::solve(1, 9, "trapezoidal", "");
    ASSERT_EQ(one.flag, 0) << "the unreplicated problem, IPOPT return code " << one.rc;

    const pad::Run counted = pad::solve(8, 9, "trapezoidal", "", false);
    EXPECT_EQ(counted.rc, -10)
        << "with the padded rows counted as equalities IPOPT should refuse this for "
        << "degrees of freedom it has: 176 counted equalities against 170 variables, "
        << "where 160 of the equalities hold any dynamics at all. "
        << "If this ever stops being -10, the count has changed and the option's "
        << "justification wants re-reading.";

    const pad::Run five = pad::solve(8, 9, "trapezoidal", "", true);
    ASSERT_EQ(five.flag, 0)
        << "eight copies at ten nodes with the rows freed, IPOPT return code "
        << five.rc;
    EXPECT_EQ(five.rc, 0);

    EXPECT_NEAR(five.J, one.J, 1.0e-9)
        << "the copies are the same trajectory and only the first one is costed, so "
        << "replication cannot change the objective";
    EXPECT_LT(pad::copy_spread(five, 8), 1.0e-9)
        << "every copy is driven by the same control from the same state";
}


// And the option leaves alone the schemes that pad nothing. Legendre holds norder+1
// conditions per state against norder+1 stored values, one of which the initial condition
// pins, so its deficit on a replicated state is real: it is refused either way, and that
// is the one place the intuition "collocation cannot carry the scenarios" was right.
TEST(DefectPadding, GlobalLobattoCollocationIsRefusedEitherWay)
{
    const pad::Run freed   = pad::solve(6, 9, "Legendre", "", true);
    const pad::Run counted = pad::solve(6, 9, "Legendre", "", false);
    EXPECT_EQ(freed.rc, -10)   << "Legendre pads nothing, so freeing nothing changes";
    EXPECT_EQ(counted.rc, -10);
}


// The same thing under every transcription that pads, so that a scheme cannot quietly lose
// the property. Legendre is included as the control: it pads nothing, so it is unaffected by
// the change and is here to show the test is not measuring the fix twice.
TEST(DefectPadding, EveryPaddingSchemeCarriesAReplicatedState)
{
    struct Case { const char* scheme; const char* transcription; };
    const std::vector<Case> cases = {
        {"trapezoidal",     ""},
        {"Hermite-Simpson", ""},
        {"Radau",           ""},
        {"Hermite-Simpson", "multiple-shooting"},
    };
    for (const Case& c : cases) {
        const pad::Run one = pad::solve(1, 11, c.scheme, c.transcription, true);
        ASSERT_EQ(one.flag, 0) << c.scheme << " " << c.transcription
                               << ", one copy, IPOPT return code " << one.rc;
        const pad::Run six = pad::solve(6, 11, c.scheme, c.transcription, true);
        ASSERT_EQ(six.flag, 0) << c.scheme << " " << c.transcription
                               << ", six copies, IPOPT return code " << six.rc;
        EXPECT_NEAR(six.J, one.J, 1.0e-7) << c.scheme << " " << c.transcription;
        EXPECT_LT(pad::copy_spread(six, 6), 1.0e-7) << c.scheme << " " << c.transcription;
    }
}
