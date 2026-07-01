#!/usr/bin/env python3
"""
psopt_codegen.py -- reusable CasADi-trace -> adouble-templated-C++ emitter for the
PSOPT hybrid wrapper (faithful route). Walks a CasADi SX Function's instruction
graph and emits straight-line scalar-templated C++ that plugs into PSOPT's fixed
dae/cost/events signatures, so CppAD tapes it exactly as hand-written code.
"""
import casadi as ca

BINARY = {
    ca.OP_ADD: "{a} + {b}", ca.OP_SUB: "{a} - {b}",
    ca.OP_MUL: "{a} * {b}", ca.OP_DIV: "{a} / {b}",
    ca.OP_POW: "pow({a}, {b})", ca.OP_CONSTPOW: "pow({a}, {b})",
    ca.OP_FMIN: "fmin({a}, {b})", ca.OP_FMAX: "fmax({a}, {b})",
    ca.OP_ATAN2: "atan2({a}, {b})",
}
UNARY = {
    ca.OP_NEG: "-({a})", ca.OP_EXP: "exp({a})", ca.OP_LOG: "log({a})",
    ca.OP_SQRT: "sqrt({a})", ca.OP_SQ: "({a}) * ({a})",
    ca.OP_SIN: "sin({a})", ca.OP_COS: "cos({a})", ca.OP_TAN: "tan({a})",
    ca.OP_ASIN: "asin({a})", ca.OP_ACOS: "acos({a})", ca.OP_ATAN: "atan({a})",
    ca.OP_TANH: "tanh({a})", ca.OP_SINH: "sinh({a})", ca.OP_COSH: "cosh({a})",
    ca.OP_FABS: "fabs({a})", ca.OP_INV: "1.0 / ({a})", ca.OP_ASSIGN: "{a}",
}
STD_FUNCS = {ca.OP_EXP:"exp",ca.OP_LOG:"log",ca.OP_SQRT:"sqrt",ca.OP_SIN:"sin",
    ca.OP_COS:"cos",ca.OP_TAN:"tan",ca.OP_ASIN:"asin",ca.OP_ACOS:"acos",
    ca.OP_ATAN:"atan",ca.OP_TANH:"tanh",ca.OP_SINH:"sinh",ca.OP_COSH:"cosh",
    ca.OP_FABS:"fabs",ca.OP_POW:"pow",ca.OP_CONSTPOW:"pow",ca.OP_FMIN:"fmin",
    ca.OP_FMAX:"fmax",ca.OP_ATAN2:"atan2"}

def emit_instructions(f, input_names, out_writer, indent="    "):
    lines, used = [], set()
    lines.append(f"{indent}T w[{max(f.sz_w(),1)}];")
    for k in range(f.n_instructions()):
        op = f.instruction_id(k)
        ii = list(f.instruction_input(k)); oo = list(f.instruction_output(k))
        if op == ca.OP_INPUT:
            lines.append(f"{indent}w[{oo[0]}] = {input_names[ii[0]]}[{ii[1]}];")
        elif op == ca.OP_OUTPUT:
            lines.append(f"{indent}{out_writer(oo[0], oo[1])} = w[{ii[0]}];")
        elif op == ca.OP_CONST:
            lines.append(f"{indent}w[{oo[0]}] = T({float(f.instruction_constant(k))!r});")
        elif op in UNARY:
            if op in STD_FUNCS: used.add(STD_FUNCS[op])
            lines.append(f"{indent}w[{oo[0]}] = {UNARY[op].format(a=f'w[{ii[0]}]')};")
        elif op in BINARY:
            if op in STD_FUNCS: used.add(STD_FUNCS[op])
            lines.append(f"{indent}w[{oo[0]}] = {BINARY[op].format(a=f'w[{ii[0]}]', b=f'w[{ii[1]}]')};")
        else:
            raise NotImplementedError(f"CasADi opcode {op} not mapped")
    return lines, used

def _using(used, indent="    "):
    return "".join(f"{indent}using std::{fn};\n" for fn in sorted(used))

def emit_header(name, dims, dae_f, L_f, phi_f, ev_f, obs_f=None):
    """dims = dict(nx, nu, npar, npath, nevents[, nobs]). Returns the header text.
    If obs_f is given, also emits gen_observation (parameter-estimation problems)."""
    nx, npath, nev = dims["nx"], dims["npath"], dims["nevents"]
    G = f"GEN_{name.upper()}_HPP"
    H = [f"// {name}_generated.hpp -- AUTO-GENERATED from a CasADi trace. Do not edit.",
         f"#ifndef {G}", f"#define {G}", "#include <cmath>", ""]

    b, u = emit_instructions(dae_f, ["states","controls","parameters","time"],
        lambda a,e: f"derivatives[{e}]" if a==0 else f"path[{e}]")
    H += ["template<class T> inline void gen_dae(T* derivatives, T* path,",
          "        const T* states, const T* controls, const T* parameters, const T* time) {",
          _using(u) + f"    for(int _i=0;_i<{nx};_i++) derivatives[_i]=T(0.0);"
          + (f"\n    for(int _i=0;_i<{npath};_i++) path[_i]=T(0.0);" if npath else "")]
    H += b + ["}", ""]

    b, u = emit_instructions(L_f, ["states","controls","parameters","time"],
        lambda a,e: f"_o[{e}]")
    H += ["template<class T> inline T gen_integrand_cost(",
          "        const T* states, const T* controls, const T* parameters, const T* time) {",
          _using(u) + "    T _o[1]; _o[0]=T(0.0);"]
    H += b + ["    return _o[0];", "}", ""]

    b, u = emit_instructions(phi_f, ["initial_states","final_states","parameters","t0","tf"],
        lambda a,e: f"_o[{e}]")
    H += ["template<class T> inline T gen_endpoint_cost(const T* initial_states,",
          "        const T* final_states, const T* parameters, const T* t0, const T* tf) {",
          _using(u) + "    T _o[1]; _o[0]=T(0.0);"]
    H += b + ["    return _o[0];", "}", ""]

    b, u = emit_instructions(ev_f, ["initial_states","final_states","parameters","t0","tf"],
        lambda a,e: f"e[{e}]")
    H += ["template<class T> inline void gen_events(T* e, const T* initial_states,",
          "        const T* final_states, const T* parameters, const T* t0, const T* tf) {",
          _using(u) + f"    for(int _i=0;_i<{nev};_i++) e[_i]=T(0.0);"]
    H += b + ["}", ""]

    if obs_f is not None:
        b, u = emit_instructions(obs_f, ["states","controls","parameters","time"],
            lambda a,e: f"observations[{e}]")
        H += ["template<class T> inline void gen_observation(T* observations,",
              "        const T* states, const T* controls, const T* parameters, const T* time) {",
              _using(u) + f"    for(int _i=0;_i<{dims['nobs']};_i++) observations[_i]=T(0.0);"]
        H += b + ["}", ""]

    H += ["#endif"]
    return "\n".join(H)


def emit_multiphase_header(name, phases):
    """Multi-phase emission. `phases` is a list (phase 1..N) of dicts with keys
    nx, npath, nevents, dae_f, L_f, phi_f, ev_f. Emits per-phase gen_*_pK functions
    plus iphase dispatchers gen_dae / gen_integrand_cost / gen_endpoint_cost / gen_events."""
    G = f"GEN_{name.upper()}_MP_HPP"
    H = [f"// {name} multi-phase -- AUTO-GENERATED from CasADi traces. Do not edit.",
         f"#ifndef {G}", f"#define {G}", "#include <cmath>", ""]
    N = len(phases)
    for k, ph in enumerate(phases, start=1):
        nx, npath, nev = ph["nx"], ph["npath"], ph["nevents"]
        b, u = emit_instructions(ph["dae_f"], ["states", "controls", "parameters", "time"],
            lambda a, e: f"derivatives[{e}]" if a == 0 else f"path[{e}]")
        init = f"    for(int _i=0;_i<{nx};_i++) derivatives[_i]=T(0.0);"
        if npath:
            init += f"\n    for(int _i=0;_i<{npath};_i++) path[_i]=T(0.0);"
        H += [f"template<class T> inline void gen_dae_p{k}(T* derivatives, T* path,",
              "        const T* states, const T* controls, const T* parameters, const T* time) {",
              _using(u) + init] + b + ["}", ""]

        b, u = emit_instructions(ph["L_f"], ["states", "controls", "parameters", "time"],
            lambda a, e: f"_o[{e}]")
        H += [f"template<class T> inline T gen_integrand_cost_p{k}(",
              "        const T* states, const T* controls, const T* parameters, const T* time) {",
              _using(u) + "    T _o[1]; _o[0]=T(0.0);"] + b + ["    return _o[0];", "}", ""]

        b, u = emit_instructions(ph["phi_f"], ["initial_states", "final_states", "parameters", "t0", "tf"],
            lambda a, e: f"_o[{e}]")
        H += [f"template<class T> inline T gen_endpoint_cost_p{k}(const T* initial_states,",
              "        const T* final_states, const T* parameters, const T* t0, const T* tf) {",
              _using(u) + "    T _o[1]; _o[0]=T(0.0);"] + b + ["    return _o[0];", "}", ""]

        if nev > 0:
            b, u = emit_instructions(ph["ev_f"], ["initial_states", "final_states", "parameters", "t0", "tf"],
                lambda a, e: f"e[{e}]")
            H += [f"template<class T> inline void gen_events_p{k}(T* e, const T* initial_states,",
                  "        const T* final_states, const T* parameters, const T* t0, const T* tf) {",
                  _using(u) + f"    for(int _i=0;_i<{nev};_i++) e[_i]=T(0.0);"] + b + ["}", ""]

    def branch(k, first):
        return ("if" if first else "else if") + f"(iphase=={k})"

    H += ["template<class T> inline void gen_dae(int iphase, T* derivatives, T* path,",
          "        const T* states, const T* controls, const T* parameters, const T* time) {"]
    for k in range(1, N + 1):
        H.append(f"    {branch(k, k==1)} gen_dae_p{k}(derivatives,path,states,controls,parameters,time);")
    H += ["}", ""]

    H += ["template<class T> inline T gen_integrand_cost(int iphase, const T* states,",
          "        const T* controls, const T* parameters, const T* time) {"]
    for k in range(1, N + 1):
        H.append(f"    {branch(k, k==1)} return gen_integrand_cost_p{k}(states,controls,parameters,time);")
    H += ["    return T(0.0);", "}", ""]

    H += ["template<class T> inline T gen_endpoint_cost(int iphase, const T* initial_states,",
          "        const T* final_states, const T* parameters, const T* t0, const T* tf) {"]
    for k in range(1, N + 1):
        H.append(f"    {branch(k, k==1)} return gen_endpoint_cost_p{k}(initial_states,final_states,parameters,t0,tf);")
    H += ["    return T(0.0);", "}", ""]

    H += ["template<class T> inline void gen_events(int iphase, T* e, const T* initial_states,",
          "        const T* final_states, const T* parameters, const T* t0, const T* tf) {"]
    for k in range(1, N + 1):
        if phases[k - 1]["nevents"] > 0:
            H.append(f"    if(iphase=={k}) gen_events_p{k}(e,initial_states,final_states,parameters,t0,tf);")
    H += ["}", "", "#endif"]
    return "\n".join(H)
