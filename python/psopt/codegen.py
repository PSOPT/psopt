"""codegen.py -- Layer B: CasADi-traced user maths -> adouble-templated C++ ->
per-problem shared library exporting the PSOPT user functions as extern "C" symbols.
The .so is cached by a hash of its source so re-solving the same problem skips the
compile."""
import hashlib, os, subprocess
from . import _emitter
from ._toolchain import CXX, CXXFLAGS, INCLUDES, LIBS, CACHE_DIR

# Thin shims at the fixed PSOPT signatures; bodies call the generated gen_* templates.
_SHIMS = r'''
extern "C" void psopt_dae(adouble* derivatives, adouble* path, adouble* states,
        adouble* controls, adouble* parameters, adouble& time, adouble* xad,
        int iphase, Workspace* workspace) {
    gen_dae<adouble>(derivatives, path, states, controls, parameters, &time);
}
extern "C" adouble psopt_integrand_cost(adouble* states, adouble* controls,
        adouble* parameters, adouble& time, adouble* xad, int iphase, Workspace* workspace) {
    return gen_integrand_cost<adouble>(states, controls, parameters, &time);
}
extern "C" adouble psopt_endpoint_cost(adouble* initial_states, adouble* final_states,
        adouble* parameters, adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace) {
    return gen_endpoint_cost<adouble>(initial_states, final_states, parameters, &t0, &tf);
}
extern "C" void psopt_events(adouble* e, adouble* initial_states, adouble* final_states,
        adouble* parameters, adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace) {
    gen_events<adouble>(e, initial_states, final_states, parameters, &t0, &tf);
}
extern "C" void psopt_linkages(adouble* linkages, adouble* xad, Workspace* workspace) { }
'''

_OBS_SHIM = r'''
extern "C" void psopt_observation(adouble* observed_variable, adouble* states,
        adouble* controls, adouble* parameters, adouble& time, int k, adouble* xad,
        int iphase, Workspace* workspace) {
    gen_observation<adouble>(observed_variable, states, controls, parameters, &time);
}
'''

def compile_single_phase(name, dims, dae_f, L_f, phi_f, ev_f, obs_f=None):
    """Emit + JIT-compile the math library for a single-phase problem; return .so path.
    If obs_f is given (parameter estimation), also emit gen_observation + its shim."""
    header = _emitter.emit_header(name, dims, dae_f, L_f, phi_f, ev_f, obs_f=obs_f)
    source = '#include "psopt.h"\nusing namespace PSOPT;\n' + header + "\n" + _SHIMS
    if obs_f is not None:
        source += _OBS_SHIM
    key = hashlib.sha1(source.encode()).hexdigest()[:16]
    os.makedirs(CACHE_DIR, exist_ok=True)
    so = os.path.join(CACHE_DIR, f"{name}_{key}.so")
    if os.path.exists(so):
        return so
    cpp = os.path.join(CACHE_DIR, f"{name}_{key}.cpp")
    with open(cpp, "w") as f:
        f.write(source)
    cmd = [CXX, *CXXFLAGS, "-fPIC", "-shared", *INCLUDES, cpp, "-o", so, *LIBS]
    subprocess.run(cmd, check=True)
    return so


_MULTI_SHIMS = r'''
extern "C" void psopt_dae(adouble* derivatives, adouble* path, adouble* states,
        adouble* controls, adouble* parameters, adouble& time, adouble* xad,
        int iphase, Workspace* workspace) {
    gen_dae<adouble>(iphase, derivatives, path, states, controls, parameters, &time);
}
extern "C" adouble psopt_integrand_cost(adouble* states, adouble* controls,
        adouble* parameters, adouble& time, adouble* xad, int iphase, Workspace* workspace) {
    return gen_integrand_cost<adouble>(iphase, states, controls, parameters, &time);
}
extern "C" adouble psopt_endpoint_cost(adouble* initial_states, adouble* final_states,
        adouble* parameters, adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace) {
    return gen_endpoint_cost<adouble>(iphase, initial_states, final_states, parameters, &t0, &tf);
}
extern "C" void psopt_events(adouble* e, adouble* initial_states, adouble* final_states,
        adouble* parameters, adouble& t0, adouble& tf, adouble* xad, int iphase, Workspace* workspace) {
    gen_events<adouble>(iphase, e, initial_states, final_states, parameters, &t0, &tf);
}
extern "C" void psopt_linkages(adouble* linkages, adouble* xad, Workspace* workspace) {
/*LINKAGES*/
}
'''

def compile_multiphase(name, phases, links, nstates_by_phase):
    """phases: list of dicts (nx, npath, nevents, dae_f, L_f, phi_f, ev_f).
    links: list of dicts {a, b, jumps:{state_index: delta}} -> auto_link + jump adjust.
    Returns the cached .so path."""
    header = _emitter.emit_multiphase_header(name, phases)
    lines = ["    int index=0;"]
    for lk in links:
        a, b = lk["a"], lk["b"]
        lines.append(f"    auto_link(linkages,&index,xad,{a},{b},workspace);")
        ns = nstates_by_phase[b - 1]
        for s, delta in lk.get("jumps", {}).items():
            lines.append(f"    linkages[index - {ns + 1} + {s}] -= {float(delta):.17g};")
    shims = _MULTI_SHIMS.replace("/*LINKAGES*/", "\n".join(lines))
    source = '#include "psopt.h"\nusing namespace PSOPT;\n' + header + "\n" + shims
    key = hashlib.sha1(source.encode()).hexdigest()[:16]
    os.makedirs(CACHE_DIR, exist_ok=True)
    so = os.path.join(CACHE_DIR, f"{name}_mp_{key}.so")
    if os.path.exists(so):
        return so
    cpp = os.path.join(CACHE_DIR, f"{name}_mp_{key}.cpp")
    with open(cpp, "w") as f:
        f.write(source)
    subprocess.run([CXX, *CXXFLAGS, "-fPIC", "-shared", *INCLUDES, cpp, "-o", so, *LIBS], check=True)
    return so
