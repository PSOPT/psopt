"""rv2oe transcribed into CasADi (USE_SMOOTH_HEAVISIDE path; oe[0..4]=a,e,i,Om,om)."""
import casadi as ca, math
EPS = 2.220446049250313e-16   # std::numeric_limits<double>::epsilon()
A_EPS = 0.1
def _cross(x,y):
    return ca.vertcat(x[1]*y[2]-x[2]*y[1], x[2]*y[0]-x[0]*y[2], x[0]*y[1]-x[1]*y[0])
def _dot(x,y): return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]
def _sh(x):    return 0.5*(1.0+ca.tanh(x/A_EPS))   # smooth_heaviside
def rv2oe(rv, vv, mu):
    K  = ca.vertcat(0.0,0.0,1.0)
    hv = _cross(rv,vv)
    nv = _cross(K,hv)
    n  = ca.sqrt(_dot(nv,nv))
    h2 = _dot(hv,hv)
    v2 = _dot(vv,vv)
    r  = ca.sqrt(_dot(rv,rv))
    rvvv = _dot(rv,vv)
    ev = ca.vertcat(*[(1.0/mu)*((v2-mu/r)*rv[j] - rvvv*vv[j]) for j in range(3)])
    p  = h2/mu
    e  = ca.sqrt(_dot(ev,ev))
    a  = p/(1.0-e*e)
    i  = ca.acos(hv[2]/ca.sqrt(h2))
    ac_Om = ca.acos(nv[0]/n)
    Om = _sh(nv[1]+EPS)*ac_Om + _sh(-(nv[1]+EPS))*(2*math.pi-ac_Om)
    ac_om = ca.acos(_dot(nv,ev)/n/e)
    om = _sh(ev[2])*ac_om + _sh(-ev[2])*(2*math.pi-ac_om)
    return ca.vertcat(a,e,i,Om,om)
