"""launch (4-phase Delta-III ascent) via the PSOPT Python API (B-2 acceptance test).
Target objective (matches hand-written C++ / CI reference): -7.529661e+03."""
import sys, math, numpy as np, casadi as ca
sys.path.insert(0, "/tmp/psopt_py")
sys.path.insert(0, "/tmp/psopt_py/examples")
import psopt
import rv2oe_casadi as R

# ---------------- constants (folded from launch.cxx user_data) ----------------
omega=7.29211585e-5; mu=3.986012e14; cd=0.5; sa=4*math.pi; rho0=1.225
H=7200.0; Re=6378145.0; g0=9.80665
thrust_srb=628500.0; thrust_first=1083100.0; thrust_second=110094.0
bt_srb=75.2; bt_first=261.0; bt_second=700.0
m_prop_srb=17010.0; m_prop_first=95550.0; m_prop_second=16820.0
mdot_srb=m_prop_srb/bt_srb; mdot_first=m_prop_first/bt_first; mdot_second=m_prop_second/bt_second
ISP_srb=thrust_srb/(g0*mdot_srb); ISP_first=thrust_first/(g0*mdot_first); ISP_second=thrust_second/(g0*mdot_second)
m_tot_srb=19290.0; m_tot_first=104380.0; m_tot_second=19300.0; m_payload=4164.0
m_dry_srb=m_tot_srb-m_prop_srb; m_dry_first=m_tot_first-m_prop_first
t1,t2,t3,t4=75.2,150.4,261.0,961.0
lat0=28.5*math.pi/180.0
x0=Re*math.cos(lat0); z0=Re*math.sin(lat0); y0=0.0
v0=[-omega*y0, omega*x0, 0.0]                       # omega_matrix * r0
m10=m_payload+m_tot_second+m_tot_first+9*m_tot_srb
m1f=m10-(6*mdot_srb+mdot_first)*t1
m20=m1f-6*m_dry_srb;  m2f=m20-(3*mdot_srb+mdot_first)*(t2-t1)
m30=m2f-3*m_dry_srb;  m3f=m30-mdot_first*(t3-t2)
m40=m3f-m_dry_first;  m4f=m_payload
af,ef,incf,Omf,omf=24361140.0,0.7308,28.5*math.pi/180.0,269.8*math.pi/180.0,130.5*math.pi/180.0
rmin,rmax,vmin,vmax=-2*Re,2*Re,-10000.0,10000.0

def oe2rv(oe, mu):
    a,e,i,Om,om,nu=oe
    p=a*(1-e*e); r=p/(1+e*math.cos(nu))
    rv=np.array([r*math.cos(nu), r*math.sin(nu), 0.0])
    vv=np.array([-math.sin(nu), e+math.cos(nu), 0.0])*math.sqrt(mu/p)
    cO,sO=math.cos(Om),math.sin(Om); co,so=math.cos(om),math.sin(om); ci,si=math.cos(i),math.sin(i)
    Rm=np.array([[cO*co-sO*so*ci, -cO*so-sO*co*ci, sO*si],
                 [sO*co+cO*so*ci, -sO*so+cO*co*ci, -cO*si],
                 [so*si, co*si, ci]])
    return Rm@rv, Rm@vv
rout,vout=oe2rv([af,ef,incf,Omf,omf,0.0], mu)

def phase_T_mdot(K):
    if K==1: return 6*thrust_srb+thrust_first, -(6*thrust_srb/(g0*ISP_srb)+thrust_first/(g0*ISP_first))
    if K==2: return 3*thrust_srb+thrust_first, -(3*thrust_srb/(g0*ISP_srb)+thrust_first/(g0*ISP_first))
    if K==3: return thrust_first, -thrust_first/(g0*ISP_first)
    return thrust_second, -thrust_second/(g0*ISP_second)

def make_dynamics(T_tot, mdot):
    def f(x,u,p,t):
        r=x[0:3]; v=x[3:6]; m=x[6]
        rad=ca.sqrt(r[0]**2+r[1]**2+r[2]**2)
        vrel=ca.vertcat(v[0]+omega*r[1], v[1]-omega*r[0], v[2])
        speedrel=ca.sqrt(vrel[0]**2+vrel[1]**2+vrel[2]**2)
        rho=rho0*ca.exp(-(rad-Re)/H)
        bcspeed=(rho/(2*m))*sa*cd*speedrel
        Drag=ca.vertcat(*[-(vrel[j]*bcspeed) for j in range(3)])
        muoverradcubed=mu/(rad**3)
        grav=ca.vertcat(*[-muoverradcubed*r[j] for j in range(3)])
        Toverm=T_tot/m
        vdot=ca.vertcat(*[Toverm*u[j]+Drag[j]+grav[j] for j in range(3)])
        return ca.vertcat(v[0],v[1],v[2], vdot[0],vdot[1],vdot[2], mdot)
    return f

path_fn=lambda x,u,p,t: ca.vertcat(u[0]**2+u[1]**2+u[2]**2)

prob=psopt.Problem(name="launch")
nodes={1:[15,18],2:[15,18],3:[15,18],4:[20,25]}
m0={1:m10,2:m20,3:m30,4:m40}; mf={1:m1f,2:m2f,3:m3f,4:m4f}
t0b={1:(0.0,0.0),2:(t1,t1),3:(t2,t2),4:(t3,t3)}
tfb={1:(t1,t1),2:(t2,t2),3:(t3,t3),4:(t3,t4)}
tspan={1:(0.0,t1),2:(t1,t2),3:(t2,t3),4:(t3,t4)}

for K in range(1,5):
    T,md=phase_T_mdot(K)
    nev=7 if K==1 else (5 if K==4 else 0)
    ph=prob.add_phase(nstates=7, ncontrols=3, nevents=nev, npath=1)
    ph.nodes=nodes[K]; ph.dynamics=make_dynamics(T,md); ph.path=path_fn
    if K==1:
        ph.events=lambda xi,xf,p,t0,tf: ca.vertcat(xi[0],xi[1],xi[2],xi[3],xi[4],xi[5],xi[6])
    if K==4:
        ph.endpoint=lambda xi,xf,p,t0,tf: -xf[6]
        ph.events=lambda xi,xf,p,t0,tf: R.rv2oe(xf[0:3], xf[3:6], mu)
    ph.bounds.lower.states=[rmin,rmin,rmin,vmin,vmin,vmin,mf[K]]
    ph.bounds.upper.states=[rmax,rmax,rmax,vmax,vmax,vmax,m0[K]]
    ph.bounds.lower.controls=[-1.0,-1.0,-1.0]; ph.bounds.upper.controls=[1.0,1.0,1.0]
    ph.bounds.lower.path=[1.0]; ph.bounds.upper.path=[1.0]
    if K==1:
        ph.bounds.lower.events=[x0,y0,z0,v0[0],v0[1],v0[2],m10]
        ph.bounds.upper.events=[x0,y0,z0,v0[0],v0[1],v0[2],m10]
    if K==4:
        ph.bounds.lower.events=[af,ef,incf,Omf,omf]; ph.bounds.upper.events=[af,ef,incf,Omf,omf]
    ph.bounds.t0=t0b[K]; ph.bounds.tf=tfb[K]
    gs=np.zeros((7,5))
    if K==4:
        gs[0:3,:]=np.array(rout).reshape(3,1); gs[3:6,:]=np.array(vout).reshape(3,1)
    else:
        gs[0,:]=x0; gs[1,:]=y0; gs[2,:]=z0; gs[3,:]=v0[0]; gs[4,:]=v0[1]; gs[5,:]=v0[2]
    gs[6,:]=np.linspace(m0[K],mf[K],5)
    ph.guess.states=gs
    gc=np.zeros((3,5)); gc[0,:]=1.0; ph.guess.controls=gc
    ph.guess.time=np.linspace(tspan[K][0],tspan[K][1],5).reshape(1,5)

prob.link_phases(1,2,jumps={6:6*m_dry_srb})
prob.link_phases(2,3,jumps={6:3*m_dry_srb})
prob.link_phases(3,4,jumps={6:m_dry_first})

alg=psopt.Algorithm(collocation_method="Chebyshev", nlp_iter_max=1000, nlp_tolerance=1e-6)
sol=prob.solve(alg)
print("\n================ RESULT ================")
print("objective    :", repr(sol.objective))
print("target (C++) : -7.529661e+03")
print("phases       :", len(sol.states), "| phase-4 states", sol.states[3].shape)
