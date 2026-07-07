/* mioc_scip.c -- Mixed-Integer Optimal Control via SCIP.
 *
 * Demonstrates solving an optimal-control problem with INTEGER controls, which
 * PSOPT's continuous NLP backends (IPOPT/CasADi) cannot do. A double integrator
 * is transcribed with explicit Euler over N stages into a MILP and solved to
 * global optimality by SCIP:
 *
 *   states  pos_k, vel_k        (k = 0..N, continuous)
 *   control u_k in {-1,0,1}     (k = 0..N-1, INTEGER thrust)
 *   aux     a_k >= |u_k|        (for the actuation objective)
 *
 *   min  sum_k a_k                              (minimum total actuation)
 *   s.t. pos_{k+1} = pos_k + dt*vel_k           (dynamics)
 *        vel_{k+1} = vel_k + dt*u_k
 *        a_k >= u_k,  a_k >= -u_k
 *        pos_0=0, vel_0=0, pos_N=1, vel_N=0     (boundary conditions)
 *
 * SCIP returns the globally optimal integer (bang-off-bang) control schedule.
 *
 * Build: see build.sh (links the SCIP built into /opt/scip).
 */
#include <stdio.h>
#include <stdlib.h>
#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#define N   30
#define DT  0.10   /* horizon T = N*DT = 3.0; max rest-to-rest reach = T^2/4 = 2.25 > 1 */

static SCIP_RETCODE run(void)
{
   SCIP* scip;
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "mioc_double_integrator") );
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   SCIP_VAR *pos[N+1], *vel[N+1], *u[N], *a[N];
   char name[64];
   int k;

   /* state variables (continuous) with boundary conditions fixed via bounds */
   for( k = 0; k <= N; ++k )
   {
      double plo=-10, phi=10, vlo=-10, vhi=10;
      if( k==0 ) { plo=phi=0.0; vlo=vhi=0.0; }      /* pos_0=0, vel_0=0 */
      if( k==N ) { plo=phi=1.0; vlo=vhi=0.0; }      /* pos_N=1, vel_N=0 */
      sprintf(name,"pos_%d",k);
      SCIP_CALL( SCIPcreateVarBasic(scip,&pos[k],name,plo,phi,0.0,SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip,pos[k]) );
      sprintf(name,"vel_%d",k);
      SCIP_CALL( SCIPcreateVarBasic(scip,&vel[k],name,vlo,vhi,0.0,SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip,vel[k]) );
   }
   /* integer controls u_k in {-1,0,1}; actuation aux a_k>=0 with objective coeff 1 */
   for( k = 0; k < N; ++k )
   {
      sprintf(name,"u_%d",k);
      SCIP_CALL( SCIPcreateVarBasic(scip,&u[k],name,-1.0,1.0,0.0,SCIP_VARTYPE_INTEGER) );
      SCIP_CALL( SCIPaddVar(scip,u[k]) );
      sprintf(name,"a_%d",k);
      SCIP_CALL( SCIPcreateVarBasic(scip,&a[k],name,0.0,1.0,1.0,SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip,a[k]) );
   }

   /* dynamics: pos_{k+1} - pos_k - dt*vel_k = 0 ; vel_{k+1} - vel_k - dt*u_k = 0 */
   for( k = 0; k < N; ++k )
   {
      SCIP_CONS* c;
      SCIP_VAR*  vp[3]; SCIP_Real cp[3];
      vp[0]=pos[k+1]; cp[0]=1.0;  vp[1]=pos[k]; cp[1]=-1.0;  vp[2]=vel[k]; cp[2]=-DT;
      sprintf(name,"posdyn_%d",k);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip,&c,name,3,vp,cp,0.0,0.0) );
      SCIP_CALL( SCIPaddCons(scip,c) ); SCIP_CALL( SCIPreleaseCons(scip,&c) );

      vp[0]=vel[k+1]; cp[0]=1.0;  vp[1]=vel[k]; cp[1]=-1.0;  vp[2]=u[k]; cp[2]=-DT;
      sprintf(name,"veldyn_%d",k);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip,&c,name,3,vp,cp,0.0,0.0) );
      SCIP_CALL( SCIPaddCons(scip,c) ); SCIP_CALL( SCIPreleaseCons(scip,&c) );

      /* a_k >= u_k  and  a_k >= -u_k  =>  a_k - u_k >= 0 ,  a_k + u_k >= 0 */
      vp[0]=a[k]; cp[0]=1.0; vp[1]=u[k]; cp[1]=-1.0;
      sprintf(name,"absp_%d",k);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip,&c,name,2,vp,cp,0.0,SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip,c) ); SCIP_CALL( SCIPreleaseCons(scip,&c) );
      vp[0]=a[k]; cp[0]=1.0; vp[1]=u[k]; cp[1]=1.0;
      sprintf(name,"absn_%d",k);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip,&c,name,2,vp,cp,0.0,SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip,c) ); SCIP_CALL( SCIPreleaseCons(scip,&c) );
   }

   SCIP_CALL( SCIPsolve(scip) );

   SCIP_SOL* sol = SCIPgetBestSol(scip);
   if( sol != NULL )
   {
      printf("\n=== Mixed-integer optimal control (SCIP) ===\n");
      printf("optimal total actuation = %.4f\n", SCIPgetSolOrigObj(scip,sol));
      printf(" k   t      pos      vel     u(int)\n");
      for( k = 0; k < N; ++k )
         printf("%2d  %.2f  %7.4f  %7.4f   %+.0f\n", k, k*DT,
                SCIPgetSolVal(scip,sol,pos[k]), SCIPgetSolVal(scip,sol,vel[k]),
                SCIPgetSolVal(scip,sol,u[k]));
      printf("%2d  %.2f  %7.4f  %7.4f\n", N, N*DT,
             SCIPgetSolVal(scip,sol,pos[N]), SCIPgetSolVal(scip,sol,vel[N]));
   }

   for( k = 0; k <= N; ++k ){ SCIP_CALL(SCIPreleaseVar(scip,&pos[k])); SCIP_CALL(SCIPreleaseVar(scip,&vel[k])); }
   for( k = 0; k <  N; ++k ){ SCIP_CALL(SCIPreleaseVar(scip,&u[k]));   SCIP_CALL(SCIPreleaseVar(scip,&a[k]));   }
   SCIP_CALL( SCIPfree(&scip) );
   return SCIP_OKAY;
}

int main(void)
{
   if( run() != SCIP_OKAY ){ fprintf(stderr,"SCIP error\n"); return 1; }
   return 0;
}
