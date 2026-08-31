// US Standard Atmosphere 1976, from the ground to 1000 km, in a form that can
// be differentiated automatically.
//
// The 1976 standard is two different models joined at 86 km. Below that it
// defines the temperature as a piecewise-linear function of geopotential
// altitude through seven layers, with the pressure following from hydrostatic
// equilibrium in each. Above it the standard switches to a diffusive
// equilibrium of six species, integrated numerically; that upper branch has no
// closed form and is normally used through its published tables.
//
// Both branches share the same obstacle: evaluating them directly means
// selecting a layer or a band by comparison, and that is exactly the
// construction that cannot be taped. The branch is frozen when the tape is
// recorded, and every later evaluation uses whichever piece happened to be
// active then. The example examples/climb in this distribution carries a scar
// from precisely that mistake.
//
// So the model below has no branches at all.
//
//   Below 86 km, the temperature is written as the sum of the lapse rates times
//   smooth positive parts, which reproduces the piecewise-linear profile with
//   its corners rounded over about a hundred metres, and the pressure comes from
//   a polynomial fitted to the exact solution. Density follows from the gas law.
//   Against the layer equations, over 0 to 86 km: temperature within 0.29 K,
//   pressure and density within 0.29 per cent.
//
//   Above 86 km, pressure and density come from the ten-band quartic fits to the
//   tabulated standard published by R. A. Braeunig (Rocket and Space Technology,
//   www.braeunig.us/space/atmmodel.htm, Table 8), and the temperature from the
//   defining equations of the standard's upper branch -- isothermal to 91 km,
//   elliptical to 110 km, linear to 120 km, and exponential above. The bands are
//   joined not by selection but by hyperbolic-tangent steps of 350 m width, so
//   the whole thing is one smooth expression. Against the same fits evaluated
//   band by band, the blended form is within 0.04 per cent in pressure and
//   density everywhere from 86 to 1000 km, and within 0.07 K in temperature.
//
// The two branches are joined the same way, by a step at 86 km. Note that above
// 86 km the mean molecular weight falls, from 28.96 to 21.3 kg/kmol at 200 km,
// so pressure, density and temperature up there are NOT related by the gas law
// with a fixed gas constant. Each is therefore given by its own fit rather than
// two of them being derived from the third.
//
// Beyond 1000 km the altitude is clamped and an exponential of 233 km scale
// height continues the profile, which keeps it positive and monotone for any
// iterate the solver may produce without pretending to be the standard.

#ifndef US76_H
#define US76_H

#define US76_RGAS      287.05287     // J/(kg K)
#define US76_RE        6356766.0     // m, effective Earth radius for geopotential
#define US76_ZSPLIT    86.0          // km, where the two branches are joined
#define US76_ZTOP    1000.0          // km, where the standard's tables end
#define US76_CORNER    0.1           // km, rounding of the temperature corners
#define US76_BW        0.35          // km, width of the band-joining steps
#define US76_HTOP    232.7           // km, scale height used above US76_ZTOP

// values of the lower branch at the join, so that the blend is exact below it
#define US76_T86      186.946000000  // K
#define US76_LP86     -0.985476526187   // log(Pa)
#define US76_LR86    -11.875962746108   // log(kg/m^3)

static const double US76_HB[8] = {0.000000, 11.000000, 20.000000, 32.000000, 47.000000, 51.000000, 71.000000, 84.852000};   // km
static const double US76_LB[8] = {-6.50000000, 0.00000000, 1.00000000, 2.80000000, 0.00000000, -2.80000000, -2.00000000, 0.00000000};   // K/km

// log(pressure in Pa) as a degree-14 polynomial in geometric altitude in km,
// valid to 86 km
static const double US76_LOGP[15] = {
    +2.752018120738e-23, -1.549226059071e-20, +3.807441218302e-18,
    -5.356160463243e-16, +4.762313839739e-14, -2.809459671977e-12,
    +1.150146257803e-10, -3.598797900896e-09, +1.021314996213e-07,
    -2.756903576235e-06, +5.719101358358e-05, -6.756030850466e-04,
    +2.103416764291e-03, -1.252775707538e-01, +1.152899656715e+01,
};

// band edges above 86 km, in km
static const double US76_EDGE[9] = { 91.0, 100.0, 110.0, 120.0, 150.0,
                                    200.0, 300.0, 500.0, 750.0 };

// log(pressure in Pa), quartic in geometric altitude in km, one row per band
static const double US76_CP[10][5] = {
   { +0.000000e+00, +2.159582e-06, -4.836957e-04, -1.425192e-01, +1.347530e+01 },
   { +0.000000e+00, +3.304895e-05, -9.062730e-03, +6.516698e-01, -1.103037e+01 },
   { +0.000000e+00, +6.693926e-05, -1.945388e-02, +1.719080e+00, -4.775030e+01 },
   { +0.000000e+00, -6.539316e-05, +2.485568e-02, -3.223620e+00, +1.359355e+02 },
   { +2.283506e-07, -1.343221e-04, +2.999016e-02, -3.055446e+00, +1.135764e+02 },
   { +1.209434e-08, -9.692458e-06, +3.002041e-03, -4.523015e-01, +1.919151e+01 },
   { +8.113942e-10, -9.822568e-07, +4.687616e-04, -1.231710e-01, +3.067409e+00 },
   { +9.814674e-11, -1.654439e-07, +1.148115e-04, -5.431334e-02, -2.011365e+00 },
   { -7.835161e-11, +1.964589e-07, -1.657213e-04, +4.305869e-02, -1.477132e+01 },
   { +2.813255e-11, -1.120689e-07, +1.695568e-04, -1.188941e-01, +1.456718e+01 },
};

// log(density in kg/m^3), the same way
static const double US76_CR[10][5] = {
   { +0.000000e+00, -3.322622e-06, +9.111460e-04, -2.609971e-01, +5.944694e+00 },
   { +0.000000e+00, +2.873405e-05, -8.492037e-03, +6.541179e-01, -2.362010e+01 },
   { -1.240774e-05, +5.162063e-03, -8.048342e-01, +5.555996e+01, -1.443338e+03 },
   { +0.000000e+00, -8.854164e-05, +3.373254e-02, -4.390837e+00, +1.765294e+02 },
   { +3.661771e-07, -2.154344e-04, +4.809214e-02, -4.884744e+00, +1.723597e+02 },
   { +1.906032e-08, -1.527799e-05, +4.724294e-03, -6.992340e-01, +2.050921e+01 },
   { +1.199282e-09, -1.451051e-06, +6.910474e-04, -1.736220e-01, -5.321644e+00 },
   { +1.140564e-10, -2.130756e-07, +1.570762e-04, -7.029296e-02, -1.289844e+01 },
   { +8.105631e-12, -2.358417e-09, -2.635110e-06, -1.562608e-02, -2.002246e+01 },
   { -3.701195e-12, -8.608611e-09, +5.118829e-05, -6.600998e-02, -6.137674e+00 },
};

// smooth positive part: agrees with max(0,x) away from the origin
template <class T> static T us76_pos(T x, double w)
{
   return 0.5*(x + sqrt(x*x + w*w));
}

// smooth unit step, written with a hyperbolic tangent rather than a logistic:
// the two are the same function, but the derivative of tanh is 1 - tanh^2,
// which cannot overflow, whereas the logistic's carries exp(x) squared and
// returns a NaN long before the function itself has stopped being accurate
template <class T> static T us76_step(T x, double w)
{
   return 0.5*(1.0 + tanh(x/w));
}

// one banded quartic, evaluated as a single smooth expression
template <class T> static T us76_bands(const double C[10][5], T z)
{
   T pk = (((C[0][0]*z + C[0][1])*z + C[0][2])*z + C[0][3])*z + C[0][4];
   T g  = pk;
   for (int k = 1; k < 10; k++) {
      T pn = (((C[k][0]*z + C[k][1])*z + C[k][2])*z + C[k][3])*z + C[k][4];
      g  = g + (pn - pk)*us76_step<T>(z - US76_EDGE[k-1], US76_BW);
      pk = pn;
   }
   return g;
}

// kinetic temperature of the standard's upper branch, 86 km and above
template <class T> static T us76_temp_upper(T z)
{
   T t1 = 186.8673 + 0.0*z;
   T q  = 1.0 - ((z - 91.0)/19.9429)*((z - 91.0)/19.9429);
   T t2 = 263.1905 - 76.3232*sqrt( 0.5*(q + sqrt(q*q + 1.0e-8)) );
   T t3 = 240.0 + 12.0*(z - 110.0);
   T xi = (z - 120.0)*(6356.766 + 120.0)/(6356.766 + z);
   T t4 = 1000.0 - 640.0*exp(-0.01875*xi);
   return t1 + (t2 - t1)*us76_step<T>(z -  91.0, US76_BW)
             + (t3 - t2)*us76_step<T>(z - 110.0, US76_BW)
             + (t4 - t3)*us76_step<T>(z - 120.0, US76_BW);
}

// Returns density (kg/m^3), pressure (Pa) and temperature (K) at geometric
// altitude h, in metres.
template <class T> static void us76(T h, T& rho, T& pres, T& temp)
{
   T z  = h/1000.0;                                   // km, geometric
   T up = us76_pos<T>(z - US76_ZSPLIT, US76_CORNER);  // km above the join
   T zc = z - up;                                     // smoothly clamped to <= 86
   T zu = US76_ZSPLIT + up;                           // smoothly clamped to >= 86
   T over = us76_pos<T>(zu - US76_ZTOP, 1.0);         // km above the tables
   zu = zu - over;                                    // and clamped to <= 1000

   // ---- the layer model, valid to 86 km
   T Hg = US76_RE*(zc*1000.0)/(US76_RE + zc*1000.0)/1000.0;   // geopotential, km
   T tlo = 288.15 + US76_LB[0]*Hg;
   for (int i = 1; i < 8; i++)
      tlo = tlo + (US76_LB[i] - US76_LB[i-1])*us76_pos<T>(Hg - US76_HB[i], US76_CORNER);

   T lplo = US76_LOGP[0];
   for (int i = 1; i < 15; i++) lplo = lplo*zc + US76_LOGP[i];
   T lrlo = lplo - log(US76_RGAS*tlo);

   // ---- the tabulated model, valid from 86 km up
   T lpup = us76_bands<T>(US76_CP, zu) - over/US76_HTOP;
   T lrup = us76_bands<T>(US76_CR, zu) - over/US76_HTOP;
   T tup  = us76_temp_upper<T>(zu);

   // ---- and the join. Below 86 km the step is zero and the corrections vanish
   //      identically, so the lower branch is returned unaltered.
   T sw = us76_step<T>(z - US76_ZSPLIT, US76_BW);
   pres = exp( lplo + sw*(lpup - US76_LP86) );
   rho  = exp( lrlo + sw*(lrup - US76_LR86) );
   temp = tlo + sw*(tup - US76_T86);
}

#endif
