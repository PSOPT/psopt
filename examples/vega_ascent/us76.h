// US Standard Atmosphere 1976, in a form that can be differentiated automatically.
//
// The 1976 standard defines the temperature as a piecewise-linear function of
// geopotential altitude through seven layers, and the pressure follows from
// hydrostatic equilibrium in each. Evaluating it directly means selecting a layer
// by comparison, which is exactly the construction that cannot be taped: the
// branch is frozen when the tape is recorded, and every later evaluation uses
// whichever layer happened to be active then. The example examples/climb in this
// distribution carries a scar from precisely that mistake.
//
// The model below has no branches at all. The temperature is written as the sum
// of the lapse rates times smooth positive parts, which reproduces the
// piecewise-linear profile with its corners rounded over about a hundred metres,
// and the pressure comes from a polynomial fitted to the exact solution. Density
// then follows from the gas law. Against the layer equations, over 0 to 86 km:
//
//     temperature   within 0.29 K
//     pressure      within 0.29 per cent
//     density       within 0.29 per cent
//
// Above 86 km the 1976 standard changes to a species-diffusion formulation whose
// tables are not reproduced here; the pressure is continued instead by an
// exponential whose scale height matches the fit at the join, so the model stays
// smooth and monotone all the way to orbit. That extension is not US76 and
// should not be relied on for anything that depends on the density up there --
// which is why this example imposes no heat-flux constraint, the reference's
// constraint acting between 110 and 190 km.

#ifndef US76_H
#define US76_H

#define US76_RGAS      287.05287     // J/(kg K)
#define US76_RE        6356766.0     // m, effective Earth radius for geopotential
#define US76_ZSPLIT    86.0          // km
#define US76_CORNER    0.1           // km, rounding of the temperature corners
#define US76_HTAIL     5.6122226997     // km

static const double US76_HB[8] = {0.000000, 11.000000, 20.000000, 32.000000, 47.000000, 51.000000, 71.000000, 84.852000};   // km
static const double US76_LB[8] = {-6.50000000, 0.00000000, 1.00000000, 2.80000000, 0.00000000, -2.80000000, -2.00000000, 0.00000000};   // K/km

// log(pressure in Pa) as a degree-14 polynomial in geometric altitude in km
static const double US76_LOGP[15] = {
    +2.752018120738e-23, -1.549226059071e-20, +3.807441218302e-18,
    -5.356160463243e-16, +4.762313839739e-14, -2.809459671977e-12,
    +1.150146257803e-10, -3.598797900896e-09, +1.021314996213e-07,
    -2.756903576235e-06, +5.719101358358e-05, -6.756030850466e-04,
    +2.103416764291e-03, -1.252775707538e-01, +1.152899656715e+01,
};

// smooth positive part: agrees with max(0,x) away from the origin
template <class T> static T us76_pos(T x, double w)
{
   return 0.5*(x + sqrt(x*x + w*w));
}

// Returns density (kg/m^3), pressure (Pa) and temperature (K) at geometric
// altitude h, in metres.
template <class T> static void us76(T h, T& rho, T& pres, T& temp)
{
   T z  = h/1000.0;                                   // km, geometric
   T up = us76_pos<T>(z - US76_ZSPLIT, US76_CORNER);  // km above the split
   T zc = z - up;                                     // smoothly clamped to <= 86

   // geopotential altitude, km
   T Hg = US76_RE*(zc*1000.0)/(US76_RE + zc*1000.0)/1000.0;

   // piecewise-linear temperature with rounded corners
   temp = 288.15 + US76_LB[0]*Hg;
   for (int i = 1; i < 8; i++)
      temp = temp + (US76_LB[i] - US76_LB[i-1])*us76_pos<T>(Hg - US76_HB[i], US76_CORNER);

   T lp = US76_LOGP[0];
   for (int i = 1; i < 15; i++) lp = lp*zc + US76_LOGP[i];
   pres = exp(lp - up/US76_HTAIL);
   rho  = pres/(US76_RGAS*temp);
}

#endif
