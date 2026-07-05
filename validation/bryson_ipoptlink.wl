(* ::Package:: *)
(* Wolfram IPOPTLink solve of Bryson-Denham, for cross-validation vs PSOPT and
   native CasADi. Run: wolframscript -file bryson_ipoptlink.wl [N]              *)
Needs["IPOPTLink`"];

nN = 200; dt = 1/nN;

x1[i_] := Symbol["x1n" <> ToString[i]];
x2[i_] := Symbol["x2n" <> ToString[i]];
u[i_]  := Symbol["un"  <> ToString[i]];

vars = Join[Table[x1[i], {i, 0, nN}], Table[x2[i], {i, 0, nN}], Table[u[i], {i, 0, nN}]];

big = 10.^19;
vb[i_] := Which[i == 0, {0., 0.}, i == nN, {0., 0.}, True, {-big, 1./9.}]; (* x1 path/bcs *)
wb[i_] := Which[i == 0, {1., 1.}, i == nN, {-1., -1.}, True, {-big, big}]; (* x2 bcs *)
varBounds = Join[Table[vb[i], {i, 0, nN}], Table[wb[i], {i, 0, nN}], Table[{-big, big}, {i, 0, nN}]];

x0 = ConstantArray[0.05, Length[vars]];

cons = Flatten[Table[{
    x1[i + 1] - x1[i] - dt/2 (x2[i] + x2[i + 1]),
    x2[i + 1] - x2[i] - dt/2 (u[i]  + u[i + 1])}, {i, 0, nN - 1}]];
conBounds = Table[{0., 0.}, {Length[cons]}];

f = Sum[dt/2 (0.5 u[i]^2 + 0.5 u[i + 1]^2), {i, 0, nN - 1}];

res = IPOPTMinimize[f, vars, x0, varBounds, cons, conBounds];
Print["return code : ", IPOPTReturnCode[res], "  (", IPOPTStringStatus[res], ")"];
Print["IPOPTLink   Bryson-Denham  N=", nN, "  objective = ", NumberForm[IPOPTMinValue[res], 8]];
