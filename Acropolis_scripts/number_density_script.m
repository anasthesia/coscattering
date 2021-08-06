#!/usr/bin/env wolframscript


(*Reading the parameter file*)
Get["~/prj/acropolis/param.m"];
{m1, m2, mZp, tant, eps, ax} = ToExpression[StringReplace[#, {"e+" :> "*^", "e-" :> "*^-"}]] & /@ {m1, m2, mZp, tant, eps, ax};
(*I am rounding the numbers to get exact number instead of floats. This is necesarry to obtain enough precision for the Boltzmann code*)
ModelParam = {aX -> Round[ax*10^15]*10^-15, epsilon -> Round[eps*10^15]*10^-15, ttheta -> Round[tant*10^15]*10^-15, MZp -> Round[mZp*10^15]*10^-15, MD1 -> Round[m1*10^15]*10^-15, MD2 -> Round[m2*10^15]*10^-15, gD1 -> 2, gD2 -> 2};
Print["ModelParam=",ModelParam];

SetDirectory["~/prj/coscattering/Sam/MathematicaNotebooks/"];
Get["Coscattering.wl"];

res=DMabundance[ModelParam,xmin,xmax];
Print["res=",res];
(*Exporting the mediator number density to a file*)
xval = 10^Range[Log10[xmin], Log10[xmax], (Log10[xmax] - Log10[xmin])/49]//N;
nMed = Map[Ymed, xval]*Map[Entr, MD1/xval //. ModelParam]//N;
Export["~/prj/acropolis/number_density.dat", Transpose[{xval, nMed}]];