#!/usr/bin/env wolframscript


(*Reading the parameter file*)
Get["param.m"];

m1=ToExpression[StringReplace[m1,"e"->" 10^"]];
m2=ToExpression[StringReplace[m2,"e"->" 10^"]];
mZp=ToExpression[StringReplace[mZp,"e"->" 10^"]];
ax=ToExpression[StringReplace[ax,"e"->" 10^"]];
eps=ToExpression[StringReplace[eps,"e"->" 10^"]];
tant=ToExpression[StringReplace[tant,"e"->" 10^"]];

(*I am rounding the numbers to get exact number instead of floats. This is necesarry to obtain enough precision for the Boltzmann code*)
ModelParam = {aX -> Round[ax*10^5]*10^-5, epsilon -> Round[eps*10^5]*10^-5, ttheta -> Round[tant*10^5]*10^-5, MZp -> Round[mZp*10^5]*10^-5, MD1 -> Round[m1*10^5]*10^-5, MD2 -> Round[m2*10^5]*10^-5, gD1 -> 2, gD2 -> 2};

SetDirectory["/home/sam/OneDrive/Coscattering/"];
Get["Coscattering.wl"];

res=DMabundance[ModelParam,xmin,xmax];
Print[res];
(*Exporting the mediator number density to a file*)
xval = N[10^Range[Log10[xmin], Log10[xmax], (Log10[xmax] - Log10[xmin])/49],10];
nMed = Map[Ymed, xval]*Map[Entr, MD1/xval //. ModelParam]//N;
Export["/home/sam/Software/acropolis/number_density.dat", Transpose[{xval, nMed}]];