(* ::Package:: *)

(* ::Text:: *)
(*Initializing replacements to calculate the cross sections + SM input variables in units of MeV*)


SMMasses={MB->4200,MC->1270,MT->17300,MU->22*10^-1,MD->47*10^-1,MS->96,ME->511*10^-3,MMU->106,MTA->1777,MW->80400,MZ->91200,MH->125000,Ma->0,Mg->0,WZ->2495};
Relations={cw->Sqrt[1-swsq],sw->Sqrt[swsq],swsq->225*10^-3,aEW->10/1279,EE->2*Sqrt[aEW*Pi],gX->2*Sqrt[Pi aX]
,ta->-(-1+DZ+eta^2*sw^2+SignAux*Sqrt[4*eta^2*sw^2+(-1+DZ+eta^2*sw^2)^2])/(2*eta*sw),ca->1/Sqrt[1+ta^2],sa->ta/Sqrt[1+ta^2],eta->(epsilon/cw)/Sqrt[1-(epsilon/cw)^2]
,chi->epsilon/cw,DZ->(MZ^-2*MZp^-2*(MZ^4+MZp^4+MZp^2*(-2*eta^2*MZ^2*sw^2-SignAux*Sqrt[MZ^4+MZp^4-2*MZ^2*MZp^2*(1+2*eta^2*sw^2)])-MZ^2*SignAux*Sqrt[MZ^4+MZp^4-2*MZ^2*MZp^2*(1+2*eta^2*sw^2)]))/2,
SignAux->(MZ-MZp)/Sqrt[(MZ-MZp)^2],ctheta->1/Sqrt[1+ttheta^2],stheta->ttheta/Sqrt[1+ttheta^2],WZp->(EE^2*Sqrt[(-4*ME^2 + MZp^2)]*(cw^4*sa^2 - 2*cw^2*sa*sw*(ca*eta + sa*sw) 
+ 5*sw^2*(ca*eta + sa*sw)^2))/(96*cw^2*Pi*sw^2) +ca^2 ctheta^4 eta^2 gX^2 Sqrt[-4 MD2^2+MZp^2] (2 MD2^2+MZp^2)/(12 chi^2 MZp^2 \[Pi])HeavisideTheta[MZp-2*MD2]
+ca^2 stheta^4 eta^2 gX^2 Sqrt[-4 MD1^2+MZp^2] (2 MD1^2+MZp^2)/(12 chi^2 MZp^2 \[Pi])HeavisideTheta[MZp-2*MD1]
-(ca^2*ctheta^2*eta^2*gX^2*Sqrt[(-(MD1-MD2)^2+MZp^2)*(-(MD1+MD2)^2+MZp^2)]*(MD1^4+MD2^4-6*MD1*MD2*MZp^2+MD2^2*MZp^2-2*MZp^4+MD1^2*(-2*MD2^2+MZp^2))*stheta^2)/(12*chi^2*MZp^5*Pi)HeavisideTheta[MZp-(MD1+MD2)]};


(* ::Text:: *)
(*Initializing cosmological input (Hubble rate, DM eq density, ...)*)


<<NumericalCalculus`

{Temp,hEff,gEff}=Transpose[Import["std_thg.tab","Table"]];
gEffSTemp=SetPrecision[Transpose[{Temp*10^3,hEff}],21];
gEffTemp=SetPrecision[Transpose[{Temp*10^3,gEff}],21];

gEffSFunc[T_]:=If[T>gEffSTemp[[Length[gEffSTemp],1]],gEffSTemp[[Length[gEffSTemp],2]],Interpolation[gEffSTemp,T,InterpolationOrder->SetPrecision[1,Infinity]]];
gEffFunc[T_]:=If[T>gEffTemp[[Length[gEffTemp],1]],gEffTemp[[Length[gEffTemp],2]],Interpolation[gEffSTemp,T,InterpolationOrder->SetPrecision[1,Infinity]]];

Mpl=2435*10^15;(*Gev*)
rhoR[T_]:=gEffFunc[T]*Pi^2/30*T^4/.gEffc->10675/100;
H[T_]:=Sqrt[rhoR[T]/3]/Mpl;
Hbar[T_]:=H[T]*(1+ND[Log[gEffSFunc[x]],x,T]/(3T))^-1;

Entr[T_]:=2Pi^2/45 gEffSFunc[T]*T^3/.gEffc->10675/100;

neqDM[x_]:=gDM*MD1^3/(2Pi^2x)BesselK[2,x]/.gDM->2//.ModelParam;
YeqDM[x_]:=neqDM[x]/Entr[MD1/x]//.ModelParam;

neqMed[x_]:=gMed*MD2^2*MD1/(2Pi^2x)BesselK[2,MD2*x/MD1]/.gMed->2//.ModelParam;
YeqMed[x_]:=neqMed[x]/Entr[MD1/x]//.ModelParam;

neq[x_,m_]:=m^2MD1/(2Pi^2x)BesselK[2,m*x/MD1]//.ModelParam;


(* ::Text:: *)
(*Calculates the thermally averaged cross section for a list of parameters (ModelParam) which contains only the free variables, a certain temperature MD1/T (x) and a specific process denoted by a flag (1= Chi1 SM -> Chi1 SM, 2 = Chi2 Chi2 -> SM SM, 3 = Chi1 SM -> Chi2 SM). The function uses the files saved in folder m_files, where the 0 temperature cross sections are stored only for the first generation of fermions, while assuming they are massless. *)


CalcGaScat[ModelParam_,x_,Process_]:=Module[{IntM,fileName,file,Nf,itt,IntMt,IntMi,integrand,Param},
IntM=0;
fileName="Int_symb"<>ToString[Process]<>"*";
file=FileNames[fileName,NotebookDirectory[]<>"m_files"];
Nf=Length[file];

Param = Join[SMMasses,Relations,ModelParam];

For[itt=1,itt<Nf+1,itt++,
Get[file[[itt]]];
smin=SetPrecision[smin//.Param,31];
integrand= (E1cm^2-Mass1^2)\[Sigma] Sqrt[s]*BesselK[1,Sqrt[s]*x/MD1];

(*Print[Table[{10^(y/10),N[\[Sigma]//.Join[Param]/.{s\[Rule]10^(y/10)*smin},31]},{y,0,20}]];
Print[ListLogLogPlot[Table[{10^(y/10),N[integrand//.Join[Param]/.{s\[Rule]10^(y/10)*smin},101]},{y,0,20}],Joined\[Rule]True]];*)

IntMi=NIntegrate[N[integrand*10^-50//.Join[Param]/.{s->y*smin},31],{y,1,Infinity},WorkingPrecision->20]*10^50*smin;
IntM=IntM+IntMi;
];
Return[MD1/(8 Pi^4 x)*IntM//.ModelParam];
 ]


(* ::Text:: *)
(*Calculates the thermally averaged  cross section times n1 n2 for a specific set of parameters and a specific process denoted by a flag (as above) for the temperatures stored in the list xval*)


AveragedXsec[ModelParam_,process_,xval_]:=Module[{points,svScat,x,y,i,smax,fileName,file},

fileName="Int_symb*";
file=FileNames[fileName,NotebookDirectory[]<>"m_files"];

svScat[x_]:=CalcGaScat[Join[ModelParam],x,process];
points={};
For[i=1,i<=Length[xval],i++,
y=Chop[svScat[xval[[i]]],10^-400];
If[y==0,Break[]];
points=Append[points,{xval[[i]],y}]];

Return[points];

]


(* ::Text:: *)
(*Solves the BE for a given set of parameters. It starts integrating from xmin to xmax*)


SolveCoupledBE[ModelParam_,xmin_,xmax_]:=Module[{l,dYDMinit,sol,sol1,solInit,YDMintermediate,dYDM,dYMed,YDM,YMed,consts,initCond,xInt},

l[x_]:=1/(Hbar[MD1/x]*Entr[MD1/x])/.ModelParam;

(*dYDMinit[x_]:=-x^4*l*CoscatRate[x]*(YDM[x]/YeqDM[x]-1);
xInt=10;
sol1=NDSolve[{YDM'[x]==dYDMinit[x],YDM[xmin]==YeqDM[xmin]},{YDM},{x,xmin,xInt},Method->"StiffnessSwitching"];
YDMintermediate=YDM[xInt]/.sol1[[1]];

xInt=xmin;
YDMintermediate=YeqDM[xInt];*)

dYDM[x_]:= -l[x]*CoscatRate[x]*(YDM[x]/YeqDM[x]-YMed[x]/YeqMed[x]);(*Should also include decays?*)
dYMed[x_]:=-l[x]*(MedAnnRate[x]*(YMed[x]^2/YeqMed[x]^2-1)-CoscatRate[x]*(YDM[x]/YeqDM[x]-YMed[x]/YeqMed[x]));

sol= NDSolve[{YDM'[x]==dYDM[x],YMed'[x]==dYMed[x],YDM[xmin]==YeqDM[xmin],YMed[xmin]==YeqMed[xmin]},{YMed,YDM},{x,xmin,xmax},Method->"StiffnessSwitching",WorkingPrecision->20,AccuracyGoal->25];

Ydm[x_]:=YDM[x]/.sol[[1]];
Ymed[x_]:=YMed[x]/.sol[[1]];

Return[(YDM[xmax]+YMed[xmax])/.sol[[1]]];
]


(* ::Text:: *)
(*Calculates the DM abundance for a specific set of model parameters*)


DMabundance[ModelParam_]:=Module[{xmin,xmax,xval,consts,DataCoscatRate,DataMedAnnRate,DatasvMedAnn,DataDMScatRate},
xmax=100;
xmin=1;
xval=Round[10^Range[Log10[xmin],Log10[xmax],0.04]*10^5]*10^-5;

DataCoscatRate=AveragedXsec[ModelParam,3,xval];
CoscatRate[x_]:=If[x<DataCoscatRate[[-1,1]],Exp[Interpolation[Log[DataCoscatRate],Log[x],InterpolationOrder->1]],0];

DataMedAnnRate=AveragedXsec[ModelParam,2,xval];
MedAnnRate[x_]:=If[x<DataMedAnnRate[[-1,1]],Exp[Interpolation[Log[DataMedAnnRate],Log[x],InterpolationOrder->1]],0];

DataDMScatRate=AveragedXsec[ModelParam,1,xval];
DMScatRate[x_]:=If[x<DataDMScatRate[[-1,1]],Exp[Interpolation[Log[DataDMScatRate],Log[x],InterpolationOrder->1]],0];

consts={Mpl->2.435*10^18,T0->2.348*10^-13,H0->2.131*10^-42,s0->(1+21/22)*2*Pi^2/45*2*T0^3,CritDens->3*H0^2*Mpl^2 };
Return[MD1*s0*SolveCoupledBE[ModelParam,xmin,xmax]/CritDens/.ModelParam//.consts];


]


DataForPlot[ModelParam_,x_]:=Module[{xmax2,xFO,CoscatRate,xmin,xmax,xval,DataCoscatRate},

(*xmax=400;
xmin=10^-3;
xval=Round[10^Range[Log10[xmin],Log10[xmax],0.05]*10^5]*10^-5;

DataCoscatRate=AveragedXsec[ModelParam,3,xval];
CoscatRate[x_]:=If[x<DataCoscatRate[[-1,1]],Exp[Interpolation[Log[DataCoscatRate],Log[x],InterpolationOrder->1]],0];
(*xmax2=pointlist[[3,-1,1]];*)

If[CoscatRate[xmin]/neqDM[xmin]<H[MD1/xmin]/.ModelParam,Return[{-1,-1,-1}]];
If[CoscatRate[xmax]/neqDM[xmax]>H[MD1/xmax]/.ModelParam,Return[{-2,-2,-2}]];

xFO=x/.FindRoot[CoscatRate[x]/(neqDM[x]*H[MD1/x])-1/.ModelParam,{x,xmin,xmax}];

If[!(xmin<xFO<xmax),Return[{-3,-3,-3}]];*)

Return[{CalcGaScat[ModelParam,x,2]/(neqDM[x]*H[MD1/x])(*,CalcGaScat[ModelParam,1,3]/(neqDM[1]*H[MD1])*),CalcGaScat[ModelParam,x,1]/(neqDM[x]*H[MD1/x])}/.ModelParam]

]


script:=Module[{},

data10={};
data20={};
data30={};
For[LogMZp=3,LogMZp<=3,LogMZp=LogMZp+1/2,
 For[LogTanT=-3,LogTanT<=0,LogTanT=LogTanT+1/4,
  For[LogEpsM=-2,LogEpsM<=1,LogEpsM=LogEpsM+1/4,
   ModelParam={aX->1,epsilon->10^-4,ttheta->10^LogTanT,MZp->10^LogMZp,MD1->10^2,MD2->MD1*(1+10^LogEpsM)};
   AppendTo[data10,Join[{LogMZp//N,LogEpsM//N,LogTanT//N},DataForPlot[ModelParam,10]]];
   AppendTo[data20,Join[{LogMZp//N,LogEpsM//N,LogTanT//N},DataForPlot[ModelParam,20]]];
   AppendTo[data30,Join[{LogMZp//N,LogEpsM//N,LogTanT//N},DataForPlot[ModelParam,30]]];
];
Print["Finished with LogTanT="<>ToString[LogTanT//N]];
]]

Export[NotebookDirectory[]<>"CoscatteringRegion2_x10.dat",data10];
Export[NotebookDirectory[]<>"CoscatteringRegion2_x20.dat",data20];
Export[NotebookDirectory[]<>"CoscatteringRegion2_x30.dat",data30];

]

