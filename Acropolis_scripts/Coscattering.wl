(* ::Package:: *)

(* ::Text:: *)
(*Initializing replacements to calculate the cross sections + SM input variables in units of MeV*)


SMMasses={MB->4200,MC->1270,MT->17300,MU->22*10^-1,MD->47*10^-1,MS->96,ME->511*10^-3,MMU->106,MTA->1777,MW->80400,MZ->91200,MH->125000,Ma->0,Mg->0,WZ->2495};
Relations={cw->Sqrt[1-swsq],sw->Sqrt[swsq],swsq->225*10^-3,aEW->10/1279,EE->2*Sqrt[aEW*Pi],gX->2*Sqrt[Pi*aX]
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
LogDgEffSFunc=Interpolation[SetPrecision[Table[{gEffSTemp[[i,1]],Log10[ND[Log[gEffSFunc[x]],x,gEffSTemp[[i,1]]]]},{i,1,Length[gEffSTemp]}],21],InterpolationOrder->SetPrecision[1,Infinity]]

Mpl=2.435*10^21;(*Mev*)
rhoR[T_]:=gEffFunc[T]*Pi^2/30*T^4/.gEffc->1075/100;
H[T_]:=SetPrecision[Sqrt[rhoR[T]/3]/Mpl,6];
Hbar[T_]:=H[T]*(1+T*10^LogDgEffSFunc[T]/3)^-1;

Entr[T_]:=2Pi^2/45 gEffc(*SFunc[T]*)*T^3/.gEffc->1075/100;

neqDM[x_]:=N[gDM*MD1^3/(2Pi^2x)BesselK[2,x]/.gDM->2//.ModelParam,21];
YeqDM[x_]:=N[neqDM[x]/Entr[MD1/x]//.ModelParam,21];

neqMed[x_]:=N[gMed*MD2^2*MD1/(2Pi^2x)BesselK[2,MD2*x/MD1]/.gMed->2//.ModelParam,21];
YeqMed[x_]:=N[neqMed[x]/Entr[MD1/x]//.ModelParam,21];

neq[x_,m_]:=m^2MD1/(2Pi^2x)BesselK[2,m*x/MD1]//.ModelParam;


(* ::Text:: *)
(*Calculates the thermally averaged cross section for a list of parameters (ModelParam) which contains only the free variables, a certain temperature MD1/T (x) and a specific process denoted by a flag (1= Chi1 SM -> Chi1 SM, 2 = Chi2 Chi2 -> SM SM, 3 = Chi1 SM -> Chi2 SM). The function uses the files saved in folder m_files, where the 0 temperature cross sections are stored only for the first generation of fermions, while assuming they are massless. *)


CalcGaScat[ModelParam_,x_,Process_]:=Module[{IntM,fileName,file,Nf,itt,IntMt,IntMi,integrand,Param},
 IntM=0;
 fileName="Int_symb"<>ToString[Process]<>"*";
 file=FileNames[fileName,"m_files"];
 Nf=Length[file];

 Param = Join[SMMasses,Relations,ModelParam];

 For[itt=1,itt<Nf+1,itt++,
  Get[file[[itt]]];
  smin=SetPrecision[smin//.Param,21];
  integrand= (E1cm^2-Mass1^2)\[Sigma] Sqrt[s]*BesselK[1,Sqrt[s]*x/MD1];
  (*Print[LogLogPlot[integrand//.Param/.{s->y*smin},{y,1,1.1}]];*)
  IntMi=NIntegrate[N[integrand*10^-50/.WD2->MD2//.Join[Param]/.{s->y*smin},21],{y,1,Infinity},WorkingPrecision->20,AccuracyGoal->20]*10^50*smin;
  IntM=IntM+IntMi;
 ];

 Return[MD1/(8 Pi^4 x)*IntM//.ModelParam];
]


(* ::Text:: *)
(*Calculates the 3 body decay width for a certain process (1 = Chi2 -> Chi1 e e)*)


Calc3BodyDecayWidth[ModelParam_,Process_]:=Module[{IntM,fileName,file,Nf,itt,IntMi,integrand,Param},
 IntM=0;
 fileName="Int3_symb"<>ToString[Process]<>"*";
 file=FileNames[fileName,"m_files"];
 Nf=Length[file];

 Param = Join[SMMasses,Relations,ModelParam];

 (*For[itt=1,itt<Nf+1,itt++,*)
 Get[file[[1]]];
 integrand= SetPrecision[MsqInt23/(2Pi)^3/(32 Mmother^3),21];

 IntM=NIntegrate[N[integrand//.Join[Param],21],{MM12,MM12min//.Param,MM12max//.Param},WorkingPrecision->20]*HeavisideTheta[MD2-MD1-2*ME]//.Param;

 Return[IntM];
]


(* ::Text:: *)
(*Calculates the thermally averaged  cross section times n1 n2 for a specific set of parameters and a specific process denoted by a flag (as above) for the temperatures stored in the list xval*)


InitRs:=Module[{data,data2,aux,finalData,rem},
data=Import[NotebookDirectory[]<>"Rratio.dat"];
data=Drop[data,4];
data=Select[data,Length[#]>10&];
rem=Select[Position[data,"the"],#[[2]]==11&];
data=Delete[data,Table[{rem[[i,1]]},{i,1,Length[rem]}]];
rem=Select[Position[data,"background"],#[[2]]==11&];
data=Delete[data,Table[{rem[[i,1]]},{i,1,Length[rem]}]];
aux=DeleteDuplicates[Table[data[[i,8]],{i,1,Length[data]}]];
data2=Table[Select[data,#[[8]]==aux[[i]]&],{i,1,Length[aux]}];
data2=Table[{data2[[i,j,1]],data2[[i,j,4]]},{i,1,Length[data2]},{j,1,Length[data2[[i]]]}];
finalData=DeleteDuplicatesBy[Join[data2[[1]],data2[[3]],data2[[7]],data2[[9]],data2[[18]],data2[[25]],data2[[33]]],First];
R[s_]:=Piecewise[{{finalData[[1,2]],s<finalData[[1,1]]},{Exp[Interpolation[Log[finalData],Log[s],InterpolationOrder->1]],finalData[[1,1]]<=s<=finalData[[-1,1]]},{finalData[[-1,2]],s>finalData[[-1,1]]}}]
]


AveragedXsec[ModelParam_,process_,xval_]:=Module[{points,svScat,x,y,i,smax,fileName,file},
  
 If[process==2||process==4,
  dataRs=Import["hadronRatio.dat"];
  R[s_]:=Piecewise[{{dataRs[[1,2]],s<dataRs[[1,1]]},{Exp[Interpolation[Log[dataRs],Log[s],InterpolationOrder->1]],dataRs[[1,1]]<=s<=dataRs[[-1,1]]},{dataRs[[-1,2]],s>dataRs[[-1,1]]}}]
 ];
 
 fileName="Int_symb*";
 file=FileNames[fileName,"m_files"];

 svScat[x_]:=CalcGaScat[Join[ModelParam],x,process];
 points={};
 For[i=1,i<=Length[xval],i++,
  y=svScat[xval[[i]]];
  (*If[y<10^-300,Break[]];*)
  points=Append[points,{xval[[i]],y}]
 ];

 Return[points];

]


(* ::Text:: *)
(*Solves the BE for a given set of parameters. It starts integrating from xmin to xmax*)


SolveCoupledBE[ModelParam_,xmin_,xmax_]:=Module[{l,dYDMinit,sol,sol1,YDMintermediate,dYDM,dYMed,YDM,YMed,consts,initCond,xInt,x2,YmedLargeX},

 l[x_]:=SetPrecision[1/(Hbar[MD1/x]*Entr[MD1/x]*x)//.ModelParam,21];

 dYDMinit[x_]:=-SetPrecision[l[x]*(CoAnnRate[x]*(YDM[x]/(YeqDM[x])-1)+(CoscatRate[x]+DecayRate[x])*(YDM[x]/(YeqDM[x])-1)),21];
 
 xInt=10;
 sol1=NDSolve[{YDM'[x]==dYDMinit[x],YDM[xmin]==YeqDM[xmin]},{YDM},{x,xmin,xInt},Method->"StiffnessSwitching",WorkingPrecision->20];
 YDMintermediate=SetPrecision[YDM[xInt]/.sol1[[1]],21];

 (*xInt=xmin;
 YDMintermediate=YeqDM[xInt];*)
 dYDM[x_]:= -SetPrecision[l[x]*(CoAnnRate[x]*(YMed[x]*YDM[x]/(YeqMed[x]*YeqDM[x])-1)+(CoscatRate[x]+DecayRate[x])*(YDM[x]/(YeqDM[x])-YMed[x]/(YeqMed[x]))),21];
 dYMed[x_]:=-SetPrecision[l[x]*(MedAnnRate[x]*(YMed[x]^2/(YeqMed[x])^2-1)+CoAnnRate[x]*(YMed[x]*YDM[x]/(YeqMed[x]*YeqDM[x])-1)-(CoscatRate[x]+DecayRate[x])*(YDM[x]/(YeqDM[x])-YMed[x]/(YeqMed[x]))),21];

 sol= NDSolve[{YDM'[x]==dYDM[x],YMed'[x]==dYMed[x],YDM[xInt]==YDMintermediate,YMed[xInt]==YeqMed[xInt]},{YMed,YDM},{x,xInt,xmax},Method->"StiffnessSwitching",WorkingPrecision->20,AccuracyGoal->35,MaxSteps->Infinity];
 (*sol= NDSolve[{YMed'[x]==dYMed[x],YMed[xInt]==YeqMed[xInt]},{YMed},{x,xInt,xmax},Method\[Rule]"StiffnessSwitching",WorkingPrecision->20,AccuracyGoal\[Rule]20];*)

 (*Ydm[x_]:=YDM[x]/.sol[[1]];
 Ymed[x_]:=YMed[x]/.sol[[1]];*)
 
 x2=SetPrecision[Max[y/.NSolve[Ga3BodyDecay==20*H[100/y],y,Reals]],6];
 YmedLargeX[x_]:=SetPrecision[Ymed[x2]*Exp[-Ga3BodyDecay/H[MD1/x]/2+Ga3BodyDecay/H[MD1/x2]/2]//.ModelParam,6];
 
 Ydm[x_]:=If[x<xInt,YDM[x]/.sol1[[1]],YDM[x]/.sol[[1]]];
 Ymed[x_]:=Which[x<xInt,YeqMed[x],x>x2,YmedLargeX[x],True,YMed[x]/.sol[[1]]];
 (*Ymed[x_]:=Which[x<xInt,YeqMed[x],True,YMed[x]/.sol[[1]]];*)

 consts={Mpl->2.435*10^21,T0->2.348*10^-10,H0->2.131*10^-39,s0->(1+21/22)*2*Pi^2/45*2*T0^3,CritDens->3*H0^2*Mpl^2 };
 Return[2*MD1*s0*(YDM[xmax]+YMed[xmax])/CritDens/.sol[[1]]//.ModelParam//.consts](*Factor of 2 to account for anti particle*);
]


SolveSingleBE[ModelParam_,xmin_,xmax_]:=Module[{l,sol,gEff,sv,Yeq,YDMintermediate,dYDM,dYMed,YDM,YMed,consts,initCond,xInt,Delta},

 l[x_]:=1/(Hbar[MD1/x]*Entr[MD1/x]*x)//.ModelParam;

 (*Factors of 2 appear due to summation of Chi and ChiBar*)
 xInt=5;
 Yeq[x_]:=2*(YeqDM[x]+YeqMed[x]);
 YDMintermediate=Yeq[xInt];
 (*Delta=(MD2-MD1)/MD1//.ModelParam;
 gEff[x_]:=4+4*(1+Delta)^(3/2)*Exp[-Delta*x];
 sv[x_]:=2*((CoAnnRate[x]/(neqMed[x]*neqDM[x]))(2/gEff[x])^2*(1+Delta)^(3/2)*Exp[-Delta*x]+(MedAnnRate[x]/neqMed[x]^2)(2/gEff[x])^2*(1+Delta)^3*Exp[-2*Delta*x]);*)

 dYDM[x_]:= -SetPrecision[l[x]*(4*CoAnnRate[x]+2*MedAnnRate[x])*((YDM[x]/Yeq[x])^2-1),21];

 sol= NDSolve[{YDM'[x]==dYDM[x],YDM[xInt]==YDMintermediate},{YDM},{x,xInt,xmax},Method->"StiffnessSwitching",WorkingPrecision->20,AccuracyGoal->10];
 (*sol= NDSolve[{YMed'[x]==dYMed[x],YMed[xInt]==YeqMed[xInt]},{YMed},{x,xInt,xmax},Method\[Rule]"StiffnessSwitching",WorkingPrecision->20,AccuracyGoal\[Rule]20];*)

 (*Ydm[x_]:=YDM[x]/.sol[[1]];
 Ymed[x_]:=YMed[x]/.sol[[1]];*)

 Ydm[x_]:=If[x<xInt,Yeq[x],YDM[x]/.sol[[1]]];

 consts={Mpl->2.435*10^21,T0->2.348*10^-10,H0->2.131*10^-39,s0->(1+21/22)*2*Pi^2/45*2*T0^3,CritDens->3*H0^2*Mpl^2 };
 Return[MD1*s0*YDM[xmax]/CritDens/.sol[[1]]//.ModelParam//.consts];
]


SolveCoscatBE[ModelParam_,xmin_,xmax_]:=Module[{l,sol,dYDM,dYMed,YDM,consts},

 l[x_]:=1/(Hbar[MD1/x]*Entr[MD1/x]*x)//.ModelParam;

 dYDM[x_]:= -SetPrecision[l[x]*(CoscatRate[x]+DecayRate[x])*(YDM[x]/(YeqDM[x])-1),21];

 sol= NDSolve[{YDM'[x]==dYDM[x],YDM[xmin]==YeqDM[xmin]},{YDM},{x,xmin,xmax},Method->"StiffnessSwitching",WorkingPrecision->20,AccuracyGoal->35,MaxSteps->Infinity];
 
 (*Ydm[x_]:=YDM[x]/.sol[[1]];
 Ymed[x_]:=YMed[x]/.sol[[1]];*)

 Ydm[x_]:=If[x<xmin,YeqDM[x],YDM[x]/.sol[[1]]];

 consts={Mpl->2.435*10^21,T0->2.348*10^-10,H0->2.131*10^-39,s0->(1+21/22)*2*Pi^2/45*2*T0^3,CritDens->3*H0^2*Mpl^2};
 Return[2*MD1*s0*YDM[xmax]/CritDens/.sol[[1]]//.ModelParam//.consts](*Factor of 2 to account for anti particle*);
]


CrossSections[ModelParam_,xmin_,xmax_]:=Module[{xval,consts,DataCoscatRate,DataCoAnnRate,DataMedAnnRate,DatasvMedAnn,DataDMScatRate},
 (*xval=Round[10^Range[Log10[xmin],Log10[xmax],(Log10[xmax]-Log10[xmin])/49]*10^5]*10^-5;50 points in the range of xmin-xmax*)
 xval=N[10^Range[Log10[xmin],Log10[xmax],(Log10[xmax]-Log10[xmin])/49],21];
 DataCoscatRate=AveragedXsec[ModelParam,3,xval];
 CoscatRateTemp[x_]:=If[x<DataCoscatRate[[-1,1]],Exp[Interpolation[Log[DataCoscatRate],Log[x]]],0];

 DataMedAnnRate=AveragedXsec[ModelParam,2,xval];
 MedAnnRateTemp[x_]:=If[x<DataMedAnnRate[[-1,1]],Exp[Interpolation[Log[DataMedAnnRate],Log[x]]],0];
 
 DataCoAnnRate=AveragedXsec[ModelParam,4,xval];
 CoAnnRateTemp[x_]:=If[x<DataCoAnnRate[[-1,1]],Exp[Interpolation[Log[DataCoAnnRate],Log[x]]],0];
 
 (*
 DataDMScatRate=AveragedXsec[ModelParam,1,xval];
 DMScatRateTemp[x_]:=If[x<DataDMScatRate[[-1,1]],Exp[Interpolation[Log[DataDMScatRate],Log[x]]],0];
 *)
 Ga3BodyDecay=Calc3BodyDecayWidth[ModelParam,1];
 DecayRateTemp[x_]:=Ga3BodyDecay*BesselK[1,x*MD2/MD1]/BesselK[2,x*MD2/MD1]*neqMed[x]//.ModelParam;
]


(* ::Text:: *)
(*Calculates the DM abundance for a specific set of model parameters*)


DMabundance[ModelParam_,xmin_,xmax_]:=Module[{xval,consts,DataCoAnnRate,DataCoscatRate,DataMedAnnRate,DatasvMedAnn,DataDMScatRate,Ga3BodyDecay},
 
 CrossSections[ModelParam,xmin,xmax];
 
 CoscatRate[x_]:=CoscatRateTemp[x];
 MedAnnRate[x_]:=MedAnnRateTemp[x];
 CoAnnRate[x_]:=CoAnnRateTemp[x];
 DecayRate[x_]:=DecayRateTemp[x];
 
 Return[SolveCoupledBE[ModelParam,xmin,xmax]];
]


DataForPlot[ModelParam_]:=Module[{xmax2,xFO,CoscatRate,xmin,xmax,xval,DataCoscatRate},
 
 xmin=1;
 (*xmax=400;
 xmin=10^-3;
 xval=Round[10^Range[Log10[xmin],Log10[xmax],0.05]*10^5]*10^-5; 

 DataCoscatRate=AveragedXsec[ModelParam,3,xval];
 CoscatRate[x_]:=If[x<DataCoscatRate[[-1,1]],Exp[Interpolation[Log[DataCoscatRate],Log[x],InterpolationOrder->1]],0];
 (*xmax2=pointlist[[3,-1,1]];*)

 If[CoscatRate[xmin]/neqDM[xmin]<H[MD1/xmin]/.ModelParam,Return[{-1,-1,-1}]];
 If[CoscatRate[xmax]/neqDM[xmax]>H[MD1/xmax]/.ModelParam,Return[{-2,-2,-2}]];
 *)

 xFO=x/.FindRoot[MedAnnRate[x]/(neqDM[x]*H[MD1/x])-1/.ModelParam,{x,xmin,100}];

 If[!(xmin<xFO<xmax),Return[{-3,-3,-3}]];

 Return[{CalcGaScat[ModelParam,xFO,1],CalcGaScat[ModelParam,xFO,3],CalcGaScat[ModelParam,xFO,4]}/(neqDM[xFO]*H[MD1/xFO])/.ModelParam]

]


script:=Module[{},

EpsMlist={1/100(*,3/100,5/100,1/10,2/10,3/10,4/10*)};
For[EpsMitt=1,EpsMitt<=Length[EpsMlist],EpsMitt++,
 ModelParam={aX->1,epsilon->10^-1,ttheta->1,MZp->300,MD2->MD1*(1+EpsMlist[[EpsMitt]]),MD1->10^2};
 CrossSections[ModelParam,1,100];
 For[LogTanT=-5,LogTanT<=-1,LogTanT=LogTanT+1/2,
  For[LogEps=-4,LogEps<=-1,LogEps=LogEps+1/2,
   
   ModelParam={aX->1,epsilon->10^LogEps,ttheta->10^LogTanT,MZp->300,MD2->MD1*(1+EpsMlist[[EpsMitt]]),MD1->10^2};
  
   CoscatRate[x_]:=CoscatRateTemp[x]*(cTheta*sTheta*epsilon)^2/(0.1/2)^2/.cTheta->1/Sqrt[ttheta^2+1]/.sTheta->ttheta/Sqrt[ttheta^2+1]/.ModelParam;
   MedAnnRate[x_]:=MedAnnRateTemp[x]*(cTheta^2*epsilon)^2/(0.1/2)^2/.cTheta->1/Sqrt[ttheta^2+1]/.sTheta->ttheta/Sqrt[ttheta^2+1]/.ModelParam;
   CoAnnRate[x_]:=CoAnnRateTemp[x]*(cTheta*sTheta*epsilon)^2/(0.1/2)^2/.cTheta->1/Sqrt[ttheta^2+1]/.sTheta->ttheta/Sqrt[ttheta^2+1]/.ModelParam;
   DecayRate[x_]:=DecayRateTemp[x]*(cTheta*sTheta*epsilon)^2/(0.1/2)^2/.cTheta->1/Sqrt[ttheta^2+1]/.sTheta->ttheta/Sqrt[ttheta^2+1]/.ModelParam;
   
   res=Join[{EpsMlist[[EpsMitt]]//N,LogEps//N,LogTanT//N,SolveCoupledBE[ModelParam,1,100],SolveCoscatBE[ModelParam,1,100],SolveSingleBE[ModelParam,1,100]}(*,DataForPlot[ModelParam]*)];
   stream=OpenAppend["MoneyPlot.dat"];
   Export[stream,{res},"Table"];
   Export[stream,"\n"];
   Close[stream];
   ];
   Print["Finished with LogTanT="<>ToString[LogTanT//N]];
  ];
 Print["Finished with LogEpsM="<>ToString[LogEpsM//N]];
 ];
];
