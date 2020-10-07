(*
     ==============================
     *  CalcHEP  3.8.4 *
     ==============================
  process  ~Chi1(p1)+e-(p2)->e-(p3)+~Chi1(p4)
*)

parameters={
 MZp -> 2.00000000000*10^(1)
,MD1 -> 1.00000000000*10^(0)
,epsilon -> 1.00000000000*10^(-6)
,aXM1 -> 1.27900000000*10^(2)
,swsq -> 2.25000000000*10^(-1)
,aEWM1 -> 1.27900000000*10^(2)
,MZ -> 9.11880000000*10^(1)
,ttheta -> 1.00000000000*10^(-5)
,ME -> 5.11000000000*10^(-4)
,Pi -> 3.14159265358*10^(0)
,cw -> pow[1-swsq,0.5]
,sw -> pow[swsq,0.5]
,aEW -> pow[aEWM1,-1]
,EE -> 2*pow[aEW,0.5]*pow[Pi,0.5]
,aX -> pow[aXM1,-1]
,gX -> 2*pow[aX,0.5]*pow[Pi,0.5]
,chi -> epsilon*pow[cw,-1]
,eta -> chi*pow[1-pow[chi,2],-0.5]
,SignAux -> (MZ-MZp)*pow[pow[MZ-MZp,2],-0.5]
,DZaux -> SignAux*pow[pow[MZ,4]+pow[MZp,4]-2*pow[MZ,2]*pow[MZp,2]*(1+2*
 pow[eta,2]*pow[sw,2]),0.5]
,DZ -> (pow[MZ,-2]*pow[MZp,-2]*(-(DZaux*pow[MZ,2])+pow[MZ,4]+pow[MZp,4]+
 pow[MZp,2]*(-DZaux-2*pow[eta,2]*pow[MZ,2]*pow[sw,2])))/2.
,taAux -> -1+DZ+pow[eta,2]*pow[sw,2]
,ta -> -(pow[eta,-1]*pow[sw,-1]*(taAux+SignAux*pow[4*pow[eta,2]*pow[sw,2]+
 pow[taAux,2],0.5]))/2.
,ca -> pow[1+pow[ta,2],-0.5]
,sa -> ta*pow[1+pow[ta,2],-0.5]
,stheta -> ttheta*pow[1+pow[ttheta,2],-0.5]
           };

substitutions={
x35x0 -> ca*eta*gX*pow[chi,-1]*pow[stheta,2]
,x65x0 -> -(pow[cw,-1]*pow[sw,-1])/4
,x65x1 -> 3*sw*(ca*eta+sa*sw)-sa*pow[cw,2]
,x65x2 -> sw*(ca*eta+sa*sw)+sa*pow[cw,2]
};

inParticles = {"~Chi1","e-"}
outParticles = {"e-","~Chi1"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =MD1^2;
p2/: SC[p2,p2] =ME^2;
p3/: SC[p3,p3] =ME^2;
p2/: SC[p2,p3] = -1*(MD1^2-MD1^2-ME^2-ME^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 1
*)
totFactor = ((8*x65x0^2*x35x0^2*EE^2)/(1));
numerator =(SC[p1,p3]^2*x65x2^2+SC[p1,p3]^2*x65x1^2-SC[p1,p3]*x65x2^2*ME^2+
 SC[p1,p3]*x65x2^2*MD1^2+SC[p1,p3]*x65x1^2*ME^2+SC[p1,p3]*x65x1^2*MD1^2+
 SC[p1,p2]^2*x65x2^2+SC[p1,p2]^2*x65x1^2+SC[p1,p2]*x65x2^2*ME^2-SC[p1,p2]*
 x65x2^2*MD1^2-SC[p1,p2]*x65x1^2*ME^2-SC[p1,p2]*x65x1^2*MD1^2-2*x65x2^2*ME^
 2*MD1^2);
denominator =(propDen[-p1+p4,MZp,0]^2);

addToSum;

finishSum;
