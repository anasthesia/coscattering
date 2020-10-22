(*
     ==============================
     *  CalcHEP  3.7 *
     ==============================
  process  ~Chi1(p1)+d(p2)->d(p3)+~Chi1(p4)
*)

parameters={
 MZp -> 2.00000000000*10^(1)
,MD1 -> 2.00000000000*10^(1)
,epsilon -> 1.00000000000*10^(-2)
,aXM1 -> 1.27900000000*10^(2)
,swsq -> 2.25000000000*10^(-1)
,aEWM1 -> 1.27900000000*10^(2)
,MZ -> 9.11880000000*10^(1)
,ttheta -> 1.00000000000*10^(-5)
,Pi -> 3.14159265358*10^(0)
,cw -> pow[1-swsq,0.5]
,sw -> pow[swsq,0.5]
,aEW -> pow[aEWM1,-1]
,EE -> 2*pow[aEW,0.5]*pow[Pi,0.5]
,aX -> pow[aXM1,-1]
,gX -> 2*pow[aX,0.5]*pow[Pi,0.5]
,ta -> -(pow[eta,-1]*pow[sw,-1]*(-1+DZ+pow[eta,2]*pow[sw,2]+pow[4*
 pow[eta,2]*pow[sw,2]+pow[-1+DZ+pow[eta,2]*pow[sw,2],2],0.5]))/2.
,ca -> pow[1+pow[ta,2],-0.5]
,sa -> ta*pow[1+pow[ta,2],-0.5]
,eta -> epsilon*pow[cw,-1]*pow[1-pow[cw,-2]*pow[epsilon,2],-0.5]
,chi -> eta*pow[1+pow[eta,2],-0.5]
,DZ -> (pow[MZ,-2]*pow[MZp,-2]*(pow[MZ,4]+pow[MZp,4]+pow[MZp,2]*(-2*
 pow[eta,2]*pow[MZ,2]*pow[sw,2]-pow[pow[MZ,4]+pow[MZp,4]-2*pow[MZ,2]*
 pow[MZp,2]*(1+2*pow[eta,2]*pow[sw,2]),0.5])-pow[MZ,2]*pow[pow[MZ,4]+
 pow[MZp,4]-2*pow[MZ,2]*pow[MZp,2]*(1+2*pow[eta,2]*pow[sw,2]),0.5]))/2.
,stheta -> ttheta*pow[1+pow[ttheta,2],-0.5]
,x25x0 -> eta*gX*sa*pow[chi,-1]*pow[stheta,2]
,x35x0 -> ca*eta*gX*pow[chi,-1]*pow[stheta,2]
,x53x0 -> (pow[cw,-1]*pow[sw,-1])/12.
,x53x1 -> -(eta*sa*sw)+ca*(-3*pow[cw,2]+pow[sw,2])
,x53x2 -> -3*eta*sa*sw+3*ca*(pow[cw,2]+pow[sw,2])
,x56x0 -> -(pow[cw,-1]*pow[sw,-1])/12.
,x56x1 -> sw*(ca*eta+sa*sw)-3*sa*pow[cw,2]
,x56x2 -> 3*(sw*(ca*eta+sa*sw)+sa*pow[cw,2])
           };

substitutions={
 x25x0 -> eta*gX*sa*pow[chi,-1]*pow[stheta,2]
,x35x0 -> ca*eta*gX*pow[chi,-1]*pow[stheta,2]
,x53x0 -> (pow[cw,-1]*pow[sw,-1])/12
,x53x1 -> -(eta*sa*sw)+ca*(-3*pow[cw,2]+pow[sw,2])
,x53x2 -> -3*eta*sa*sw+3*ca*(pow[cw,2]+pow[sw,2])
,x56x0 -> -(pow[cw,-1]*pow[sw,-1])/12
,x56x1 -> sw*(ca*eta+sa*sw)-3*sa*pow[cw,2]
,x56x2 -> 3*(sw*(ca*eta+sa*sw)+sa*pow[cw,2])
              };

inParticles = {"~Chi1","d"}
outParticles = {"d","~Chi1"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =MD1^2;
p2/: SC[p2,p2] =0^2;
p3/: SC[p3,p3] =0^2;
p2/: SC[p2,p3] = -1*(MD1^2-MD1^2-0^2-0^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 2
*)
totFactor = ((8*x53x0^2*x25x0^2*EE^2)/(1));
numerator =(SC[p1,p3]^2*x53x2^2+SC[p1,p3]^2*x53x1^2+SC[p1,p3]*x53x2^2*MD1^2+
 SC[p1,p3]*x53x1^2*MD1^2+SC[p1,p2]^2*x53x2^2+SC[p1,p2]^2*x53x1^2-SC[p1,p2]*
 x53x2^2*MD1^2-SC[p1,p2]*x53x1^2*MD1^2);
denominator =(propDen[-p1+p4,MZ,0]^2);

addToSum;

(*
  Diagram  2 in subprocess 2
*)
totFactor = ((16*x56x0*x53x0*x35x0*x25x0*EE^2)/(1));
numerator =(SC[p1,p3]^2*x56x2*x53x2+SC[p1,p3]^2*x56x1*x53x1+SC[p1,p3]*x56x2*
 x53x2*MD1^2+SC[p1,p3]*x56x1*x53x1*MD1^2+SC[p1,p2]^2*x56x2*x53x2+SC[p1,p2]^
 2*x56x1*x53x1-SC[p1,p2]*x56x2*x53x2*MD1^2-SC[p1,p2]*x56x1*x53x1*MD1^2);
denominator =(propDen[-p1+p4,MZ,0]*propDen[-p1+p4,MZp,0]);

addToSum;

(*
  Diagram  3 in subprocess 2
*)
totFactor = ((8*x56x0^2*x35x0^2*EE^2)/(1));
numerator =(SC[p1,p3]^2*x56x2^2+SC[p1,p3]^2*x56x1^2+SC[p1,p3]*x56x2^2*MD1^2+
 SC[p1,p3]*x56x1^2*MD1^2+SC[p1,p2]^2*x56x2^2+SC[p1,p2]^2*x56x1^2-SC[p1,p2]*
 x56x2^2*MD1^2-SC[p1,p2]*x56x1^2*MD1^2);
denominator =(propDen[-p1+p4,MZp,0]^2);

addToSum;

finishSum;
