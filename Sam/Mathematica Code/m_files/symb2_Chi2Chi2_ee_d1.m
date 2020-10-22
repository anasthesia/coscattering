(*
     ==============================
     *  CalcHEP  3.7 *
     ==============================
  process  ~Chi2(p1)+~Chi2~(p2)->e-(p3)+e+(p4)
*)

parameters={
 MZp -> 2.00000000000*10^(1)
,MD2 -> 1.00000000000*10^(1)
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
,ctheta -> pow[1+pow[ttheta,2],-1/2]
,x28x0 -> eta*gX*sa*pow[chi,-1]*pow[ctheta,2]
,x38x0 -> ca*eta*gX*pow[chi,-1]*pow[ctheta,2]
,x62x0 -> (pow[cw,-1]*pow[sw,-1])/4.
,x62x1 -> -3*eta*sa*sw-ca*pow[cw,2]+3*ca*pow[sw,2]
,x62x2 -> -(eta*sa*sw)+ca*(pow[cw,2]+pow[sw,2])
,x65x0 -> -(pow[cw,-1]*pow[sw,-1])/4.
,x65x1 -> 3*sw*(ca*eta+sa*sw)-sa*pow[cw,2]
,x65x2 -> sw*(ca*eta+sa*sw)+sa*pow[cw,2]
,WZ -> 0.00000000000*10^(0)
,WZp -> 0.00000000000*10^(0)
           };

substitutions={
 x28x0 -> eta*gX*sa*pow[chi,-1]*pow[ctheta,2]
,x38x0 -> ca*eta*gX*pow[chi,-1]*pow[ctheta,2]
,x62x0 -> (pow[cw,-1]*pow[sw,-1])/4
,x62x1 -> -3*eta*sa*sw-ca*pow[cw,2]+3*ca*pow[sw,2]
,x62x2 -> -(eta*sa*sw)+ca*(pow[cw,2]+pow[sw,2])
,x65x0 -> -(pow[cw,-1]*pow[sw,-1])/4
,x65x1 -> 3*sw*(ca*eta+sa*sw)-sa*pow[cw,2]
,x65x2 -> sw*(ca*eta+sa*sw)+sa*pow[cw,2]
              };

inParticles = {"~Chi2","~Chi2~"}
outParticles = {"e-","e+"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =MD2^2;
p2/: SC[p2,p2] =MD2^2;
p3/: SC[p3,p3] =0^2;
p2/: SC[p2,p3] = -1*(0^2-MD2^2-MD2^2-0^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 1
*)
totFactor = ((8*x62x0^2*x28x0^2*EE^2)/(1));
numerator =(2*SC[p1,p3]^2*x62x2^2+2*SC[p1,p3]^2*x62x1^2-2*SC[p1,p3]*
 SC[p1,p2]*x62x2^2-2*SC[p1,p3]*SC[p1,p2]*x62x1^2-2*SC[p1,p3]*x62x2^2*MD2^2-
 2*SC[p1,p3]*x62x1^2*MD2^2+SC[p1,p2]^2*x62x2^2+SC[p1,p2]^2*x62x1^2+3*
 SC[p1,p2]*x62x2^2*MD2^2+3*SC[p1,p2]*x62x1^2*MD2^2+2*x62x2^2*MD2^4+2*x62x1^
 2*MD2^4);
denominator =(propDenSq[-p1-p2,MZ,0]);

addToSum;

finishSum;
