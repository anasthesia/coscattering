(*
     ==============================
     *  CalcHEP  3.8.4 *
     ==============================
  process  ~Chi1(p1)+Chi2~(p2)->e-(p3)+e+(p4)
*)

parameters={
 mZDinput -> 2.00000000000*10^(1)
,MD1input -> 2.00000000000*10^(1)
,MD2input -> 1.00000000000*10^(1)
,tthetainput -> -5.00000000000*10^(0)
,epsilon -> 1.00000000000*10^(-2)
,aXM1 -> 1.27900000000*10^(2)
,swsq -> 2.25000000000*10^(-1)
,aEWM1 -> 1.27900000000*10^(2)
,mZinput -> 9.11880000000*10^(1)
,Pi -> 3.14159265358*10^(0)
,cw -> pow[1-swsq,0.5]
,sw -> pow[swsq,0.5]
,aEW -> pow[aEWM1,-1]
,EE -> 2*pow[aEW,0.5]*pow[Pi,0.5]
,aX -> pow[aXM1,-1]
,gX -> 2*pow[aX,0.5]*pow[Pi,0.5]
,ta -> -(pow[eta,-1]*pow[sw,-1]*(-1+DZ+pow[eta,2]*pow[sw,2]-2*
 Atan[10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
 *(-1+DZ)]*pow[Pi,-1]*pow[4*pow[eta,2]*pow[sw,2]+pow[-1+DZ+pow[eta,2]*
 pow[sw,2],2],0.5]))/2.
,ca -> pow[1+pow[ta,2],-0.5]
,sa -> ta*pow[1+pow[ta,2],-0.5]
,eta -> epsilon*pow[cw,-1]*pow[1-pow[cw,-2]*pow[epsilon,2],-0.5]
,chi -> eta*pow[1+pow[eta,2],-0.5]
,DZ -> (pow[mZDinput,-2]*pow[mZinput,-2]*(pow[mZDinput,4]+pow[mZinput,4]+2*
 Atan[10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
 *(mZDinput-mZinput)]*pow[mZinput,2]*pow[Pi,-1]*pow[pow[mZDinput,4]+
 pow[mZinput,4]-2*pow[mZDinput,2]*pow[mZinput,2]*(1+2*pow[eta,2]*pow[sw,2])
 ,0.5]+pow[mZDinput,2]*(-2*pow[eta,2]*pow[mZinput,2]*pow[sw,2]+2*
 Atan[10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
 *(mZDinput-mZinput)]*pow[Pi,-1]*pow[pow[mZDinput,4]+pow[mZinput,4]-2*
 pow[mZDinput,2]*pow[mZinput,2]*(1+2*pow[eta,2]*pow[sw,2]),0.5])))/2.
,ctheta -> pow[1+pow[ttheta,2],-0.5]
,stheta -> ttheta*pow[1+pow[ttheta,2],-0.5]
,x26x0 -> -(ctheta*eta*gX*sa*stheta*pow[chi,-1])
,x27x0 -> -(ctheta*eta*gX*sa*stheta*pow[chi,-1])
,x36x0 -> -(ca*ctheta*eta*gX*stheta*pow[chi,-1])
,x37x0 -> -(ca*ctheta*eta*gX*stheta*pow[chi,-1])
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
x26x0 -> -(ctheta*eta*gX*sa*stheta*pow[chi,-1])
,x27x0 -> -(ctheta*eta*gX*sa*stheta*pow[chi,-1])
,x36x0 -> -(ca*ctheta*eta*gX*stheta*pow[chi,-1])
,x37x0 -> -(ca*ctheta*eta*gX*stheta*pow[chi,-1])
,x62x0 -> (pow[cw,-1]*pow[sw,-1])/4
,x62x1 -> -3*eta*sa*sw-ca*pow[cw,2]+3*ca*pow[sw,2]
,x62x2 -> -(eta*sa*sw)+ca*(pow[cw,2]+pow[sw,2])
,x65x0 -> -(pow[cw,-1]*pow[sw,-1])/4
,x65x1 -> 3*sw*(ca*eta+sa*sw)-sa*pow[cw,2]
,x65x2 -> sw*(ca*eta+sa*sw)+sa*pow[cw,2]
};

inParticles = {"~Chi1","Chi2~"}
outParticles = {"e-","e+"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p4 = +p1+p2-p3;
p1/: SC[p1,p1] =MD1^2;
p2/: SC[p2,p2] =MD2^2;
p3/: SC[p3,p3] =0^2;
p2/: SC[p2,p3] = -1*(0^2-MD1^2-MD2^2-0^2-2*SC[p1,p2]+2*SC[p1,p3])/2;

initSum;

(*
  Diagram  1 in subprocess 1
*)
totFactor = ((4*x62x0^2*x27x0*x26x0*EE^2)/(1));
numerator =(4*SC[p1,p3]^2*x62x2^2+4*SC[p1,p3]^2*x62x1^2-4*SC[p1,p3]*
 SC[p1,p2]*x62x2^2-4*SC[p1,p3]*SC[p1,p2]*x62x1^2-4*SC[p1,p3]*x62x2^2*MD1^2-
 4*SC[p1,p3]*x62x1^2*MD1^2+2*SC[p1,p2]^2*x62x2^2+2*SC[p1,p2]^2*x62x1^2+
 SC[p1,p2]*x62x2^2*MD2^2+2*SC[p1,p2]*x62x2^2*MD2*MD1+3*SC[p1,p2]*x62x2^2*
 MD1^2+SC[p1,p2]*x62x1^2*MD2^2+2*SC[p1,p2]*x62x1^2*MD2*MD1+3*SC[p1,p2]*
 x62x1^2*MD1^2+x62x2^2*MD2^3*MD1+x62x2^2*MD2^2*MD1^2+x62x2^2*MD2*MD1^3+
 x62x2^2*MD1^4+x62x1^2*MD2^3*MD1+x62x1^2*MD2^2*MD1^2+x62x1^2*MD2*MD1^3+
 x62x1^2*MD1^4);
denominator =(propDen[-p1-p2,MZ,WZ]^2);

addToSum;

(*
  Diagram  2 in subprocess 1
*)
totFactor = ((8*x65x0*x62x0*x36x0*x27x0*EE^2)/(1));
numerator =(4*SC[p1,p3]^2*x65x2*x62x2+4*SC[p1,p3]^2*x65x1*x62x1-4*SC[p1,p3]*
 SC[p1,p2]*x65x2*x62x2-4*SC[p1,p3]*SC[p1,p2]*x65x1*x62x1-4*SC[p1,p3]*x65x2*
 x62x2*MD1^2-4*SC[p1,p3]*x65x1*x62x1*MD1^2+2*SC[p1,p2]^2*x65x2*x62x2+2*
 SC[p1,p2]^2*x65x1*x62x1+SC[p1,p2]*x65x2*x62x2*MD2^2+2*SC[p1,p2]*x65x2*
 x62x2*MD2*MD1+3*SC[p1,p2]*x65x2*x62x2*MD1^2+SC[p1,p2]*x65x1*x62x1*MD2^2+2*
 SC[p1,p2]*x65x1*x62x1*MD2*MD1+3*SC[p1,p2]*x65x1*x62x1*MD1^2+x65x2*x62x2*
 MD2^3*MD1+x65x2*x62x2*MD2^2*MD1^2+x65x2*x62x2*MD2*MD1^3+x65x2*x62x2*MD1^4+
 x65x1*x62x1*MD2^3*MD1+x65x1*x62x1*MD2^2*MD1^2+x65x1*x62x1*MD2*MD1^3+x65x1*
 x62x1*MD1^4);
denominator =(propDen[-p1-p2,MZ,WZ]*propDen[-p1-p2,MZp,WZp]);

addToSum;

(*
  Diagram  3 in subprocess 1
*)
totFactor = ((4*x65x0^2*x37x0*x36x0*EE^2)/(1));
numerator =(4*SC[p1,p3]^2*x65x2^2+4*SC[p1,p3]^2*x65x1^2-4*SC[p1,p3]*
 SC[p1,p2]*x65x2^2-4*SC[p1,p3]*SC[p1,p2]*x65x1^2-4*SC[p1,p3]*x65x2^2*MD1^2-
 4*SC[p1,p3]*x65x1^2*MD1^2+2*SC[p1,p2]^2*x65x2^2+2*SC[p1,p2]^2*x65x1^2+
 SC[p1,p2]*x65x2^2*MD2^2+2*SC[p1,p2]*x65x2^2*MD2*MD1+3*SC[p1,p2]*x65x2^2*
 MD1^2+SC[p1,p2]*x65x1^2*MD2^2+2*SC[p1,p2]*x65x1^2*MD2*MD1+3*SC[p1,p2]*
 x65x1^2*MD1^2+x65x2^2*MD2^3*MD1+x65x2^2*MD2^2*MD1^2+x65x2^2*MD2*MD1^3+
 x65x2^2*MD1^4+x65x1^2*MD2^3*MD1+x65x1^2*MD2^2*MD1^2+x65x1^2*MD2*MD1^3+
 x65x1^2*MD1^4);
denominator =(propDen[-p1-p2,MZp,WZp]^2);

addToSum;

finishSum;
