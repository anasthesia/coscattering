(*
     ==============================
     *  CalcHEP  3.7 *
     ==============================
  process  Zp(p1)->e-(p2)+e+(p3)
*)

parameters={
 MZp -> 2.00000000000*10^(1)
,swsq -> 2.25000000000*10^(-1)
,aEWM1 -> 1.27900000000*10^(2)
,MZ -> 9.11880000000*10^(1)
,Pi -> 3.14159265358*10^(0)
,cw -> pow[1-swsq,0.5]
,sw -> pow[swsq,0.5]
,aEW -> pow[aEWM1,-1]
,EE -> 2*pow[aEW,0.5]*pow[Pi,0.5]
,SignAux -> (MZ-MZp)*pow[pow[MZ-MZp,2],-0.5]
,taAux -> -1+DZ+pow[eta,2]*pow[sw,2]
,ta -> -(pow[eta,-1]*pow[sw,-1]*(taAux+SignAux*pow[4*pow[eta,2]*pow[sw,2]+
 pow[taAux,2],0.5]))/2.
,ca -> pow[1+pow[ta,2],-0.5]
,sa -> ta*pow[1+pow[ta,2],-0.5]
,eta -> chi*pow[1-pow[chi,2],-0.5]
,chi -> epsilon*pow[cw,-1]
,DZ -> (pow[MZ,-2]*pow[MZp,-2]*(-(DZaux*pow[MZ,2])+pow[MZ,4]+pow[MZp,4]+
 pow[MZp,2]*(-DZaux-2*pow[eta,2]*pow[MZ,2]*pow[sw,2])))/2.
,x65x0 -> -(pow[cw,-1]*pow[sw,-1])/4.
,x65x1 -> 3*sw*(ca*eta+sa*sw)-sa*pow[cw,2]
,x65x2 -> sw*(ca*eta+sa*sw)+sa*pow[cw,2]
           };

substitutions={
 x65x0 -> -(pow[cw,-1]*pow[sw,-1])/4
,x65x1 -> 3*sw*(ca*eta+sa*sw)-sa*pow[cw,2]
,x65x2 -> sw*(ca*eta+sa*sw)+sa*pow[cw,2]
              };

inParticles = {"Zp"}
outParticles = {"e-","e+"}

SetAttributes[ SC, Orderless ];

SC[ a_ , b_ + c_ ] := SC[a,b]+SC[a,c];

SC[ x_?NumberQ * a_ , b_ ] := x * SC[ a, b ]



p3 = +p1-p2;
p1/: SC[p1,p1] =MZp^2;
p2/: SC[p2,p2] =0^2;
p1/: SC[p1,p2] = -1*(0^2-MZp^2-0^2)/2;

initSum;

(*
  Diagram  1 in subprocess 1
*)
totFactor = ((4*x65x0^2*EE^2*MZp^2)/(3));
numerator =(x65x2^2+x65x1^2);
denominator =1;

addToSum;

finishSum;
