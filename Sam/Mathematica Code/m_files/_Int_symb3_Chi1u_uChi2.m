\[Sigma] = (ctheta^2*EE^2*eta^2*gX^2*stheta^2*
      (sa^2*(-17*eta^2*sa^2*(-1 + sw^2) + 2*ca*eta*sa*Sqrt[1 - sw^2]*
          (-17 + 20*sw^2) + ca^2*(17 - 40*sw^2 + 32*sw^4))*
        (((MD1^2 - s)*(MD2^2 - s)*(2*MD1*MD2*MZ^2*s + 
            s*(2*MZ^4 + 3*MZ^2*s + 2*s^2 - 2*MD2^2*(MZ^2 + s)) + 
            MD1^2*(-2*s*(MZ^2 + s) + MD2^2*(MZ^2 + 2*s))))/
          (MZ^2*s*(MD1^2*(MD2^2 - s) + s*(-MD2^2 + MZ^2 + s))) + 
         (MD1^2 - 2*MD1*MD2 + MD2^2 - 2*(MZ^2 + s))*
          Log[(MD1^2*(MD2^2 - s) + s*(-MD2^2 + MZ^2 + s))/(MZ^2*s)]) + 
       ca^2*(2*ca*eta*sa*(17 - 20*sw^2)*Sqrt[1 - sw^2] - 
         17*ca^2*eta^2*(-1 + sw^2) + sa^2*(17 - 40*sw^2 + 32*sw^4))*
        (((MD1^2 - s)*(MD2^2 - s)*(2*MD1*MD2*MZp^2*s + 
            s*(2*MZp^4 + 3*MZp^2*s + 2*s^2 - 2*MD2^2*(MZp^2 + s)) + 
            MD1^2*(-2*s*(MZp^2 + s) + MD2^2*(MZp^2 + 2*s))))/
          (MZp^2*s*(MD1^2*(MD2^2 - s) + s*(-MD2^2 + MZp^2 + s))) + 
         (MD1^2 - 2*MD1*MD2 + MD2^2 - 2*(MZp^2 + s))*
          Log[(MD1^2*(MD2^2 - s) + s*(-MD2^2 + MZp^2 + s))/(MZp^2*s)]) - 
       2*ca*sa*(eta*sa^2*(17 - 20*sw^2)*Sqrt[1 - sw^2] + 
         ca^2*eta*Sqrt[1 - sw^2]*(-17 + 20*sw^2) + 
         ca*sa*(-17 + 40*sw^2 - 32*sw^4 - 17*eta^2*(-1 + sw^2)))*
        (((MD1^2 - s)*(-MD2^2 + s))/s + 
         ((2*MD1*MD2*MZ^2 + MZ^4 + MD1^2*(2*MD2^2 - MZ^2 - 2*s) + 2*MZ^2*s + 
             2*s^2 - MD2^2*(MZ^2 + 2*s))*Log[(MD1^2*(MD2^2 - s) + s*
                (-MD2^2 + MZ^2 + s))/(MZ^2*s)] - (2*MD1*MD2*MZp^2 + MZp^4 + 
             MD1^2*(2*MD2^2 - MZp^2 - 2*s) + 2*MZp^2*s + 2*s^2 - 
             MD2^2*(MZp^2 + 2*s))*Log[(MD1^2*(MD2^2 - s) + s*(-MD2^2 + 
                 MZp^2 + s))/(MZp^2*s)])/(MZ^2 - MZp^2))))/
     (12*chi^2*(MD1^2 - s)^2*sw^2*(Pi - Pi*sw^2))
smin = Max[(MD1 + MU)^2, (MD2 + MU)^2]
E1cm = (MD1^2 + s)/(2*Sqrt[s])
Mass1 = MD1
