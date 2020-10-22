\[Sigma] = (ca^2*ctheta^2*EE^2*eta^2*gX^2*stheta^2*
      (-5*ca^2*eta^2*(-1 + sw^2) - 2*ca*eta*sa*Sqrt[1 - sw^2]*(-5 + 6*sw^2) + 
       sa^2*(5 - 12*sw^2 + 8*sw^4))*
      (((MD1^2 - s)*(MD2^2 - s)*(2*MD1*MD2*MZp^2*s + 
          s*(2*MZp^4 + 3*MZp^2*s + 2*s^2 - 2*MD2^2*(MZp^2 + s)) + 
          MD1^2*(-2*s*(MZp^2 + s) + MD2^2*(MZp^2 + 2*s))))/
        (MZp^2*s*(MD1^2*(MD2^2 - s) + s*(-MD2^2 + MZp^2 + s))) + 
       (MD1^2 - 2*MD1*MD2 + MD2^2 - 2*(MZp^2 + s))*
        Log[(MD1^2*(MD2^2 - s) + s*(-MD2^2 + MZp^2 + s))/(MZp^2*s)]))/
     (4*chi^2*(MD1^2 - s)^2*sw^2*(Pi - Pi*sw^2))
smin = Max[(MD1 + ME)^2, (MD2 + ME)^2]
E1cm = (MD1^2 + s)/(2*Sqrt[s])
Mass1 = MD1
