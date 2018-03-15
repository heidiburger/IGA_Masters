function N_dXi = basis_funct_deriv(p,Xi,i,xi)
%Where s is the start of the knot vector range for a particular element
      %p is the polynomial order
      %Xi is the knot vector
      %nr is the shape function number
      %xi is the position in the parametric domain
       
if p == 0
    N = 0;    
else
      
    a_num = p*basis_funct(p-1,Xi,i,xi);
    a_denom = (Xi(i+p)-Xi(i));
    
    if a_num == 0 && a_denom == 0
        a = 0;
    else
        a = a_num/a_denom;
    end
    
    b_num = p*basis_funct(p-1,Xi,i+1,xi);
    b_denom = (Xi(i+p+1)-Xi(i+1));
    
    if b_num == 0 && b_denom == 0
        b = 0;
    else 
        b = b_num/b_denom;
    end
    
    N_dXi = a - b;
    
%       N_dXi = (p/(Xi(s+p)-Xi(s)))*basis_funct(p-1,Xi,nr,xi) - (p/(Xi(s+p+1)-Xi(s+1)))*basis_funct(p-1,Xi,nr+1,xi);
%       a = (p/(Xi(s+p)-Xi(s)))*basis_funct(p-1,Xi,nr,xi);
%       b = (p/(Xi(s+p+1)-Xi(s+1)))*basis_funct(p-1,Xi,nr+1,xi);
%       N_dXi= a - b;
%       if N_dXi==inf
%           N_dXi=0;
%       end
end
