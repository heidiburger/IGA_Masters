function N = basis_funct(p,Xi,i,xi)
%Where p is the order of the Shape Function
%      Xi is the Knot vector
%      i is the Shape function number
%      xi is the point at which you want to evaluate
if p == 0

   if xi >= Xi(i) && xi <= Xi(i+1)
       N = 1;
   else
       N = 0;
   end
    
else
    
    a_num = (xi - Xi(i))*basis_funct(p-1,Xi,i,xi);
    a_denom = (Xi(i+p)-Xi(i));
    if a_num == 0 && a_denom == 0
        a = 0;
    else
        a = a_num/a_denom;
    end

    b_num = (-xi + Xi(i+p+1))*basis_funct(p-1,Xi,i+1,xi);
    b_denom = (Xi(i+p+1)-Xi(i+1));
    
    if b_num == 0 && b_denom == 0
        b = 0;
    else 
        b = b_num/b_denom;
    end
    
    N = a + b;
    
end