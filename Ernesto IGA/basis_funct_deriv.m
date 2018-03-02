function N = basis_funct_deriv(p,Xi,i,xi)

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
    
    N = a - b;
    
end