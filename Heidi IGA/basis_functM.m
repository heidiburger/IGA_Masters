function M = basis_functM(q,Eta,i,eta)

if q == 0
    
   if eta >= Eta(i) && eta <= Eta(i+1)
       M = 1;
   else
       M = 0;
   end
    
else
    
    a_num = (eta - Eta(i))*basis_funct(q-1,Eta,i,eta);
    a_denom = (Eta(i+q)-Eta(i));
    if a_num == 0 && a_denom == 0
        a = 0;
    else
        a = a_num/a_denom;
    end
    
    b_num = (-eta + Eta(i+q+1))*basis_funct(q-1,Eta,i+1,eta);
    b_denom = (Eta(i+q+1)-Eta(i+1));
    
    if b_num == 0 && b_denom == 0
        b = 0;
    else 
        b = b_num/b_denom;
    end
    
    M = a + b;
    
end