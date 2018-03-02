function N = basis_funct2(p,Xi,Eta,i,xi,eta)

if p == 0
    
   if xi >= Xi(i) && xi <= Xi(i+1)
       N = 1;
   else
       N = 0;
   end
   if eta >= Eta(i) && eta <= Eta(i+1)
       M = 1;
   else
       M = 0;
   end
    
else
    
    Na_num = (xi - Xi(i))*basis_funct(p-1,Xi,i,xi);
    Na_denom = (Xi(i+p)-Xi(i));
    if Na_num == 0 && Na_denom == 0
        a = 0;
    else
        a = Na_num/Na_denom;
    end
    
    Nb_num = (-xi + Xi(i+p+1))*basis_funct(p-1,Xi,i+1,xi);
    Nb_denom = (Xi(i+p+1)-Xi(i+1));
    
    if Nb_num == 0 && Nb_denom == 0
        b = 0;
    else 
        b = Nb_num/Nb_denom;
    end
    
    N = a + b;
    
    Na_num = (xi - Xi(i))*basis_funct(p-1,Xi,i,xi);
    Na_denom = (Xi(i+p)-Xi(i));
    if Na_num == 0 && Na_denom == 0
        a = 0;
    else
        a = Na_num/Na_denom;
    end
    
    Mb_num = (-xi + Xi(i+p+1))*basis_funct(p-1,Xi,Eta,i+1,xi,eta);
    Mb_denom = (Xi(i+p+1)-Xi(i+1));
    
    if Mb_num == 0 && Mb_denom == 0
        b = 0;
    else 
        b = Mb_num/Mb_denom;
    end
    
    M = a + b;
    
end