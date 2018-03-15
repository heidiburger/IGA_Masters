function dR_dXi = basis_funct_derivR(p,Xi,I,n,xi,w)
%Where p is the order of the basis function
%      Xi is the Knot vector
%      I is the basis functin number
%      n is the number of basis functions
%      xi_vec is the point or vector at which you want to evaluate
%      w is the NURBS constants 

        Wdash=0;
        W=0;
        
        for ii = 1:n         %Loop through basis functions
            W = W + ( basis_funct(p,Xi,ii,xi)*w(ii) );
        end
        
        for ii = 1:n         %Loop through basis functions
            Wdash = Wdash + ( basis_funct_deriv(p,Xi,ii,xi)*w(ii) );
        end
        
        
        

        dR_dXi=w(I)*(basis_funct_deriv(p,Xi,I,xi)*W-basis_funct(p,Xi,I,xi)*Wdash)/W^2;