function x_vec = basis_functR1(p,Xi,n,xi_vec,const_w,B)
%Where p is the order of the basis function
%      Xi is the Knot vector
%      n is the number of basis functions
%      xi is the point at which you want to evaluate
%      const_w is the NURBS constants 

ppoints=length(xi_vec);

    x_vec=zeros(1,ppoints);
    for jj = 1:ppoints          %Full points to plot NURBS loop
        c=0;
        W=0;
        
        for j = 1:n         %Loop through basis functions
            W = W + ( basis_funct(p,Xi,j,xi_vec(jj))*const_w(j) );
        end

        for i = 1:n         %Loop through basis functions
            c = c + basis_funct(p,Xi,i,xi_vec(jj))*const_w(i)*B(i);
        end
        x_vec(jj) = c/W;    %NURBS vector
    end
   