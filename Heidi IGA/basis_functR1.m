function x_vec = basis_functR1(p,Xi,n,xi_vec,w,B)
%Where p is the order of the basis function
%      Xi is the Knot vector
%      n is the number of basis functions
%      xi_vec is the point or vector at which you want to evaluate
%      w is the NURBS constants 
%      B is the control point polygon

ppoints=length(xi_vec);

    x_vec=zeros(1,ppoints);
    for ii_points = 1:ppoints          %Full points to plot NURBS loop
        c=0;
        W=0;
        
        for ii = 1:n         %Loop through basis functions
            W = W + ( basis_funct(p,Xi,ii,xi_vec(ii_points))*w(ii) );
        end

        for ii = 1:n         %Loop through basis functions
            c = c + basis_funct(p,Xi,ii,xi_vec(ii_points))*w(ii)*B(ii);
        end
        x_vec(ii_points) = c/W;    %NURBS vector
    end
   
