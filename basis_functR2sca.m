function [temp_vec] = basis_functR2sca(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,d)
%Where p is the order of the basis function in the xi direction
%      q is the order of the basis functin in the eta direction
%      Xi is the Knot vector in the xi direction
%      Eta is the Knot vector in the eta direction
%      n is the number of basis functions in the xi direction
%      m is the number of basis functions in the eta direction
%      xi_vec is the point or vector at which you want to evaluate in the xi dirction
%      eta_vec is the point or vector at which you want to evaluate in the eta dirction
%      w is the NURBS constants 
%      d is the control point polygon (in vector form)
%% Notes
% Control point polygon must be in vector form. 
%%

    ppointsN=length(xi_vec);
    ppointsM=length(eta_vec);
    temp_vec=zeros(1,ppointsN*ppointsM);

    count=0;
    for jj_points = 1:ppointsM
        for ii_points = 1:ppointsN          %Full points to plot NURBS loop
            c=0;
            W=0;
            for jj = 1:m         %Loop through basis functions
                for ii = 1:n
                    W = W + ( basis_funct(p,Xi,ii,xi_vec(ii_points)) )*( basis_funct(q,Eta,jj,eta_vec(jj_points)) )*w(ii,jj) ;
                end
            end

            ij=0;
            for jj = 1:m         %Loop through basis functions
                for ii =1:n
                    ij=ij+1;
                    c = c + ( basis_funct(p,Xi,ii,xi_vec(ii_points)) )*( basis_funct(q,Eta,jj,eta_vec(jj_points)) )*w(ii,jj)*d(ij,:); 
                end
            end
            count=count+1;
            temp_vec(count) = c/W;    %NURBS vector
        end
    end