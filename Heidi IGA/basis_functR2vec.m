function [x_vec, y_vec,c_vec] = basis_functR2vec(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,Bvec)
%Where p is the order of the basis function in the xi direction
%      q is the order of the basis functin in the eta direction
%      Xi is the Knot vector in the xi direction
%      Eta is the Knot vector in the eta direction
%      n is the number of basis functions in the xi direction
%      m is the number of basis functions in the eta direction
%      xi_vec is the point or vector at which you want to evaluate in the xi dirction
%      eta_vec is the point or vector at which you want to evaluate in the eta dirction
%      w is the NURBS constants 
%      Bvec is the control point polygon (in vector form)
%% Notes
% Control point polygon must be in vector form. 
%%

    ppointsN=length(xi_vec);
    ppointsM=length(eta_vec);

    x_vec=zeros(1,ppointsN);
    y_vec=zeros(1,ppointsM);
    c_vec=zeros(ppointsN*ppointsM,2);

    count=0;
    for jj_points = 1:ppointsM
        for ii_points = 1:ppointsN          %Full points to plot NURBS loop
            c=0;
            W=0;
            for jj = 1:m         %Loop through basis functions
                for ii = 1:n
                    W = W + ( basis_funct(p,Xi,ii,xi_vec(ii_points)) )*( basis_funct(q,Eta,jj,eta_vec(jj_points)) )*w(ii,jj);
                end
            end

            ij=0;
            for jj = 1:m         %Loop through basis functions
                for ii =1:n
                    ij=ij+1;
                    c = c + ( basis_funct(p,Xi,ii,xi_vec(ii_points)) )*( basis_funct(q,Eta,jj,eta_vec(jj_points)) )*w(ii,jj)*Bvec(ij,:);
                end
            end

            count=count+1;
            c_vec(count,:) = c/W;    %NURBS vector
            x_vec(ii_points)= c_vec(count,1);
        end
        y_vec(jj_points)= c_vec(count,2);
    end

