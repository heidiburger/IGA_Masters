function [flux_x, flux_y,flux] = basis_funct_derivR2vec(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,Bvec,d)
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
%% For testing
% p=2;
% q=2;
% Xi=[0 0 0 0.5 1 1 1 ];
% Eta=[0 0 0 1 1 1 ];
% n=4;
% m=3;
% xi_vec=0.1;
% eta_vec=0.1;
% B=[-1 0,          -2.5 0,     -4 0
%        -1 sqrt(2)-1,  -2.5 0.75,  -4 4
%         1-sqrt(2) 1,  -0.75 2.5,  -4 4
%         0 1,           0 2.5,      0 4];
%     w=[ 1 1 1
%         (1+1/sqrt(2))/2 1 1 
%         (1+1/sqrt(2))/2 1 1 
%         1 1 1]; 
%     for j=1:m                   %Loop to convert B into a vector of pairs
%         Bvec(n*(j-1)+1:n*j,[1;2])=B(1:n,[2*j-1, 2*j]);   
%     end
%%
    ppointsN=length(xi_vec);
    ppointsM=length(eta_vec);

    flux_x=zeros(1,ppointsN);
    flux_y=zeros(1,ppointsM);
    flux=zeros(ppointsN*ppointsM,2);
    
    szN=n;
    szM=m;
    x=Bvec(:,1);
    y=Bvec(:,2);

    count=0;
    for jj_points = 1:ppointsM
        for ii_points = 1:ppointsN          %Full points to plot NURBS loop
            
                	N=zeros(1,szN);             %Elemental basis function matrix for xi direction
                    N_dXi=zeros(1,szN);         %Elemental basis function derivative matrix for xi direction
                    for nr=1:n
                        I = nr;
                        N(nr)=basis_functR(p,Xi,I,n,xi_vec(ii_points),w);
                        N_dXi(nr)=basis_funct_derivR(p,Xi,I,n,xi_vec(ii_points),w);
                    end
                    
                    M=zeros(1,szM);             %Elemental basis function matrix for eta direction
                    M_dEta=zeros(1,szM);        %Elemental basis function derivative matrix for eta direction
                    for nr=1:m
                        I = nr;
                        M(nr)=basis_functR(q,Eta,I,m,eta_vec(jj_points),w);
                        M_dEta(nr)=basis_funct_derivR(q,Eta,I,m,eta_vec(jj_points),w);
                    end

                    nr=0;
                    R=zeros(1,n*m);         %Vector of basis function products
                    dR_dXi=zeros(1,n*m);    %Derivative of basis function product vector wrt to xi
                    dR_dEta=zeros(1,n*m);   %Derivative of basis function product vector wrt to eta
                    for nrM=1:m
                        for nrN=1:n
                            nr=nr+1;
                            R(nr)=N(nrN)*M(nrM);
                            dR_dXi(nr)=N_dXi(nrN)*M(nrM);
                            dR_dEta(nr)=N(nrN)*M_dEta(nrM);
                        end
                    end
                    
                    JXi=[dR_dXi*x, dR_dEta*x;                           %Parametric Jacobian
                         dR_dXi*y, dR_dEta*y];
                    
                    dR_dxy=JXi\[dR_dXi;dR_dEta];

            count=count+1;
            flux(count,:) = dR_dxy*d;    %NURBS vector
            flux_x(ii_points)= flux(count,1);
        end
        flux_y(jj_points)= flux(count,2);
    end