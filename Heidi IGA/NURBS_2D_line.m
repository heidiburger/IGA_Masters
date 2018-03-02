%%2D NURBS code - 2D Line
%%Heidi Burger
%%Chapter 2-3, Nguyen 2012
% close all
%% Notes
%This code can draw 2D or 3D lines. 
%% Initialisation
    D=2;    %2D or 3D
%%Knot Vectors
    p=3;                                        %Polynomial order
    Xi= [0 0 0 1 1 2 2 3 3 4 4 4];        %Knot vector, Xi
    Xi_plot= (Xi(p+1:end-p+1)+Xi(p:end-p))/2;   %Knot vector, Xi, without multiplicity
%     Eta= [0 0 0 0 0.25 0.5 0.75 1 1 1 1];       %Knot vector, Eta
%     Eta_plot= (Eta(p+1:end-p+1)+Eta(p:end-p))/2;%Knot vector, Eta, without multiplicity
    n=length(Xi)-p-1;           %Number of basis functions for Xi, 1D
%     m=length(Eta)-p-1;          %Number of basis functions for Eta, 1D
    ppoints=400; %Number of points used to plot spline
%To get points at control points:
    a=min(Xi);                      %Lower bound of Xi-vector
    b=max(Xi);                      %Upper bound of Xi-vector
%     aa=min(Eta);                    %Lower bound of Eta-vector
%     bb=max(Eta);                    %Upper bound of Eta-vector
    xi_vec = linspace(a,b,ppoints);     %Parametric vector of spline, xi
%     eta_vec = linspace(aa,bb,ppoints);  %Parametric vector of spline, eta

    
    c_vec=zeros(ppoints,D);           %Initialise vector for spline
%% Control Points
    B=zeros(n,3);     %Control Points
%     const_w=zeros(n,1);     %Control Weights
    const_w=[1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1];%[1 1 1 1 1 1 1 1 1];%
    switch (D)
        case 2
            B=[1 0;1 1;0 1;-1 1;-1 0;-1 -1;0 -1;1 -1;1 0];
        case 3
            B=[1 0 0.1;1 1 0.2;0 1 0.3;-1 1 0.4;-1 0 0.5;-1 -1 0.6;0 -1 0.7;1 -1 0.8;1 0 0.9];
    end
    
%% Main
    for jj = 1:ppoints          %Full points to plot NURBS loop
        c=0;
        W=0;
        
        for i = 1:n         %Loop through basis functions
            c = c + basis_funct(p,Xi,i,xi_vec(jj))*const_w(i)*B(i,:);
        end
        for j = 1:n         %Loop through basis functions
            W = W + ( basis_funct(p,Xi,j,xi_vec(jj))*const_w(j) );
        end

        c_vec(jj,:) = c/W;    %NURBS vector
    end 
    %Define polynomial here
    func = @(x) sqrt(1 - x.^2); %3*x.^4 -16*x.^3 +8*x.^2;    %Define function here
    poly = func(c_vec(1:100,1));                    %Full points on polynomial
%% Plot Graphs

    figure(1)
%     subplot(2,1,1)
    plot(c_vec(:,1),c_vec(:,2)) %Uncomment to show full circle

%     plot(c_vec(1:100,1),c_vec(1:100,2)) %Uncomment to show quarter circle
%     hold on
%     plot(c_vec(1:100),poly,'k:')       %Plot full polynomial
%     hold off
    title('2D NURBS line')
    xlabel('X-value')
    ylabel('Y-value')
    axis equal
    if D==3
        figure(3)
        plot3(c_vec(:,1),c_vec(:,2),c_vec(:,3))
    end
    
%     figure(2)
% %     subplot(2,1,2)
%     plot(c_vec(1:100,1),c_vec(1:100,2)-poly, c_vec,zeros(size(xi_vec)))      %Plot difference between spline and polynomial, as well as the zero line. 
%     title('Error between NURBS and Polynomial')
%     ylabel('Error')
%     xlabel('x-value (c-vec)')