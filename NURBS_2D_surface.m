%%2D NURBS code - 2D Surface
%%Heidi Burger
%%Chapter 2-3, Nguyen 2012
% close all

%% Initialisation
%%Setup
    el_numN=2;       %Number of elements per side for xi side. Depends on size of B
    el_numM=1;       %Number of elements per side for eta side.
    tot_el_num=el_numN*el_numM;   %Total number of elements in mesh (forces mesh to be square).
    xs=0;           %x start value
    xe=1;           %x end value
    ys=0;           %y start value
    ye=1;           %y end value
    
%Knot Vectors
    p=2;                            %Polynomial order (Xi)
    q=2;                            %Polynomial order (Eta)
    xs_vec=ones(1,p)*xs;
    ys_vec=ones(1,q)*ys;
    xe_vec=ones(1,p)*xe;
    ye_vec=ones(1,q)*ye;
    Xi=[xs_vec linspace(xs,xe,el_numN+1) xe_vec];         %Knot Vector
    Eta=[ys_vec linspace(ys,ye,el_numM+1) ye_vec];         %Knot Vector
    
%     p=2;                                        %Polynomial order
%     q=2;                                        %Polynomial order
%     Xi=  [0 0 0 0.5 1 1 1];              %Knot vector, Xi
%     Xi_plot= (Xi(p+1:end-p+1)+Xi(p:end-p))/2;   %Knot vector, Xi, without multiplicity
%     Eta= [0 0 0 1 1 1];             %Knot vector, Eta (If you change this length, you need to change the control points and const_w)
%     Eta_plot= (Eta(q+1:end-q+1)+Eta(q:end-q))/2;%Knot vector, Eta, without multiplicity
    n=length(Xi)-p-1;           %Number of basis functions for Xi, 1D
    m=length(Eta)-q-1;          %Number of basis functions for Eta, 1D
    ppoints=25;%50/2; %Number of points used to plot spline
%To get points at control points:
    a=min(Xi);                      %Lower bound of Xi-vector
    b=max(Xi);                      %Upper bound of Xi-vector
    aa=min(Eta);                    %Lower bound of Eta-vector
    bb=max(Eta);                    %Upper bound of Eta-vector
    xi_vec = linspace(a,b,ppoints);     %Parametric vector of spline, xi
    eta_vec = linspace(aa,bb,ppoints);  %Parametric vector of spline, eta
    c_vec=zeros(ppoints,2);             %Initialise vector for spline
%% Control Points
%For a dougnut shape (Toroidal)
%     B=zeros(n,3*m);     %Control Points
%     const_w=zeros(n,m);     %Control Weights
%     const_w=[1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1
%              1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2)
%              1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1
%              1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2)
%              1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1
%              1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2)
%              1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1
%              1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2) 0.5 1/sqrt(2)
%              1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1 1/sqrt(2) 1];
%          
%     B=[5 0 -1, 6 0 -1, 6 0 0, 6 0 1, 5 0 1, 4 0 1, 4 0 0, 4 0 -1, 5 0 -1
%        5 5 -1, 6 6 -1, 6 6 0, 6 0 1, 5 5 1, 4 4 1, 4 4 0, 4 4 1, 5 5 1
%        0 5 -1, 0 6 -1, 0 6 0, 0 6 1, 0 5 1, 0 4 1, 0 4 0, 0 4 -1 0 5 -1
%        -5 5 -1, -6 6 -1, -6 6 0, -6 6 1, -5 5 1, -4 4 1, -4 4 0, -4 4 -1, -5 5 -1
%        -5 0 -1, -6 0 -1, -6 0 0, -6 0 1, -5 0 1, -4 0 1, -4 0 0, -4 0 -1, -5 0 -1
%        -5 -5 -1, -6 -6 -1, -6 -6 0, -6 -6 1, -5 -5 1, -4 -4 1, -4 -4 0, -4 -4 -1, -5 -5 -1
%        0 -5 -1, 0 -6 -1, 0 -6 0, 0 -6 1, 0 -5 1, 0 -4 1, 0 -4 0, 0 -4 -1, 0 -5 -1
%        5 -5 -1, 6 -6 -1, 6 -6 0, 6 -6 1, 5 -5 1, 4 -4 1, 4 -4 0, 4 -4 -1, 5 -5 -1
%        5 0 -1, 6 0 -1, 6 0 0, 6 0 1, 5 0 1, 4 0 1, 4 0 0 , 4 0 -1, 5 0 -1];     

%For plate with 1/4 hole on edge
B=[-1 0, -2.5 0, -4 0
   -1 sqrt(2)-1, -2.5 0.75, -4 4
   1-sqrt(2) 1, -0.75 2.5, -4 4
   0 1, 0 2.5, 0 4];
const_w=[ 1 1 1
          (1+1/sqrt(2))/2 1 1
          (1+1/sqrt(2))/2 1 1
          1 1 1];
% %For plate with 1/4 hole on edge
%     B=[-4 0, -2.5 0, -1 0
%        -4 4, -2.5 0.75,-1 sqrt(2)-1
%        -4 4, -0.75 2.5, 1-sqrt(2) 1
%        0 4, 0 2.5, 0 1];
%     w=[ 1 1 1
%         1 1 (1+1/sqrt(2))/2
%         1 1 (1+1/sqrt(2))/2
%         1 1 1];      
% % For 1D Heat plate (check 3.25 and -1.5 coordinate) 
% B=[-2 4, -5 4, -8 4
%    -2 4-2*(sqrt(2)-1), -5 3.25, -8 0
%    2*(1-sqrt(2)) 2, -1.5 1, -8 0
%    0 2, 0 1, 0 0];
% const_w=[ 1 1 1
%           (1+1/sqrt(2))/2 1 1
%           (1+1/sqrt(2))/2 1 1
%           1 1 1];
% % For 1D Heat plate (check 3.25 and -1.5 coordinate) 
% B=[-8 4, -5 4,    -2 4
%    -8 0, -5 3.25, -2 4-2*(sqrt(2)-1)
%    -8 0, -1.5 1,  2*(1-sqrt(2)) 2
%    0 0,  0 1,     0 2];
% const_w=[ 1 1 1
%           1 1 (1+1/sqrt(2))/2
%           1 1 (1+1/sqrt(2))/2
%           1 1 1];
% % For rectangle
%    B=[0 0, 0 1, 0 2, 0 3
%       1 0, 1 1, 1 2, 1 3
%       2 0, 2 1, 2 2, 2 3
%       3 0, 3 1, 3 2, 3 3
%       4 0, 4 1, 4 2, 4 3];
%    const_w=ones(n,m);
%         
    for j=1:m                   %Loop to convert B into a vector of pairs
        Bvec(n*(j-1)+1:n*j,[1;2])=B(1:n,[2*j-1, 2*j]);   
    end
%% Main
count=0;
    for ii = 1:ppoints
        for jj = 1:ppoints          %Full points to plot NURBS loop
            c=0;
            W=0;
            %Only evaluate for where the basis functions are non-zero. 
            for i = 1:n         %Loop through basis functions
                for j = 1:m
                    W = W + ( basis_funct(p,Xi,i,xi_vec(ii)) )*( basis_funct(q,Eta,j,eta_vec(jj)) )*const_w(i,j) ;
                end
            end

            for i = 1:n         %Loop through basis functions
                for j =1:m
                    c = c + ( basis_funct(p,Xi,i,xi_vec(ii)) )*( basis_funct(q,Eta,j,eta_vec(jj)) )*const_w(i,j)*B(i,[2*j-1 2*j]);
                end
            end
            count=count+1;
            c_vec(count,:) = c/W;    %NURBS vector
        end
    end
    x_vec=c_vec(:,1);
    y_vec=c_vec(:,2);
%     [x_vec, y_vec] = basis_functR2vec(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,Bvec);
%% Plot Graphs

    figure(1)
%     plot(c_vec(:,1),c_vec(:,2),'*')
    plot(x_vec,y_vec,'*')
    axis equal