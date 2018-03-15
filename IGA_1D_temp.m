%%1D IGA temperature code in radial direction
%%Heidi Burger
%February 2018
clear all
close all

%% Notes
% Change number of elements (el_num) or polynomial order (p). Option to
% change Neumann (flux) and Dirichlet boundary conditions (select
% corresponding nodes) and source term (b) in main code. Domain defined by B and const_w

%% SETUP
    s=-101;                               %Source term
    k=10;                               %Conductance
    x_start=0;                         %Start of radius
    R=2;                               %Radius of plate
    
    el_num=3;                           %Number of elements
%     funcExact = @(r) s/(4*k)*(R^2 - r.^2);     %Define exact function here. can only be used for terms where x_start=0
    funcExact = @(r) x_start^2/2*s/k*log(r)+R.^2/4*s/k-x_start.^2/2*s/k*log(R)-r.^2/4*s/k;      %More general Exact solution. Works for Doughnut shapes. 
    
    xi_s=0;                               %Knot vector start value
    xi_e=1;                               %Knot vector end value
%Knot Vectors
    p=2;                                %Polynomial order
    xi_s_vec=ones(1,p)*xi_s;
    xi_e_vec=ones(1,p)*xi_e;
    Xi=[xi_s_vec linspace(xi_s,xi_e,el_num+1) xi_e_vec];                        %Knot Vector
%Control Polygon
    n=length(Xi)-p-1;                   %Number of basis functions for 1D. Same as number of nodes (nn)
    B=zeros(1,n);
    B(1)=x_start;
    B(end)=R;
    B = linspace(x_start, R,n);
%     B(2:end-1)=linspace(x_start+1/(p*el_num),R-(p*el_num-1)/(p*el_num),n-2);      %Control polygon
    w=2*ones(size(B));            %Weights associated with control points. Weights constant to convert to B-splines
    
    ppoints=50;                         %Number of points used to plot spline

%To get points at control points: 
    a=min(Xi);                          %Lower bound of Xi-vector
    b=max(Xi);                          %Upper bound of Xi-vector
    xi_vec = linspace(a,b,ppoints);     %Parametric vector for spline
    
%Boundary Conditions
    n_E=[ n];           %Dirichlet nodes
    d_E=[ 0]';          %Dirichlet values
    n_F=linspace(1,n,n);
    n_F(n_E)=[];         %All other nodes
    n_flux=[3];          %Elemental nodes on which flux occur. For [3], the flux is on the right side of the element

%% Gauss Quadrature
%For 5 point Gauss Quadrature    
    w_gp = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
    Xi_bar = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)),1/3*sqrt(5+2*sqrt(10/7))];

%% MAIN   
%Initialise Matrices
    K=zeros(n);                %Initialise K
    F=zeros(n,1);              %Initialise F

%Main Loops
    for el=1:el_num
        if p==2             %Length of asm is of length p+1
        	asm = [el el+1 el+2];           %Matrix to assemble into global K. Called sctr in Nguyen paper. Will change for different order elements
        elseif p==3
            asm = [el el+1 el+2 el+3];
        elseif p==4
            asm = [el el+1 el+2 el+3 el+4];
        end
        
        sz = length(asm);   %Number of basis functions in element
        x = B(asm)';    %Control points for element
        
        K_e=0;                          %Initialise elemental K
        F_e=0;                          %Initialise elemental F
        F_eN=0;                         %Initialise elemental flux F
        i=el+p;                         %Equivalence of i in Nguyen paper. Indicates Knot values to use in each element.
        
        for gp=1:length(Xi_bar)
            
            xi = 0.5*[(Xi(i+1)-Xi(i))*Xi_bar(gp) + Xi(i+1) + Xi(i)];         %Parametric domain

            N_dXi = zeros(1,sz);         %Derivative of Basis Functions Matrix
            N = zeros(1,sz);             %Basis Functions Matrix
            for nr=1:sz
                I = nr+el-1;
                N_dXi(nr) = basis_funct_deriv(p,Xi,I,xi);
                N(nr) = basis_funct(p,Xi,I,xi);
            end

            JXi = (N_dXi*x);                %Parametric jacobian
            
            JXi_bar = 0.5*(Xi(i+1)-Xi(i));  %Isoparametric jacobian
           
            N_dx = JXi\N_dXi';              %Derivatives of basis functions wrt x

        %NURBS code to determine value of Gauss point in physical domain  
            r=basis_functR1(p,Xi,n,xi,w,B);

            K_e = K_e + JXi*JXi_bar*w_gp(gp)*k*r*(N_dx*N_dx');
            F_e = F_e + JXi*JXi_bar*w_gp(gp)*s*r*N';
            
        end

        K(asm,asm) = K(asm,asm) + K_e;
        F(asm)     = F(asm) + F_e;
% No Neumann condition in this question.
    end 
    
%% Solve using partition method
    K_E=K(n_E,n_E);
    K_EF= K(n_E,n_F);
    K_F=K(n_F,n_F);
    
    F_E = F(n_E);
    F_F = F(n_F);
    
    d_F=K_F\(F_F-K_EF'*d_E);    %Solve for unknown d nodes
    d(n_E)=d_E;
    d(n_F)=d_F;            %Control polygon for field variables          
    
%% DRAW
%Draw field variable vector
    c_vec=basis_functR1(p,Xi,n,xi_vec,w,d);    
%Draw x-vec against which to plot field variable
    x_vec=basis_functR1(p,Xi,n,xi_vec,w,B);
    
%% PLOT
figure(1)
    hold on
    plot(B,d,'ro-')              %Plot control polygon

    
    exact=funcExact(x_vec);
    plot(x_vec,exact,'k')
    plot(x_vec,c_vec,'b--')           %Plot NURBS of field variable
    legend('Control Polygon','Exact','Nurbs')
    
    hold off
    title('1D Scalar Plot')
    xlabel('r - Physical domain')
    ylabel('1D Field Variable')  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    