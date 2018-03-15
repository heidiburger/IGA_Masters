%%1D IGA code
%%Heidi Burger
%November 2017
%%Chapter 2-3, Nguyen 2012
clear all
%% Notes
% Change number of elements (el_num) or polynomial order (p). Option to
% change Neumann (flux) and Dirichlet boundary conditions (select
% corresponding nodes) and source term (b) in main code. Domain defined by B and const_w

%% Setup
    el_num=2;                           %Number of elements
    funcExact = @(x) -(x.^3)/6+x/6;     %Define exact function here
    xs=0;                               %x start value
    xe=1;                               %x end value
%Knot Vectors
    p=2;                            %Polynomial order
    xs_vec=ones(1,p)*xs;
    xe_vec=ones(1,p)*xe;
    Xi=[xs_vec linspace(xs,xe,el_num+1) xe_vec];         %Knot Vector
%Control Polygon
    n=length(Xi)-p-1;               %Number of basis functions for 1D. Same as number of nodes (nn)
    B=zeros(1,n);
    B(1)=xs;
    B(end)=xe;
    B(2:end-1)=linspace(1/(p*el_num),(p*el_num-1)/(p*el_num),n-2);      %Control polygon
%     B(1:end)=linspace(xs,xe,n)        %Inaccurate B vec. Converges slower

    w=2*ones(size(B));                  %Weights associated with control points. Weights constant to convert to B-splines
    
    ppoints=50;                 %Number of points used to plot spline

%To get points at control points: 
    a=min(Xi);                  %Lower bound of Xi-vector
    b=max(Xi);                  %Upper bound of Xi-vector
    xi_vec = linspace(a,b,ppoints);     %Parametric vector for spline

%% Control Points (use to calculate from a given polynomial)
%     %Solve for Weights (B) again, given values, C which are here obtained from creating the B-spline, but would be read off a given B-spline in future.
%     C=zeros(size(Xi_plot));            
%     func = @(x) x;    %Define function here
%     poly = func(xi_vec);                    %Full points on polynomial
%     C = func(Xi_plot);                      %Sampled points on polynomial
%     C=C';
% %Weights of a given polynomial/NURBS
%     A=zeros(n);
%     for i=1:n
%         for j=1:n
%             A(i,j) = basis_funct(p,Xi,j,Xi_plot(i));        %Basis function matrix
%         end
%     end
%     B=A\C;
%     B=B';
%     B
    
%     B=[0 0.25 0.75 1];    %Control Polygon for 2 elements
%     B=[0 1/6 1/2 5/6 1];  %Control Polygon for 3 elements

%Boundary Conditions
    n_E=[1 n];           %Dirichlet nodes
    d_E=[0 0]';          %Dirichlet values
    n_F=linspace(1,n,n);
    n_F(n_E)=[];         %All other nodes
    n_flux=[p+1];          %Elemental nodes on which flux occur. For [3], the flux is on the right side of the element

%% Main for geometry to get field variables polygon

% For 5 point Gauss Quadrature    
    w_gp = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
    Xi_bar = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)),1/3*sqrt(5+2*sqrt(10/7))];
%     Xi_bar = [sqrt(0.6) 0 -sqrt(0.6)];    %Gauss integration points/Parent Domain for 3 point Gauss
%     w = [5/9 8/9 5/9];                    %Gauss integration weights for 3 point Gauss

%Initialise Matrices
    K=zeros(n);                %Initialise K
    F=zeros(n,1);              %Initialise F
%% Main    
    for el=1:el_num
        if p==2
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
            
            xi = 0.5*( (Xi(i+1)-Xi(i))*Xi_bar(gp) + Xi(i+1) + Xi(i) );         %Parametric domain
          
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
         
            b = xi;     %Equation for source term   
            
            K_e = K_e + JXi*JXi_bar*w_gp(gp)*(N_dx*N_dx');
            F_e = F_e + JXi*JXi_bar*w_gp(gp)*b*N';
            
        end

        K(asm,asm) = K(asm,asm) + K_e;
        F(asm)     = F(asm) + F_e;
    %Neumann force contribution
        if el==el_num           %Last node
            sz = length(asm);
            x = B(asm)';    %Control points for element
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

                flux=0;   %-10*y_source;        %Include x_source and y_source calculations for more complicated flux term

                F_eN = F_eN + JXi*JXi_bar*w_gp(gp)*flux*N(n_flux)';
            end
        end
        F(asm(n_flux))=F(asm(n_flux))+F_eN;
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
    
%% Draw field variable polygon
    c_vec=zeros(1,ppoints);             %Initialise vector for spline
    for jj = 1:ppoints          %Full points to plot NURBS loop
        c=0;
        W=0;
        
        for j = 1:n         %Loop through basis functions
            W = W + ( basis_funct(p,Xi,j,xi_vec(jj))*w(j) );
        end

        for i = 1:n         %Loop through basis functions
            c = c + basis_funct(p,Xi,i,xi_vec(jj))*w(i)*d(i);
        end
        c_vec(jj) = c/W;    %NURBS vector
    end
%% Draw x-vec against which to plot field variable
    x_vec=zeros(1,ppoints);             %Initialise vector for spline
    for jj = 1:ppoints          %Full points to plot NURBS loop
        c=0;
        W=0;
        
        for j = 1:n         %Loop through basis functions
            W = W + ( basis_funct(p,Xi,j,xi_vec(jj))*w(j) );
        end

        for i = 1:n         %Loop through basis functions
            c = c + basis_funct(p,Xi,i,xi_vec(jj))*w(i)*B(i);
        end
        x_vec(jj) = c/W;    %NURBS vector
    end
% %% Draw x values against which to plot the control polygon
%     d_vec=zeros(1,ppoints);
%     for jj = 1:ppoints          %Full points to plot NURBS loop
%         c=0;
%         W=0;
% %         x_vec=linspace(0,1,127);
%         for j = 1:n         %Loop through basis functions
%             W = W + ( basis_funct(p,Xi,j,xi_vec(jj))*const_w(j) );
%         end
% 
%         for i = 1:n         %Loop through basis functions
%             c = c + basis_funct(p,Xi,i,xi_vec(jj))*const_w(i)*B(i);
%         end
%         d_vec(jj) = c/W;    %NURBS vector
%     end
%% Plot graphs
    hold on
    plot(B,d,'ro-')              %Plot control polygon

    plot(x_vec,c_vec,'b--')           %Plot NURBS of field variable
    exact=funcExact(x_vec);
    plot(xi_vec,exact,'k')
    legend('Control Polygon','Nurbs','Exact')
    
    hold off
    title('1D Scalar Plot')
    xlabel('Xi parameter')
    ylabel('1D Field Variable')
