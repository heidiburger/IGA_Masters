%%2D IGA vector code
%%Heidi Burger
%%March 2018
clear all
close all

%% Notes
%

%% Setup
%Grid
    el_numN=2;                  %Number of elements per side for xi side. B depends on p and q. 
    el_numM=1;                  %Number of elements per side for eta side. 
    el_tot=el_numN*el_numM;     %Total number of elements in mesh. 
    x_s=0;                      %x start value
    x_e=4;                      %x end value
    y_s=0;                      %y start value
    y_e=4;                      %y end value
%Knot Vectors
    xi_s=0;                     %xi knot vector start value
    xi_e=1;                     %xi knot vector end value
    eta_s=0;                    %eta knot vector start value
    eta_e=1;                    %eta knot vector end value
    
    p=2;                        %Polynomial order (xi direction)
    q=2;                        %Polynomial order (eta direction)
    xi_s_vec=ones(1,p)*xi_s;    
    xi_e_vec=ones(1,p)*xi_e;
    eta_s_vec=ones(1,q)*eta_s;
    eta_e_vec=ones(1,q)*eta_e;
    
    Xi=[xi_s_vec linspace(xi_s,xi_e,el_numN+1) xi_e_vec];           %Knot vector (xi direction)
    Eta=[eta_s_vec linspace(eta_s,eta_e,el_numM+1) eta_e_vec];      %Knot vector (eta direction)
%Control Polygon
    n=length(Xi)-p-1;           %Number of basis functions in xi direction
    m=length(Eta)-q-1;          %Number of basis functions in eta direction
    nn=n*m;                     %Number of nodes in mesh
    B=zeros(n,2*m);
%     B=[0 0, 0 1, 0 2, 0 3;
%        1 0, 1 1, 1 2, 1 3;
%        2 0, 2 1, 2 2, 2 3;
%        3 0, 3 1, 3 2, 3 3;
%        4 0, 4 1, 4 2, 4 3];     %Control points.

    x_span=(x_e-x_s)/(n-1);
    y_span=(y_e-y_s)/(m-1);
    for j=1:m
        for i=1:n
            B(i,2*j-1)=x_s+x_span*(i-1);
            B(i,2*j)=y_s+y_span*(j-1);
        end
    end
    for j=1:m                   %Loop to convert B into a vector of pairs
        Bvec(n*(j-1)+1:n*j,[1;2])=B(1:n,[2*j-1, 2*j]);   
    end
    w=ones(n,m);                %Matix of NURBS weights associated with each control point

%Plotting Points
    ppointsN=50;                %Number of plotting points in xi direction
    ppointsM=50;                %Number of plotting points in eta direction
    Xi_max=max(Xi);
    Xi_min=min(Xi);
    Eta_max=max(Eta);
    Eta_min=min(Eta);
    xi_vec=linspace(Xi_min,Xi_max,ppointsN);            %Parametric vector in xi direction
    eta_vec=linspace(Eta_min,Eta_max,ppointsM);         %Parametric vector in eta direction
   
%Boundary Conditions
    s=@(x,y) 0;
    flux=@(x,y) 0;
    
    n_E= [1:n (n*m-n+1):n*m];            %Numbers of Dirichlet nodes
    d_E= [zeros(1,n)  10*ones(1,n)]';    %Temperature values for Dirichlet nodes
    n_F=linspace(1,nn,nn);
    n_F(n_E)=[];                        %Numbers of Neumann nodes
    n_flux=[];
    for i=1:q+1
        n_flux=[n_flux 1+(i-1)*(p+1)];
    end
    
    el_flux=[1 4];                      %Not used in this code

    
%% Gauss Quadrature
%For 5 point Gauss Quadrature    
    w_gp = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
    Xi_bar = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)),1/3*sqrt(5+2*sqrt(10/7))];
    Eta_bar = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)),1/3*sqrt(5+2*sqrt(10/7))];
    
%Connectivity Matrix
    NODE=zeros(el_tot,(p+1)*(q+1));     %ICA for nodes. Will change for varying p and q. 
    c=0;
    numNrow=0;
    for i=1:el_tot
        numNrow=numNrow+1;
        c=c+1;
        node=[];
        for j=1:q+1
            node=[node [(c+(j-1)*(n)):((c+(j-1)*(n))+p)]];
        end
        NODE(i,:)=node;
%         NODE(i,:)=[c, c+1, c+2, c+n, c+n+1, c+2+n, c+2*n, c+1+2*n, c+2+2*n];
        if rem(numNrow,el_numN)==0
            c=c+p;
            numNrow=0;
        end
    end
    
%% Main
%Initialise Matrices
    K=zeros(nn,nn);                     %Initialise global stiffness matrix
    F=zeros(nn,1);                      %Initialise global force matrix
    d=zeros(nn,1);                      %Initialise global displacement matrix
    
%Elemental Loop
    el=0;
    for elM=1:el_numM                   %eta direction elemental loop
        for elN=1:el_numN               %xi direction elemental loop
            el=el+1;                    %Global element number
            asm=NODE(el,:);             %Elemental connectivity matrix
            
            x=Bvec(asm,1);              %x control points associated with element
            y=Bvec(asm,2);              %y control points associated with element
            szN=p+1;                    %Number of basis functions in element in xi direction
            szM=q+1;                    %Number of basis functions in element in eta direction
            
            K_e=0;                      %Initialise elemental stiffness matrix
            F_e=0;                      %Initialise elemental source force matrix
            F_eN=0;                     %Initialise elemental flux force matrix
            i=elN+p;                    %Knot span selector in xi direction. Selects first and second (i+1) knot of element
            j=elM+q;                    %Knot span selector in eta direction. Selects first and second (i+1) knot of element
            
            for gp2=1:length(Eta_bar)           %Second Gauss Quadrature loop
                for gp1=1:length(Xi_bar)        %First Gauss Quadrature loop
                    
                    xi = 0.5*[(Xi(i+1)-Xi(i))*Xi_bar(gp1) + Xi(i+1) + Xi(i)];               %Parametric value of Gauss point (xi direction)
                    eta = 0.5*[(Eta(j+1)-Eta(j))*Eta_bar(gp2) + Eta(j+1) + Eta(j)];         %Parametric value of Gauss point (eta direction)
                    
                    N=zeros(1,szN);             %Elemental basis function matrix for xi direction
                    N_dXi=zeros(1,szN);         %Elemental basis function derivative matrix for xi direction
                    for nr=1:szN
                        I = nr+elN-1;
                        N(nr)=basis_functR(p,Xi,I,n,xi,w);
                        N_dXi(nr)=basis_funct_derivR(p,Xi,I,n,xi,w);
                    end
                    
                    M=zeros(1,szM);             %Elemental basis function matrix for eta direction
                    M_dEta=zeros(1,szM);        %Elemental basis function derivative matrix for eta direction
                    for nr=1:szM
                        I = nr+elM-1;
                        M(nr)=basis_functR(q,Eta,I,m,eta,w);
                        M_dEta(nr)=basis_funct_derivR(q,Eta,I,m,eta,w);
                    end
                    
                    nr=0;
                    R=zeros(size(asm));         %Vector of basis function products
                    dR_dXi=zeros(size(asm));    %Derivative of basis function product vector wrt to xi
                    dR_dEta=zeros(size(asm));   %Derivative of basis function product vector wrt to eta
                    for nrM=1:szM
                        for nrN=1:szN
                            nr=nr+1;
                            R(nr)=N(nrN)*M(nrM);
                            dR_dXi(nr)=N_dXi(nrN)*M(nrM);
                            dR_dEta(nr)=N(nrN)*M_dEta(nrM);
                        end
                    end
                    
                    JXi=[dR_dXi*x, dR_dEta*x;                           %Parametric Jacobian
                         dR_dXi*y, dR_dEta*y];
                    JXi_bar=0.25*(Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j));     %Isoparametric jacobian
                    
                    dR_dxy=JXi\[dR_dXi;dR_dEta];                        %Derivatives of basis functions wrt to x and y
                    dR_dx=dR_dxy(1,:);
                    dR_dy=dR_dxy(2,:);
                    
                    N_dx=[dR_dx,zeros(1,9)
                          zeros(1,9),dR_dy
                          dR_dy,dR_dx];
                    N_dx=[dR_dx(1) 0, dR_dx(2) 0, 

                    [x_source, y_source]=basis_functR2vec(p,q,Xi,Eta,n,m,xi,eta,w,Bvec);        %Gauss points in physical domain
                    b=s(x_source,y_source);                                                %Source equation as a functin of x_source and y_source
                    
                    K_e = K_e + JXi_bar*det(JXi)*w_gp(gp1)*w_gp(gp2)*(N_dx'*N_dx);
                    F_e = F_e + JXi_bar*det(JXi)*w_gp(gp1)*w_gp(gp2)*b*R';
                end
            end
            %Insert flux
        end
    end