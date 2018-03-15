%%2D IGA scalar code
%%Heidi Burger
%%March 2018
clear all
close all

%% Notes
%

%% Setup
%System Properties
    k=1;                                %W/cm/C
    D=[k 0;0 k];
    s=@(x,y) 0;%-(x).*(4-y);%0.5.*(8-x).^2.*(4-y).^2;       %Source function
    flux=@(y) 0;%-10*y;                    %Flux function
% %Control Polygon
%     B=[-8 4, -5 4,    -2 4
%        -8 0, -5 3.25, -2 4-2*(sqrt(2)-1)
%        -8 0, -1.5 1,   2*(1-sqrt(2)) 2
%         0 0,  0 1,     0 2];             %Control points
%     
%     w=[ 1 1 1
%         1 1 (1+1/sqrt(2))/2
%         1 1 (1+1/sqrt(2))/2
%         1 1 1];                         %Matrix of NURBS weights associated with each control point. 
% %Control Polygon
%     B=[0 4, 3 4,     6 4
%        0 0, 3 3.25,  6 4-2*(sqrt(2)-1)
%        0 0, 6.5 1,   8+2*(1-sqrt(2)) 2
%        8 0, 8 1,     8 2];             %Control points
%     
%     w=[ 1 1 1
%         1 1 (1+1/sqrt(2))/2
%         1 1 (1+1/sqrt(2))/2
%         1 1 1];    

% %For plate with 1/4 hole on edge
%     B=[-4 0, -2.5 0,    -1 0
%        -4 4, -2.5 0.75, -1 sqrt(2)-1
%        -4 4, -0.75 2.5,  1-sqrt(2) 1
%         0 4,  0 2.5,     0 1];
%     w=[ 1 1 1
%         1 1 (1+1/sqrt(2))/2
%         1 1 (1+1/sqrt(2))/2
%         1 1 1];
    B=[-1 0,          -2.5 0,     -4 0
       -1 sqrt(2)-1,  -2.5 0.75,  -4 4
        1-sqrt(2) 1,  -0.75 2.5,  -4 4
        0 1,           0 2.5,      0 4];
    w=[ 1 1 1
        (1+1/sqrt(2))/2 1 1 
        (1+1/sqrt(2))/2 1 1 
        1 1 1];   
    w_temp=ones(4,3);
% %For plate with 1/4 hole on edge with 16 points in total, not 12
%     B=[-4 0, -3 0,    -2 0,  -1 0
%        -4 4, -3 2,    -2 2*(sqrt(2)-1),  -1 sqrt(2)-1
%        -4 4, -2 3,    2*(1-sqrt(2)) 2,   1-sqrt(2) 1
%         0 4,  0 3,     0 2,   0 1];
%     w=[ 1 1 1 1
%         1 1 (1+1/sqrt(2))/2 (1+1/sqrt(2))/2
%         1 1 (1+1/sqrt(2))/2 (1+1/sqrt(2))/2
%         1 1 1 1];
    size_B=size(B);
    n=size_B(1);                        %Number of basis functions in xi direction
    m=size_B(2)/2;                      %Number of basis functions in eta direction
    nn=n*m;                     %Number of nodes in mesh
    
    for j=1:m                   %Loop to convert B into a vector of pairs
        Bvec(n*(j-1)+1:n*j,[1;2])=B(1:n,[2*j-1, 2*j]);   
    end
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
    
    el_numN=n-p;
    el_numM=m-q;
    el_tot=el_numN*el_numM;
    
    Xi=[xi_s_vec linspace(xi_s,xi_e,el_numN+1) xi_e_vec];           %Knot vector (xi direction)
    Eta=[eta_s_vec linspace(eta_s,eta_e,el_numM+1) eta_e_vec];      %Knot vector (eta direction)  
    
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
    n_E=[1 2 3 4 9 10 11 12];                %Numbers of Dirichlet nodes
    d_E=[0 0 0 0 10 10 10 10]';               %Temperature values for Dirichlet nodes
    n_F=linspace(1,nn,nn);
    n_F(n_E)=[];                        %Numbers of Neumann nodes
    n_flux=[1 2];
%     for i=1:q+1
%         n_flux=[n_flux 1+(i-1)*(p+1)];
%     end
    n_flux_side=[1 2 3];
    
    el_flux=[1];                      %Not used in this code 
    
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
                    
                    [x_source, y_source]=basis_functR2vec(p,q,Xi,Eta,n,m,xi,eta,w,Bvec);        %Gauss points in physical domain
                    b=0;                                                %Source equation as a functin of x_source and y_source
                    
                    K_e = K_e + JXi_bar*det(JXi)*w_gp(gp1)*w_gp(gp2)*(dR_dxy'*D*dR_dxy);
                    F_e = F_e + JXi_bar*det(JXi)*w_gp(gp1)*w_gp(gp2)*s(x_source,y_source)*R';
                    
                end
            end
            
        %%For the Flux term
            if ismember(el,el_flux)          %Sides on which flux term appears.% Or use: ismember(el,el_flux) %elM==1
                sz = q+1;       %Number of basis functions required for each element
                yy = Bvec(asm(n_flux_side),2);    %Control points for element. (Similar to x and y above)
                for gp=1:length(Eta_bar)

                    eta = 0.5*[(Eta(j+1)-Eta(j))*Eta_bar(gp) + Eta(j+1) + Eta(j)];         %Parametric domain

                    N_dEta = zeros(1,sz);         %Derivative of Basis Functions Matrix
                    N = zeros(1,sz);             %Basis Functions Matrix
                    for nr=1:sz
                        I = nr+elM-1;
                        N_dEta(nr) = basis_funct_derivR(q,Eta,I,m,eta,w);
                        N(nr) = basis_functR(q,Eta,I,m,eta,w);
                    end

                    JEta = (N_dEta*yy);                 %Parametric jacobian
                    JEta_bar = 0.5*(Eta(j+1)-Eta(j));   %Isoparametric jacobian
                    N_dx = JEta\N_dEta';                %Basis functions wrt x
                    
                %NURBS code to determine value of Gauss point in physical domain  
                    [x_source, y_source] = basis_functR2vec(p,q,Xi,Eta,n,m,xi,eta,w,Bvec);
                %Flux equation
%                     flux=0;   %-10*y_source;        %Include x_source and y_source calculations for more complicated flux term

                    F_eN = F_eN + JEta*JEta_bar*w_gp(gp)*flux(y_source)*N';      %Neaumann contribution to Force matrix
                end
                F(asm(n_flux))=F(asm(n_flux))-F_eN(n_flux);     %Should it be +ve or -ve?
            end
        
        %Global assembly
            K(asm,asm) = K(asm,asm) + K_e;
            F(asm)     = F(asm) + F_e;
            
        end
    end
    
%% Partition
    K_E=K(n_E,n_E);
    K_EF= K(n_E,n_F);
    K_F=K(n_F,n_F);
    
    F_E = F(n_E);
    F_F = F(n_F);
    
    d_F=K_F\(F_F-K_EF'*d_E);    %Solve for unknown d nodes
    d(n_E)=d_E;
    d(n_F)=d_F;            %Control polygon for field variables 
%% Flux
    [flux_x, flux_y,flux] = basis_funct_derivR2vec(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,Bvec,d);
%% Draw
%Create temperature NURBS vector for points accross the domain. 
    [temp_vec] = basis_functR2sca(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,d);
%Create spatial NURBS vector for points accross the domain. 
    [x_vec, y_vec,c_vec] = basis_functR2vec(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,Bvec);
%% Plot
    temp_mat=reshape(temp_vec,ppointsN,ppointsM);       %Change from a vector to a matrix corresponding to ppoints. 

figure(1)   
    scatter3(c_vec(:,1),c_vec(:,2),temp_vec,10,temp_vec,'filled')
    title('3D mesh plot of temperature over 2D domain')
    xlabel('x-direction')
    ylabel('y-direction')
    zlabel('Temperature')
    axis equal
    
figure(2)
    scatter(Bvec(:,1),Bvec(:,2))
    axis equal
    
figure(3)
    hold on
    quiver(c_vec(:,1),c_vec(:,2),flux(:,1),flux(:,2),2)
%     quiver(c_vec(:,1),c_vec(:,2),flux(:,1),zeros(size(flux(:,2))),2,'r')
%     quiver(c_vec(:,1),c_vec(:,2),zeros(size(flux(:,1))),flux(:,2),2,'b')
    hold off
    axis equal
    
% vtkwrite('execute','polydata','lines',c_vec(:,1),c_vec(:,2),temp_vec)
% tri = delaunay(c_vec(:,1),c_vec(:,2));   
% vtkwrite('surface_2D_scalar.vtk','polydata','triangle',c_vec(:,1),c_vec(:,2),temp_vec,tri);