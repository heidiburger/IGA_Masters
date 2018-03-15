%%IGA_2D_scalar
%%Heidi Burger
clear all
%% Notes
%Change el_numM and el_numN, n_E and d_E, NODE, B and const_w. Set side of
%elements on which flux appear and which elements in loop. 
%% Setup
    el_numN=3;       %Number of elements per side for xi side. Depends on size of B
    el_numM=2;       %Number of elements per side for eta side.
    tot_el_num=el_numN*el_numM;   %Total number of elements in mesh (forces mesh to be square).
    xi_s=0;           %x start value
    xi_e=4;           %x end value
    eta_s=0;           %y start value
    eta_e=3;           %y end value
    
%Knot Vectors
    p=2;                            %Polynomial order (Xi)
    q=2;                            %Polynomial order (Eta)
    xi_s_vec=ones(1,p)*xi_s;
    eta_s_vec=ones(1,q)*eta_s;
    xi_e_vec=ones(1,p)*xi_e;
    eta_e_vec=ones(1,q)*eta_e;
    Xi=[xi_s_vec linspace(xi_s,xi_e,el_numN+1) xi_e_vec];         %Knot Vector %Check for non  1-to-1 mapping btwn Xi and X domain in future
    Eta=[eta_s_vec linspace(eta_s,eta_e,el_numM+1) eta_e_vec];         %Knot Vector
    
%Control Polygon
    n=length(Xi)-p-1;               %Number of basis functions for 2D (Xi)
    m=length(Eta)-q-1;              %Number of basis functions for 2D (Eta)
    nn=n*m;                         %Number of nodes in mesh
% %Control Points and Weights for plate with 1/4 hole on edge
%    B=[-4 0, -2.5 0, -1 0                    %Size of B depends of number of elements.
%        -1 sqrt(2)-1, -2.5 0.75, -4 4
%        1-sqrt(2) 1, -0.75 2.5, -4 4
%        0 1, 0 2.5, 0 4];
%    const_w=[ 1 1 1
%               (1+1/sqrt(2))/2 1 1
%               (1+1/sqrt(2))/2 1 1
%               1 1 1];
% % For 1D Heat plate (check 3.25 and -1.5 coordinate) 
% B=[-2 4, -5 4, -8 4
%    -2 4-2*(sqrt(2)-1), -5 3.25, -8 0
%    2*(1-sqrt(2)) 2, -1.5 1, -8 0
%    0 2, 0 1, 0 0];
% Bvec=[B(:,1:2); B(:,3:4); B(:,5:6)];
% const_w=[ 1 1 1
%           (1+1/sqrt(2))/2 1 1
%           (1+1/sqrt(2))/2 1 1
%           1 1 1];
%Control points and Weights for a rectangular plate. 
   B=[0 0, 0 1, 0 2, 0 3
      1 0, 1 1, 1 2, 1 3
      2 0, 2 1, 2 2, 2 3
      3 0, 3 1, 3 2, 3 3
      4 0, 4 1, 4 2, 4 3];                              %Control points
   Bvec=[B(:,1:2); B(:,3:4); B(:,5:6); B(:,7:8)];       %Make so that it manually converts B to a vector of pairs. 
   w=ones(n,m);                                   %Weights associated with control points
    
   ppoints=50;                 %Number of points used to plot NURBS per side of mesh. Change in future to ppoints on x and y side to be different
                            
%To get points at control points: 
    a=min(Xi);                      %Lower bound of Xi-vector
    b=max(Xi);                      %Upper bound of Xi-vector
    aa=min(Eta);                    %Lower bound of Eta-vector
    bb=max(Eta);                    %Upper bound of Eta-vector
    xi_vec = linspace(a,b,ppoints);     %Parametric vector of spline, xi
    eta_vec = linspace(aa,bb,ppoints);  %Parametric vector of spline, eta
    
    
%% Boundary Conditions    
    n_E=[1 2 3 4 5 16 17 18 19 20];            %Numbers of the nodes which are dirichlet nodes.
    d_E=[0 0 0 0 0 10 10 10 10 10]';            %Temperature values for bottom and top control points. 
    
%     d_E=[0 0 0 0 0 10/3 10/3 20/3 20/3 10  10  10  10  10]';   %Temperature values for all boundary control points all around mesh. 
%     n_E=[1 2 3 4 5 6 10 11 15 16 17 18 19 20];                 %Numbers of the nodes which are dirichlet nodes.
%     d_E=[0 0 0 0 10 10 10 10]';            %Temperature values for left adn right side control points. 
%     n_E=[1 6 11 16 5 10 15 20 ];          %Transpose temp_actual when using this node set. 

    n_F=linspace(1,nn,nn);
    n_F(n_E)=[];
    n_flux=[1 4 7];                             %Elemental nodes on which flux occur. For [1 4 7], the flux is on the left side of the element
    el_flux=[1 4];                          %Elements that contain flux sides. 
%% Gauss Quad and Matrix Initialisation

%For 5 pt Gauss Quadrature
    w_gp = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
    Xi_bar = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)),1/3*sqrt(5+2*sqrt(10/7))];
    Eta_bar = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)),1/3*sqrt(5+2*sqrt(10/7))];
%     Xi_bar = [sqrt(0.6) 0 -sqrt(0.6)];    %Gauss integration points/Parent Domain for 3 point Gauss
%     w = [5/9 8/9 5/9];                    %Gauss integration weights for 3 point Gauss

%Inititalise matrices   
    K=zeros(nn);                %Initialise K
    F=zeros(nn,1);              %Initialise F  
    d=zeros(nn,1);              %Initialise d
    
%% Nodes in Element
    NODE=zeros(tot_el_num,(p+1)*(q+1));                            %ICA for nodes
    c=0;
    three=0;
    for i=1:tot_el_num
        three=three+1;
        c=c+1;
        NODE(i,:)=[c, c+1, c+2, c+n, c+n+1, c+2+n, c+2*n, c+1+2*n, c+2+2*n];
        if rem(three,el_numN)==0
            c=c+2;
            three=0;
        end
    end

%% Main to construct K and F
    el=0;
    for elM=1:m-q
        for elN=1:n-p
            el=el+1;
            if p==2
                asm = NODE(el,:);           %Matrix to assemble into K for p and q=2. Called sctr in Nguyen paper. Will change for different order elements
    %         elseif p==3
    %             asm = [el el+1 el+2 el+3];
    %         elseif p==4
    %             asm = [el el+1 el+2 el+3 el+4];
            end

            K_e=0;
            F_e=0;
            F_eN=0;
            i=elN+p;                         %Equivalence of i in Nguyen paper. Indicates Knot values to use in each element.
            j=elM+q;                         %Equivalence of j in Nguyen paper. Indicates Knot values to use in each element.
            
            x=Bvec(asm,1);      %Location of elemental control points in the physical space
            y=Bvec(asm,2);
            
%             bfN=linspace(elN,elN+p+1,p+2);      %Basis functions required for this element in N direction. %[elN,elN+1,elN+2];
%             bfM=linspace(elM,elM+q+1,q+2);      %Basis functions required for this element in M direction. %[elN,elN+1,elN+2];
            
            for gp1=1:length(Xi_bar)
                for gp2=1:length(Eta_bar)

                    xi = 0.5*[(Xi(i+1)-Xi(i))*Xi_bar(gp1) + Xi(i+1) + Xi(i)];         %Parametric values of gauss point
                    eta = 0.5*[(Eta(j+1)-Eta(j))*Eta_bar(gp2) + Eta(j+1) + Eta(j)];
                    
                %Individual basis function components in each direction
                    N1=basis_functR(p,Xi,elN,n,xi,w);       %Change number of basis functions based on p and q in future. 
                    M1=basis_functR(q,Eta,elM,m,eta,w);
                    N2=basis_functR(p,Xi,elN+1,n,xi,w);
                    M2=basis_functR(q,Eta,elM+1,m,eta,w);
                    N3=basis_functR(p,Xi,elN+2,n,xi,w);
                    M3=basis_functR(q,Eta,elM+2,m,eta,w);
                    R=[N1*M1 N2*M1 N3*M1 N1*M2 N2*M2 N3*M2 N1*M3 N2*M3 N3*M3]; %Capital phi matrix
                %Individual derivative basis function components in each direction    
                    dN1_dXi=basis_funct_derivR(p,Xi,elN,n,xi,w);%basis_funct_deriv(p,Xi,elN,xi);%
                    dM1_dEta=basis_funct_derivR(q,Eta,elM,m,eta,w);
                    dN2_dXi=basis_funct_derivR(p,Xi,elN+1,n,xi,w);
                    dM2_dEta=basis_funct_derivR(q,Eta,elM+1,m,eta,w);
                    dN3_dXi=basis_funct_derivR(p,Xi,elN+2,n,xi,w);
                    dM3_dEta=basis_funct_derivR(q,Eta,elM+2,m,eta,w);
                %Derivative wrt Xi-Eta domain    
                    dR_dXi =[dN1_dXi*M1 dN2_dXi*M1 dN3_dXi*M1 dN1_dXi*M2 dN2_dXi*M2 dN3_dXi*M2 dN1_dXi*M3 dN2_dXi*M3 dN3_dXi*M3];       %Propogate automatically in a loop. 
                    dR_dEta=[N1*dM1_dEta N2*dM1_dEta N3*dM1_dEta N1*dM2_dEta N2*dM2_dEta N3*dM2_dEta N1*dM3_dEta N2*dM3_dEta N3*dM3_dEta];
                %Parametric jacobian    
                    JXi=[dR_dXi*x, dR_dEta*x;
                         dR_dXi*y, dR_dEta*y];      %Jacobian to move between the physical and parametric domain. 

                %Isoparametric jacobian     
                    JXi_bar=0.25*(Xi(i+1)-Xi(i))*(Eta(j+1)-Eta(j)); %Jacobian to move between the parametric and parent (gp) domain.
                %Derivative wrt x-y domain    
%                     [dR_dxy]=[dR_dXi' dR_dEta']/JXi;       %Same as GN' in traditional FEM
%                     dR_dx=dR_dxy(:,1);        % To separate out the x and y sides
%                     dR_dy=dR_dxy(:,2);
                    
                    dR_dxy=JXi\[dR_dXi;dR_dEta];
%                     dR_dxy=dR_dxy';
                    
                %NURBS code to determine value of Gauss point in physical domain    
                    [x_source, y_source] = basis_functR2vec(p,q,Xi,Eta,n,m,xi,eta,w,Bvec);
                %Equation for source term    
                    b=0;%*x_source;%(
%                     b=0.5*(x_source)^2*(y_source)^2;  %Use x-source and
%                     y_source as the coordinate to determine the source
%                     size at that point. 
                    
%                     N_dx=[dR_dx',zeros(1,9)
%                           zeros(1,9),dR_dy'
%                           dR_dy',dR_dx'];                   %Size:   3x18    Only for 2D vector

                    K_e = K_e + det(JXi)*det(JXi_bar)*w_gp(gp1)*w_gp(gp2)*(dR_dxy'*dR_dxy);        %Does it need two gp weights? 

                    F_e = F_e + det(JXi)*det(JXi_bar)*w_gp(gp1)*w_gp(gp2)*b*R';
                  
                    K(asm,asm) = K(asm,asm) + K_e;
                    F(asm)     = F(asm) + F_e;
                end
            end
            %% For the Flux term
            if ismember(el,el_flux)           %Sides on which flux term appears
                sz = q+1;       %Number of basis functions required for each element
                yy = B(1,2+[elM, elM+1, elM+2])';    %Control points for element. (Similar to x and y above)
                for gp=1:length(Eta_bar)

                    eta = 0.5*[(Eta(j+1)-Eta(j))*Eta_bar(gp) + Eta(j+1) + Eta(j)];         %Parametric domain

                    N_dEta = zeros(1,sz);         %Derivative of Basis Functions Matrix
                    N = zeros(1,sz);             %Basis Functions Matrix
                    for nr=1:sz
                        I = nr+elM-1;
                        N_dEta(nr) = basis_funct_derivR(q,Eta,I,m,eta,w);%basis_funct_deriv(q,Eta,I,eta);
                        N(nr) = basis_functR(q,Eta,I,m,eta,w);%basis_funct(q,Eta,I,eta);
                    end

                    JEta = (N_dEta*yy);                 %Parametric jacobian
                    JEta_bar = 0.5*(Eta(i+1)-Eta(i));   %Isoparametric jacobian
                    N_dx = JEta\N_dEta';                %Basis functions wrt x
                %NURBS code to determine value of Gauss point in physical domain  
                    [x_source, y_source] = basis_functR2vec(p,q,Xi,Eta,n,m,xi,eta,w,Bvec);
                %Flux equation
                    flux=0;   %-10*y_source;        %Include x_source and y_source calculations for more complicated flux term

                    F_eN = F_eN + JEta*JEta_bar*w_gp(gp)*flux*N';      %Neaumann contribution to Force matrix

                    F(asm(n_flux))=F(asm(n_flux))+F_eN;

                end
            end
        end
    end
%% Partitioning
    K_E=K(n_E,n_E);
    K_EF= K(n_E,n_F);
    K_F=K(n_F,n_F);
    
    F_E = F(n_E);
    F_F = F(n_F);
    
    d_F=K_F\(F_F-K_EF'*d_E);        %Solve for unknown d nodes
    d(n_E)=d_E;
    d(n_F)=d_F;            %Control polygon for field variables  
%     dmat=[d(1:el_numN+2)];
%     for i = 1:m-1
%         dmat=[dmat d( (i*(el_numN+2)+1) : (i+1)*(el_numN+2)) ];        %Convert d into a matrix of values. 
%     end
    dmat=d(1:n);
    for i = 1:m-1
        dmat=[dmat d( (i*n+1):(i+1)*n ) ];        %Convert d into a matrix of values. 
    end
%% Create temperature NURBS
    [temp_vec] = basis_functR2sca(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,d);
%% Create spatial NURBS
    [x_vec, y_vec] = basis_functR2vec(p,q,Xi,Eta,n,m,xi_vec,eta_vec,w,Bvec);
%% Plot Graphs

%     figure(1)
%     plot(c_vec(:,1),c_vec(:,2),'*')     %Check why it bunches up in the middle. 
%     axis equal
       
%     figure(2)
%     scatter(c_vec(:,1),c_vec(:,2),50,temp_vec,'*')      %Temperature plot over 2D domain 
%     xlabel('X-direction')
%     ylabel('Y-direction')
%     title('2D temperature plot with fixed boundary') 
%     axis equal       

temp_mat=reshape(temp_vec,ppoints,ppoints);
    
figure(1)
    contour(x_vec,y_vec,temp_mat',ppoints)
    title('2D contour plot of temperature distribution')
    xlabel('x-direction')
    ylabel('y-direction')
    
figure(2)   
    mesh(x_vec,y_vec,temp_mat')     %meshc(temp_mat')
    temp_actual=zeros(ppoints,ppoints);
    for i=1:ppoints
        temp_actual(i,:)=10/3*y_vec;
    end
    hold on
    mesh(x_vec,y_vec,temp_actual')
    hold off
    title('3D mesh plot of temperature over 2D domain')
    xlabel('x-direction')
    ylabel('y-direction')
    zlabel('Temperature')
    
figure(3)
    mesh(x_vec,y_vec,((temp_mat'-temp_actual')))%./temp_actual))
    title('3D mesh plot of error over 2D domain')
    xlabel('x-direction')
    ylabel('y-direction')
    zlabel('error')
 
% figure(4)
%     cx_mat=reshape(c_vec(:,1),ppoints,ppoints);
%     x_mat=cx_mat(1,:);
%     temp_fx=(10/4)*x_mat;
%     temp_mat_vec=temp_mat(:,1);
%     plot(x_mat,temp_fx)
%     hold on
%     plot(x_mat,temp_mat_vec)
%     hold off
%     title('Temperature corrolation plot between analytical and NURBS temperature value at y=0')
%     xlabel('x-direction')
%     ylabel('Temperature value')
%     legend('Analytical','NURBS values')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    