%%2D B-Spline code
%%Heidi Burger
%%Chapter 2-3, Nguyen 2012
% close all

%% Initialisation
%%Knot Vectors
    p=2;                                        %Polynomial order
    q=2;                                        %Polynomial order
%     Xi= [0 0 0 0 0.25 0.5 0.75 1 1 1 1];        %Knot vector, Xi
%     Xi_plot= (Xi(p+1:end-p+1)+Xi(p:end-p))/2;   %Knot vector, Xi, without multiplicity
%     Eta= [0 0 0 0 0.25 0.5 0.75 1 1 1 1];       %Knot vector, Eta
%     Eta_plot= (Eta(p+1:end-p+1)+Eta(p:end-p))/2;%Knot vector, Eta, without multiplicity
    Xi=  [0 0 0 0.5 1 1 1];              %Knot vector, Xi
    Xi_plot= (Xi(p+1:end-p+1)+Xi(p:end-p))/2;   %Knot vector, Xi, without multiplicity
    Eta= [0 0 0 1 1 1];             %Knot vector, Eta (If you change this length, you need to change the control points and const_w)
    Eta_plot= (Eta(q+1:end-q+1)+Eta(q:end-q))/2;%Knot vector, Eta, without multiplicity
    n=length(Xi)-p-1;           %Number of basis functions for Xi, 1D
    m=length(Eta)-q-1;          %Number of basis functions for Eta, 1D
    ppoints=25; %Number of points used to plot spline
%To get points at control points:
    a=min(Xi);                      %Lower bound of Xi-vector
    b=max(Xi);                      %Upper bound of Xi-vector
    aa=min(Eta);                    %Lower bound of Eta-vector
    bb=max(Eta);                    %Upper bound of Eta-vector
    xi_vec = linspace(a,b,ppoints);     %Parametric vector of spline, xi
    eta_vec = linspace(aa,bb,ppoints);  %Parametric vector of spline, eta
    c_vec=zeros(ppoints,2);           %Initialise vector for spline
%%Control Points
% B=zeros(n,m);
% for i=1:n
%     for j=1:m
%         B(i,j)=sin(2*Xi(i+2))*sin(Eta(j+2));
%     end
% end
%For plate with 1/4 hole on edge
B=[-1 0, -2.5 0, -4 0
   -1 sqrt(2)-1, -2.5 0.75, -4 4
   1-sqrt(2) 1, -0.75 2.5, -4 4
   0 1, 0 2.5, 0 4];
% const_w=[ 1 1 1
%           1 1 1
%           1 1 1
%           1 1 1];
%% Main
count=0;
for ii = 1:ppoints
    for jj = 1:ppoints
        c=[0 0];
        for i = 1:n
            for j = 1:m
                N= basis_functN(p,Xi,i,xi_vec(ii));
                M= basis_functM(q,Eta,j,eta_vec(jj));

                c=c+N*M*B(i,[2*j-1 2*j]);
            end
        end
        count=count+1;
        c_vec(count,:) = c;    %NURBS vector
    end
end
%% Main
% count=0;
%     for ii = 1:ppoints
%         for jj = 1:ppoints          %Full points to plot NURBS loop
%             c=0;
%             W=0;
%             %Only evaluate for where the basis functions are non-zero. 
%             for i = 1:n         %Loop through basis functions
%                 for j = 1:m
%                     W = W + ( basis_funct(p,Xi,i,xi_vec(ii)) )*( basis_funct(q,Eta,j,eta_vec(jj)) )*const_w(i,j) ;
%                 end
%             end
% 
%             for i = 1:n         %Loop through basis functions
%                 for j =1:m
%                     c = c + ( basis_funct(p,Xi,i,xi_vec(ii)) )*( basis_funct(q,Eta,j,eta_vec(jj)) )*const_w(i,j)*B(i,[2*j-1 2*j]);
%                 end
%             end
%             count=count+1;
%             c_vec(count,:) = c;    %NURBS vector
%         end
%     end
    %%
% %%Plot graphs
% figure(1)
% %Control Polygon
% Xi_plot= Xi_plot(2:end-1);
% Eta_plot= Eta_plot(2:end-1);
% B_plot= B(2:end-1,2:end-1);
% 
% surf(Xi_plot,Eta_plot,B_plot) %,100,'m') %Use surf instead of contour3.
% alpha(0.05)
% title('Control Polygon')
% xlabel('Xi')
% ylabel('Eta')
% % figure(2)
% hold on
% %B-spline
% contour3(xi_vec,eta_vec,c_vec,50)
% title('2D B-spline')
% xlabel('Xi')
% ylabel('Eta')
% hold off

%% Plot Graphs

    figure(1)
    plot(c_vec(:,1),c_vec(:,2),'*')
    axis equal