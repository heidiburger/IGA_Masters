%%IGA code
%%Heidi Burger
%%Chapter 2-3, Nguyen 2012
% close all

%% Initialisation
%%Knot Vectors
    p=2;                        %Polynomial order
    Xi= [0 0 0 1/5 2/5 3/5 4/5 1 1 1]; %[0 0 0 1 2 3 4 5 5 5];  %Knot vector
    Xi_plot = (Xi(p+1:end-p+1)+Xi(p:end-p))/2;  %Knot vector without multiplicity
    n=length(Xi)-p-1;           %Number of basis functions for 1D
    ppoints=100;             %Number of points used to plot spline
%To get points at control points:
    a=min(Xi);                  %Lower bound of Xi-vector
    b=max(Xi);                  %Upper bound of Xi-vector
    xi_vec = linspace(a,b,ppoints);     %Parametric vector for spline
    c_vec=zeros(1,ppoints);     %Initialise vector for spline
%%Control Points
%Solve for Weights (B) again, given values, C which are here obtained from creating the B-spline, but would be read off a given B-spline in future.
    C=zeros(size(Xi_plot));            
    func = @(x) sqrt(1 - x.^2); %3*x.^4 -16*x.^3 +8*x.^2;    %Define function here
    poly = func(xi_vec);                    %Full points on polynomial
    C = func(Xi_plot);                      %Sampled points on polynomial
    C=C';
%Weights of a given polynomial/B-spline
    A=zeros(n);
    for i=1:n
        for j=1:n
            A(i,j) = basis_funct(p,Xi,j,Xi_plot(i));        %Basis function matrix
        end
    end
    B=A\C;
    B=B';
    % B=[1 3 3 0 2 1 4];    %Nguyen paper spline, not matched to polynomial

%% Main
    for jj = 1:ppoints      %Full points to plot B-spline loop
        c=0;
        for i = 1:n         %Loop through basis functions
            c=c+basis_funct(p,Xi,i,xi_vec(jj))*B(i);
        end
        c_vec(jj) = c;      %B-spline vector
    end

%% Plot graphs
    figure(1)
%     subplot(2,1,1)
% Control Polygon
    plot(Xi_plot,B,'ro-')       %Plot control polygon
    hold on
    plot(xi_vec,poly,'k')       %Plot full polynomial
% Spline
    plot(xi_vec,c_vec,'--')     %Plot spline
    hold off
    title('1D B-spline')
    xlabel('Parameter')
    ylabel('1D Magnitude')
    legend('Control Polygon','Polynomial','B-spline')
    axis equal

    %% Error Analysis
% %Error plot for all B-spline points to Polynomial.
%     difference=zeros(size(xi_vec));         %Initialise difference vector
%     polygon=zeros(size(xi_vec));            %Initialise polygon vector
%     for c=1:length(xi_vec)
%         for i = 1:length(Xi_plot)-1         %Run for each line segment of control polygon
%             if (xi_vec(c)>=Xi_plot(i)) && (xi_vec(c)<Xi_plot(i+1))
%                 polygon(c) = B(i)+((B(i+1)-B(i))/(Xi_plot(i+1)-Xi_plot(i)))*(xi_vec(c)-Xi_plot(i));     %Function for control polygon segment
%                 difference(c) = c_vec(c)-polygon(c);            %Difference between NURBS and control polygon
%             end
%         end
%     end
%     
% %Fill in last entry of vectors
%     polygon(end)=B(end);
%     difference(end)=c_vec(end)-polygon(end);

    figure(2)
%     subplot(2,1,2)
    plot(xi_vec,c_vec-poly,xi_vec,zeros(size(xi_vec)))      %Plot difference between spline and polynomial, as well as the zero line. 
    title('Error between B-spline and Polynomial')
    ylabel('Error')
    xlabel('Knot Vector (xi-vec)')
    

