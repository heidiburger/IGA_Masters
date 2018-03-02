clear all
% close all
% clc

p = 2;

k_e = zeros(p+1);
f_e = zeros(p+1,1);
K = zeros(4);
F = zeros(4,1);
weights = [5/9,8/9,5/9];
xi_parent = [-sqrt(3/5),0,sqrt(3/5)];

Xi = [0,0,0,0.5,1,1,1]; % knot vector. Defines parameter space
e_span = zeros(2,2); %2 elements, 2 points to define span
e_support = zeros(2,2); %2 elements,
e_span(1,1) = Xi(3);
e_span(1,2) = Xi(4);
e_support(1,1) = 1;
e_support(1,2) = 3;
e_span(2,1) = Xi(4);
e_span(2,2) = Xi(5);
e_support(2,1) = 2;
e_support(2,2) = 4;

%Xi_parent =

for e = 1:2
    Xi_min = e_span(e,1);
    Xi_max = e_span(e,2);
    if e == 1
        Xi_parent = [-1,-1,-1,1,2,2,2];
    elseif e ==2
        Xi_parent = [-2,-2,-2,-1,1,1,1];
    end
    
    for I = e_support(e,1):e_support(e,2)
        for J = e_support(e,1):e_support(e,2)
            
            sum = 0;
            for j = 1:p+1
                sum = sum + 4*4*basis_funct_deriv(2,Xi_parent,I,xi_parent(j))*basis_funct_deriv(2,Xi_parent,J,xi_parent(j))*weights(j);
            end
            
            
            k_e(I+1-e_support(e,1),J+1-e_support(e,1)) = sum;
            K(I,J) = K(I,J)+ sum;
            
        end
        sum1 = 0;
        for jj = 1:p+1
            sum1 = sum1 + 0.5*((Xi_max-Xi_min)*xi_parent(jj)+ (Xi_max+Xi_min))*basis_funct(2,Xi_parent,I,xi_parent(jj))*weights(jj);
            %xi_parent(jj)
            %0.5*((Xi_max-Xi_min)*xi_parent(jj)+ (Xi_max+Xi_min))
        end
        f_e(I+1-e_support(e,1)) = sum1;
        F(I) =  F(I) + sum1;
    end
%     k_e
%     f_e
end
% K
% F
% K\F
F_un_mod = F;
K_un_mod = K;

K(1,:) =0;
K(:,1) =0;
K(end,:) =0;
K(:,end) =0;
K(1,1) = 1;
K(end,end) = 1;
F(1) = 0;
F(end) = 0;

K\F;

n = length(Xi) - p - 1;

B_cv = 1*K\F;
% B_cv(2) = 0.046;
% B_cv(3) = 0.078;
B_x = (Xi(p+1:end-p+1)+Xi(p:end-p))/2;
%B(1) = 0;
%B(end) = 0;

n_plot_points = 50;
xi_vec = linspace(Xi(1),Xi(end),n_plot_points);

C_vec = zeros(1,n_plot_points);
exact_vec = zeros(1,n_plot_points);
x_vec =  zeros(1,n_plot_points);

for jj = 1:n_plot_points
    C = 0;
    C_x = 0;
    xi = xi_vec(jj);
    for i = 1:1:n
        
        C = C + basis_funct(p,Xi,i,xi)*B_cv(i);
        C_x = C_x + basis_funct(p,Xi,i,xi)*B_x(i);
    end
    C_vec(jj) = C;
    x_vec(jj) = C_x; %end up being xi because xi and x are aligned
    exact_vec(jj) = -xi^3/6+xi/6;
end

plot(x_vec,C_vec)
hold on
plot(xi_vec,exact_vec,'m--')
Xi_to_plot = (Xi(p+1:end-p+1)+Xi(p:end-p))/2;
plot(Xi_to_plot,B_cv,'ro-')