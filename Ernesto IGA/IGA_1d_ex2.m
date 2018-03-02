clear all
close all
clc

p = 2;
x_start = 0;
x_end = 1;
n_elem = 32;
alpha = 50;

u_start = 0;
u_end = 0;

%no_GP = p+1; %NB, if b changing may want to "over" integrate
no_GP = 5;
if no_GP == 1
    weights = 2;
    xi_parent = 0;
elseif no_GP == 2
    weights = [1,1];
    xi_parent = [-sqrt(1/3),sqrt(1/3)];
elseif no_GP == 3
    weights = [5/9,8/9,5/9];
    xi_parent = [-sqrt(3/5),0,sqrt(3/5)];
elseif no_GP == 4
    weights = [(18-sqrt(30))/36,(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36];
    xi_parent = [-sqrt(3/7+2/7*sqrt(6/5)),-sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5))];
elseif no_GP == 5
    weights = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,(322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
    xi_parent = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)),0,1/3*sqrt(5-2*sqrt(10/7)),1/3*sqrt(5+2*sqrt(10/7))];
end



K = zeros(p+n_elem);
F = zeros(p+n_elem,1);

knot_vec_int = linspace(x_start,x_end,n_elem+1);
%knot_vec_int1 = linspace(x_start,x_start+(x_end-x_start)/2,floor(n_elem/2-p/2+2));
%knot_vec_int2 = ones(1,p-3)*x_start+(x_end-x_start)/2;
%knot_vec_int3 = linspace(x_start+(x_end-x_start)/2,x_end,ceil(n_elem/2-p/2+2));
%knot_vec_int = [knot_vec_int1,knot_vec_int2,knot_vec_int3];

knot_vec_front = ones(1,p)*x_start;
knot_vec_end = ones(1,p)*x_end;

knot_vec = [knot_vec_front,knot_vec_int,knot_vec_end];

%knot_vec = [0,0,0,0.5,1,1,1]; % knot vector. Defines parameter space
e_span = zeros(n_elem,2); %2 elements, 2 points to define span
e_support = zeros(n_elem,2); %2 elements,
for i = 1:n_elem
    e_span(i,1) = knot_vec(p+i);
    e_span(i,2) = knot_vec(p+i+1);
    
    e_support(i,1) = i;
    e_support(i,2) = p+i;
end


for e = 1:n_elem
    %k_e = zeros(p+1);
    %f_e = zeros(p+1,1);
    
    Xi_min = e_span(e,1);
    Xi_max = e_span(e,2);
    J2 = 0.5*(Xi_max-Xi_min);
    
    xi_GPs =  0.5*((Xi_max-Xi_min).*xi_parent+ (Xi_max+Xi_min));
    
    %J2 = 0.25*( Xi_max - Xi_min)
    for I = e_support(e,1):e_support(e,2)
        for J = e_support(e,1):e_support(e,2)
            
            sum = 0;
            for j = 1:no_GP
                
                sum = sum + J2*basis_funct_deriv(p,knot_vec,I,xi_GPs(j))*basis_funct_deriv(p,knot_vec,J,xi_GPs(j))*weights(j);
            end
            
            %k_e(I+1-e_support(e,1),J+1-e_support(e,1)) = sum;
            K(I,J) = K(I,J)+ sum;
            
        end
        sum1 = 0;
        for jj = 1:no_GP
            if xi_GPs(jj)<=0.58 && xi_GPs(jj)>=0.42
                b =  (2*alpha^2-4*(alpha^2*(xi_GPs(jj)-0.5))^2)*exp(-(alpha*(xi_GPs(jj)-0.5))^2);
            else
                b = 0;
            end
            b = xi_GPs(jj);
            sum1 = sum1 + J2*b*basis_funct(p,knot_vec,I,xi_GPs(jj))*weights(jj);
        end
        %f_e(I+1-e_support(e,1)) = sum1;
        F(I) =  F(I) + sum1;
    end
    
end

bcs = zeros(size(F));
bcs(1) = u_start;
bcs(end) = u_end;
F_w_BCs = K*bcs;

K(1,:) =0;
K(:,1) =0;
K(end,:) =0;
K(:,end) =0;
K(1,1) = 1;
K(end,end) = 1;
F = F-F_w_BCs;
F(1) = u_start;
F(end) = u_end;




n = length(F);

B_cv = 1*K\F;

n_plot_points = 200;
xi_vec = linspace(knot_vec(1),knot_vec(end),n_plot_points);

C_vec = zeros(1,n_plot_points);
exact_vec = zeros(1,n_plot_points);
x_vec =  zeros(1,n_plot_points);

for jj = 1:n_plot_points
    C = 0;
    C_x = 0;
    xi = xi_vec(jj);
    for i = 1:1:n
        C = C + basis_funct(p,knot_vec,i,xi)*B_cv(i);
    end
    C_vec(jj) = C;
    exact_vec(jj) = xi+exp(-(alpha*(xi-0.5))^2);
end


plot(xi_vec,C_vec) %IGA soln
hold on
%plot(xi_vec,exact_vec,'m--') %exact soln

if p == 1
    Xi_to_plot = linspace(x_start,x_end,p+n_elem);
elseif p == 2
    Xi_to_plot = [x_start,(knot_vec_int(2:end)-knot_vec_int(1:end-1))/2+knot_vec_int(1:end-1),x_end];
elseif p == 3
    if n_elem == 1
        Xi_to_plot = [x_start,x_start+(x_end-x_start)/3,x_end - (x_end-x_start)/3,x_end];
    else
        Xi_to_plot = [x_start,x_start+(knot_vec_int(3)-knot_vec_int(2))/3,knot_vec_int(2:end-1),x_end-((knot_vec_int(3)-knot_vec_int(2))/3),x_end];
    end
elseif p == 4
    if n_elem == 1
        Xi_to_plot = [x_start,x_start+(x_end-x_start)/4,x_start+(x_end-x_start)/2,x_end - (x_end-x_start)/4,x_end];
    else
        if mod(n_elem,2) == 0
            knot_spacing = (knot_vec_int(3)-knot_vec_int(2));
            l_kvi = length(knot_vec_int);
            Xi_to_plot = [x_start,x_start+knot_spacing/4,x_start+knot_spacing*3/4, knot_vec_int(2:floor(l_kvi/2))+knot_spacing/2,...
                knot_vec_int(floor(l_kvi/2)+1:end-2)+knot_spacing/2,x_end-knot_spacing*3/4,x_end-knot_spacing/4,x_end];
            
        else
            knot_spacing = (knot_vec_int(3)-knot_vec_int(2));
            Xi_to_plot = [x_start,x_start+knot_spacing/4,x_start+knot_spacing*3/4, knot_vec_int(3:end/2)-knot_spacing/2,...
                x_start+(x_end-x_start)/2,knot_vec_int(end/2+1:end-2)+knot_spacing/2,x_end-knot_spacing*3/4,x_end-knot_spacing/4,x_end];
        end
    end
end
%Xi_to_plot =[x_start,x_start+(knot_vec_int(3)-knot_vec_int(2))/p,knot_vec_int(2:end-1),x_end-((knot_vec_int(3)-knot_vec_int(2))/p),x_end]

%Xi_to_plot = linspace(x_start-(knot_vec_int(2)-knot_vec_int(1))/2,x_end+(knot_vec_int(2)-knot_vec_int(1))/2,p+n_elem)
Xi_to_plot(1) = x_start;
Xi_to_plot(end) = x_end;
plot(Xi_to_plot,B_cv,'ro-') %

%x_cv = [0,0.17,0.5,0.9,1]
%plot(x_cv,B_cv,'ro-')