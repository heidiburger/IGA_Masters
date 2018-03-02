clc
close all

p = 2;
Xi = [0,0,0,1,2,3,4,5,5,5];
n = length(Xi) - p - 1;

B = [1,3,3,0,2,1,4];

n_plot_points = 50;
xi_vec = linspace(0,5,n_plot_points);

C_vec = zeros(1,n_plot_points);

for jj = 1:n_plot_points
    C = 0;
    xi = xi_vec(jj);
    for i = 1:1:n
        
        C = C + basis_funct(p,Xi,i,xi)*B(i);
        
    end
    C_vec(jj) = C;
end

plot(xi_vec,C_vec)
hold on
Xi_to_plot = (Xi(p+1:end-p+1)+Xi(p:end-p))/2;
plot(Xi_to_plot,B,'ro-')