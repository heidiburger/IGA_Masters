%B-Spline Basis Functions
%Ernesto 2017
clear all
% close all

figure(4)
for p = 1:4
    %p = 0;
    for n_elem = 1:3
        
        x_start = 0;
        x_end = 1;
         
        
        knot_vec_int = linspace(x_start,x_end,n_elem+1);
        
        knot_vec_front = ones(1,p)*x_start;
        knot_vec_end = ones(1,p)*x_end;
        
        knot_vec = [knot_vec_front,knot_vec_int,knot_vec_end];
        
        
        n_basis_functions = n_elem+p;
        n_pts_to_plot = 200;
        x_plot = linspace(x_start,x_end,n_pts_to_plot);
        N_plot = zeros(n_basis_functions,n_pts_to_plot);
        for I = 1:n_basis_functions
            for k = 1:n_pts_to_plot
                N_plot(I,k) =  basis_funct(p,knot_vec,I,x_plot(k));
            end
        end
        
        subplot(3,4,(n_elem-1)*4+p)
        plot(x_plot,N_plot)
        title(['p = ',num2str(p),', e = ',num2str(n_elem)])
        legend('1','2','3','4','5','6','7')
    end
end

ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0  1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.45, 0.98,'Basis Functions')
clc