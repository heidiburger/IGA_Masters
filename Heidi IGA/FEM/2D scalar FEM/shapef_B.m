function [B,M]=shapef_B(x_e,y_e)

%Shape Functions for each element
    M=[1 1 1;x_e;y_e]';
    B=[0 1 0;0 0 1]/M;
end