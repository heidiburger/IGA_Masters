function [N]=shapefT(x,y,x_e,y_e)

%Shape Functions for each element
    P=[1;x;y];
    M=[1 1 1;x_e;y_e];
    N=inv(M)*P;
end