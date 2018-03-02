function dYdx = bvp_deriv(x,Y)

dYdx(1) = Y(2);
dYdx(2) = - x;
end