function res = bvp_bcs(ya,yb)
res = [ ya(1)-0 %dirichlet BC at start of domain. Set ya(2) for Neumann BC
    yb(1)-0];  %dirichlet BC at end of domain. Set yb(2) for Neumann BC

