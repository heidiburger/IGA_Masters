close all
clear all
clc

solinit = bvpinit(linspace(0,1,40),[0,0]); %(x_start, x_end, no_points),[initial guesses for y1, y2]

sol = bvp4c(@bvp_deriv,@bvp_bcs,solinit);

plot(sol.x,sol.y(1,:),'m--');

