function [resx resy] = load_linear_d(geometry, concentration,randx,randy)

%the function assigns spatial coorditanes to particles;
%it uses 2xNbig values between 0 and 1 - randx, randy
%geometry and spatial concentration distribution - concentration

%the function works only with a linear spatial distribution

n0 = concentration.param(1);
dn = concentration.param(2) - concentration.param(1);

x_size = geometry.x_size;
y_size = geometry.y_size;

resx = x_size/dn*((n0^2 + randx*(2*n0*dn + dn^2)).^0.5 - n0);
resy = y_size*randy;
