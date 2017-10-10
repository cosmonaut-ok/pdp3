function [ex ey] = e_from_fi(fi,geometry, bc);

%26.08
%the function returns two arrays which contains values of x- and
%y-components of electrical field;
%input parameters of function - 
%      fi - two-dimensional array of potential values;
%      geometry, bc
% 

ngx = geometry.ngx;
ngy = geometry.ngy;

dx = geometry.dx;
dy = geometry.dy;

x_type = bc.x_type;
y_type = bc.y_type;

%A - is the matrix of differentiation operator: 
%ey(i) = -(fi(i+1) - fi(i-1))/2/dy
%  (1 -1  0  0 ... 0)
%  (1  0 -1  0 ... 0)
%  (0  1  0 -1 ... 0)
%  (... ...  ... ...)
%  (0 ... 0  0  1 -1)

A = zeros(ngy);
dy_inv_half = 1/dy/2;
for i = 1:ngy-1
    A(i,i+1) = (-1)*dy_inv_half;
    A(i+1,i) = dy_inv_half;
end
A(1,1) = 1/dy;
A(1,2) = (-1)/dy;
A(end,end-1) = 1/dy;
A(end,end) = -1/dy;
ey = A*fi;

if strcmp(y_type, 'periodic')
    ey(1,:) = (fi(end,:) - fi(2,:))/2/dy;
    ey(end,:) = (fi(end-1,:) - fi(1,:))/2/dy;
end

%the same operations to calculate ex

A = zeros(ngx);

dx_inv_half = 1/dx/2;

for i = 1:ngx-1
    A(i,i+1) = (-1)*dx_inv_half;
    A(i+1,i) = dx_inv_half;
end
A(1,1) = 1/dx;
A(1,2) = (-1)/dx;
A(end,end-1) = 1/dx;
A(end,end) = (-1)/dx;
ex = (A*fi')';

if strcmp(x_type, 'periodic')
    ex(:,1) = (fi(:,end) - fi(:,2))/2/dx;
    ex(:,end) = (fi(:,end-1) - fi(:,1))/2/dx;
end