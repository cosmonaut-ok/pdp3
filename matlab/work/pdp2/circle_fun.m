function z = circle_fun(x,y,geometry,param)

%this function defines spatial distribution of plasma concentration
%x, y - vectors, that define 2D coordinates where concentration value
%should be calculated
%z - 2D matrix of concentration values
%x = {x1, x2, ..., xN}
%y = {y1, y2, ..., yN}
%z = {n(x1,y1) n(x2,y1) ... n(xN,y1)}
%    {n(x1,y2) n(x2,y2) ... n(xN,y2)}
%    {  ...      ...    ...    ...  }
%    {n(x1,yN) n(x2,yN) ... n(xN,yN)}

Nx = length(x);
Ny = length(y);

z = zeros(Ny,Nx);
[X Y] = meshgrid(x,y);

% center_x = param(1);
% center_y = param(2);
% radius = param(3);
% n0 = param(4);

center_x = 0.03;
center_y = 0.05;
radius = 0.01;
n0 = 1*10^11;

dx = geometry.x_size/Nx;
dy = geometry.y_size/Ny;


for i = 1:Nx
    for j = 1:Ny
        xi = i*dx;
        yi = j*dy;
        if ( (xi - center_x)^2 + (yi - center_y)^2 )^0.5 < radius
            z(j,i) = n0;
        end 
    end;
end;