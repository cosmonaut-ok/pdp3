function z = homogene_fun(x,y,geometry,param)

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


n0 = param(1);

z = ones(Ny,Nx)*n0;