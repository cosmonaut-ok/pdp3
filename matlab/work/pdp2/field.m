function [ex ey fi] = field(rho, bcx_type, bcy_type, bcl, bcr, bcu, bcd, dx, dy)

%  10.01.2006
%  the function solves discrete 2D Poisson's equation;
%  rho - 2D matrix with charge density values, [Cl*m^-3]
%


%
%  [ngy ngx] = size(rho) - size of two-dimensional 
%                          rectangular grid(including borders)
%  dx, dy - length and width of elementary grid cell, [m]
%
%  bcx_type, bcy_type - type of boundary conditions on
%                         left, right(x_bc_type) and
%                         upper, lower(y_bc_type) boundaries
%  boundary conditions can be of three types:
%           0 - periodic
%           1 - Dirichlet
%           2 - Neumann

%  NOTE! In case of periodic or Dirichlet boundary conditions
%  rho values are given only in the inner region of a grid
%  i.e. exluding the boundaries (it's natural for periodic 
%  boundary conditions because there are no boundaries at all).
%  In case of Neumann boundary conditions one must specify 
%  the charge density values also on the boundaries.

%  bcl, bcr, bcu, bcd -
%    vectors of ngy(bcl, bcr) or ngx(bcu,bcd) length;
%    represent the values of potential (in case of Dirichlet
%    boundary conditions), [V], or the values of first spatial
%    derivative in the direction perpendicular to the boundary
%    if the boundary conditions are of Neumann type, [V*m^-1];
%
%    if the bc are periodical bcl, bcr, bcu, bcd can be any 
%    number or array of numbers
%
%
%
%

%                               bcu
%       ------------------------------------------------------
%       |  #  |  #  |  #  |  #  |  #  |  #  |  #  |  #  |  #  |
%       ------------------------------------------------------
%       |  #  |     |     |     |     |     |     |     |  #  |
%       ------------------------------------------------------
%       |  #  |     |     |     |     |i,j-1|     |     |  #  |
%       ------------------------------------------------------
%       |  #  |     |     |     |i-1,j| i,j |i+1,j|     |  #  |
%   bcl ------------------------------------------------------   bcr
%       |  #  |     |     |     |     |i,j+1|     |     |  #  |
%       ------------------------------------------------------
%       |  #  |     |     |     |     |     |     |     |  #  |
%       ------------------------------------------------------
%       |  #  |     |     |     |     |     |     |     |  #  |
%       ------------------------------------------------------
%       |  #  |  #  |  #  |  #  |  #  |  #  |  #  |  #  |  #  |
%       ------------------------------------------------------
%                                bcd


eps0 = 8.85e-12; % [Farad/m]

[ngy ngx] = size(rho);
rho = rho*(-1)/eps0;

% calculate temp - [ngy ngx] array which is used on the second step of
% algorithm; fourier image of potential u(m,n) is equal to fourier image
% of charge dencity rho(m,n) divided by temp(m,n);

% temp(m,n) = (-4)*( (sin(pi/2*m/(M+1))/dx)^2 + (sin(pi/2*n/(N+1))/dy)^2 ),
% m = 1..M, n = 1..N  -  in case of Dirichlet bc

% temp(m,n) = (-4)*( (sin(pi/2*m/(M-1))/dx)^2 + (sin(pi/2*n/(N-1))/dy)^2 ),
% m = 0..M-1, n = 1..N-1  -  in case of Neumann bc

% temp(m,n) = (-4)*( (sin(pi*m/M)/dx)^2 + (sin(pi*n/N)/dy)^2 ),
% m = 0..M-1, n = 0..N-1  - in case of periodic bc

if bcx_type == 1
    rho = rho(:, 2:end-1);
    ngx = ngx - 2;
end
if bcy_type == 1
    rho = rho(2:end-1, :);
    ngy = ngy - 2;
end

switch bcx_type
    case 0
        hor_temp = sin(pi*(0:ngx-1)/ngx)/dx;
        hor_temp = hor_temp.*hor_temp;
        hor_temp = ones(ngy,1)*hor_temp;
        
    case 1
        hor_temp = sin(pi*(1:ngx)/2/(ngx+1))/dx;
        hor_temp = hor_temp.*hor_temp;
        hor_temp = ones(ngy,1)*hor_temp;  
        
    case 2
        hor_temp = sin(pi*(0:ngx-1)/2/(ngx-1))/dx;
        hor_temp = hor_temp.*hor_temp;
        hor_temp = ones(ngy,1)*hor_temp;          
end

switch bcy_type
    case 0
        vert_temp = sin(pi*(0:ngy-1)'/ngy)/dy;
        vert_temp = vert_temp.*vert_temp;
        vert_temp = vert_temp*ones(1,ngx);
        
    case 1
        vert_temp = sin(pi*(1:ngy)'/2/(ngy+1))/dy;
        vert_temp = vert_temp.*vert_temp;
        vert_temp = vert_temp*ones(1,ngx);    
        
    case 2
        vert_temp = sin(pi*(0:ngy-1)'/2/(ngy-1))/dy;
        vert_temp = vert_temp.*vert_temp;
        vert_temp = vert_temp*ones(1,ngx);  
end
    
temp = (-4)*(hor_temp + vert_temp);
if temp(1,1) == 0
    temp(1,1) = 1;
end

clear hor_temp vert_temp 

% adjust bcl, bcr, bcu, bcd to the purposes of the code
[t1 t2] = size(bcl);
if t1 == 1
    bcl = bcl';
end
[t1 t2] = size(bcr);
if t1 == 1
    bcr = bcr';
end
[t1 t2] = size(bcu);
if t2 == 1
    bcu = bcu';
end
[t1 t2] = size(bcd);
if t2 == 1
    bcd = bcd';
end

clear t1 t2

switch bcx_type
    case 1
        rho(:,1)  = rho(:,1)  - bcl/dx/dx;
        rho(:,ngx) = rho(:,ngx) - bcr/dx/dx;    
        
    case 2
        rho(:,1) = rho(:,1) + 2*bcl/dx;
        rho(:,ngx) = rho(:,ngx) - 2*bcr/dx;
end

switch bcy_type
    case 1
        rho(1,:)  = rho(1,:) - bcu/dy/dy;
        rho(ngy,:) = rho(ngy,:) - bcd/dy/dy;
        
    case 2
        rho(1,:)  = rho(1,:) + 2*bcu/dy;
        rho(ngy,:) = rho(ngy,:) - 2*bcd/dy;        
end

% main routine

% (fourier, sine, cosine) transform of each row of fi
fi = rho;

switch bcx_type
    case 0
        for i = 1:ngy
            fi(i,:) = fft(fi(i,:));
        end
        
    case 1
        for i = 1:ngy
            fi(i,:) = m_dst(fi(i,:),0);
        end
    
    case 2
        for i = 1:ngy
            fi(i,:) = m_dct(fi(i,:),0);
        end      
end

% (fourier, sine, cosine) transform of each column of fi
switch bcy_type
    case 0
        for i = 1:ngx
            fi(:,i) = fft(fi(:,i));
        end
        
    case 1
        for i = 1:ngx
            fi(:,i) = m_dst(fi(:,i),0);
        end
    
    case 2
        for i = 1:ngx
            fi(:,i) = m_dct(fi(:,i),0);
        end
end

fi = fi./temp;
if (bcx_type ~= 1)&(bcy_type ~= 1)
    fi(1,1) = 0;
end

% inverse transforms of rows and columns
switch bcx_type
    case 0
        for i = 1:ngy
            fi(i,:) = ifft(fi(i,:));
        end     
    case 1
        for i = 1:ngy
            fi(i,:) = m_dst(fi(i,:),1);
        end
    
    case 2
        for i = 1:ngy
            fi(i,:) = m_dct(fi(i,:),1);
        end
end

switch bcy_type
    case 0
        for i = 1:ngx
            fi(:,i) = ifft(fi(:,i));
        end
        
    case 1
        for i = 1:ngx
            fi(:,i) = m_dst(fi(:,i),1);
        end     
        
    case 2
        for i = 1:ngx
            fi(:,i) = m_dct(fi(:,i),1);
        end
end

fi = real(fi);

if bcx_type == 1
    fi = cat(2,bcl,fi,bcr);
    ngx = ngx + 2;
end

if bcy_type == 1
    fi = cat(1,bcu,fi,bcd);
    ngy  = ngy + 2;
end

A = zeros(ngy);
for i = 1:ngy-1
    A(i,i+1) = -1;
    A(i+1,i) = 1;
end
A(1,1) = 1;
A(end,end) = -1;
ey = A*fi/2/dy;

if bcy_type == 0
    ey(1,:) = (fi(end,:) - fi(2,:))/2/dy;
    ey(end,:) = (fi(end-1,:) - fi(1,:))/2/dy;
else
    ey(1,:) = ey(1,:)*2;
    ey(end,:) = ey(end,:)*2;
end

A = zeros(ngx);
for i = 1:ngx-1
    A(i,i+1) = -1;
    A(i+1,i) = 1;
end
A(1,1) = 1;
A(end,end) = -1;
ex = (A*fi'/2/dx)';

if bcx_type == 0
    ex(:,1) = (fi(:,end) - fi(:,2))/2/dx;
    ex(:,end) = (fi(:,end-1) - fi(:,1))/2/dx;
else
    ex(:,1) = ex(:,1)*2;
    ex(:,end) = ex(:,end)*2;
end

