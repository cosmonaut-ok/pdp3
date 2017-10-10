function res = field_5(rho, geometry, bc)

%  10.01.2006
%  the function solves discrete 2D Poisson's equation;
%  rho - 2D matrix with charge density values, [Cl*m^-3]
%


%
%  [ngy ngx] = size(rho) - size of two-dimensional 
%                          rectangular grid(including borders)
%  geometry.dx, geometry.dy - length and width of elementary grid cell, [m]
%
%  bc.x_type, bc.y_type - type of boundary conditions on
%                         left, right(x_type) and
%                         upper, lower(y_type) boundaries
%  boundary conditions can be of three types: periodic, Dirichlet, Neumann

%  NOTE! In case of periodic or Dirichlet boundary conditions
%  rho values are given only in the inner region of a grid
%  i.e. exluding the boundaries (it's natural for periodic 
%  boundary conditions because there are no boundaries at all).
%  In case of Neumann boundary conditions one must specify 
%  the charge density values also on the boundaries.

%  bcl, bcr, bct, bcb -
%    vectors of ngy(bcl, bcr) or ngx(bct,bcb) length;
%    represent the values of potential (in case of Dirichlet
%    boundary conditions), [V], or the values of first spatial
%    derivative in the direction perpendicular to the boundary
%    if the boundary conditions are of Neumann type, [V*m^-1];
%
%    if the bc are periodical bcl, bcr, bct, bcb can be any 
%    number or array of numbers
%
%
%
%

%                               bct
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
%                                bcb


%20.07.2006
%The function didn't work in case when bc of all electrodes are Dirichlet
%or when on the one pair of electrodes are Dirichlet bc and on the other
%are Neumann bc

%26.08.2006
%small changes: from now on boundary conditions and geometry are put into
%the function as a whole structures; fields of bc are .x_type, .y_type,
%.top_value, .bottom_value, .left_value, .right_value, .left_vector, 
%.right_vector, .top_vector, .bottom_vector;
%geometry - .dx, .dy, .ngx, .ngy;

%10.03.2007
%The function is slightly simplified. Such parameters as bc.right_vector,
%bc.left_vector, bc.top_vector, bc.bottom_vector is not specified in the 
%input bc structure. These vectors are constructed in the function. 

%12.08.2007
%Introduce global variables ex, ey, fi

global fi

eps0 = 8.85e-12; % [Farad/m]

rho = rho*(-1)/eps0;

x_type = bc.x_type;
y_type = bc.y_type;

ngx = geometry.ngx;
ngy = geometry.ngy;
dx = geometry.dx;
dy = geometry.dy;

bcl = bc.left_value*ones(geometry.ngy,1);
bcr = bc.right_value*ones(geometry.ngy,1);
bct = bc.top_value*ones(1,geometry.ngx);
bcb = bc.bottom_value*ones(1,geometry.ngx);

%analisys of the boundary conditions correctness
%more detailed: if the bc are of the same type on the all electrodes (both
%on the left/right and on the top/bottom then the certain elements of the
%bc vectors must coincide, e.g. the first value of the bcl vector and the 
%first value of the bct vector must be equal as soon as they represent the 
%value of the potential (or normal component of the electric field) in the 
%top left corner of the grid

if (strcmp(x_type,'dirichlet')&strcmp(y_type,'dirichlet'))|(strcmp(x_type,'neumann')&strcmp(y_type,'neumann'))
    if (bct(1) ~= bcl(1))
        disp('Warning: the fist value of bct vector and the first value of bcl vector should coincide!');
    end
    up_left_value = (bct(1) + bcl(1))/2;
    
    if (bct(end) ~= bcr(1))
        disp('Warning: the last value of bct vector and the first value of bcr vector should coincide!');
    end 
    up_right_value = (bct(end) + bcr(1))/2;
    
    if (bcb(end) ~= bcr(end))
        disp('Warning: the last value of bcb vector and the last value of bcr vector should coincide!');
    end 
    down_right_value = (bcb(end) + bcr(end))/2;   
    
    if (bcb(1) ~= bcl(end))
        disp('Warning: the first value of bcb vector and the last value of bcl vector should coincide!');
    end 
    down_left_value = (bcb(1) + bcl(end))/2;  
end

%in case of Dirichlet bc (the values of potential are given on the boundary)
%we don't have a need to know the values of charge density on the boundary
%to solve the Poisson equation

%the length of vectors which represent the values on the boundary depends on
%the types of boundary conditions;

if strcmp(x_type, 'dirichlet')
    rho = rho(:, 2:end-1);
    ngx = ngx - 2;
    bct = bct(2:end-1);
    bcb = bcb(2:end-1);
end

if strcmp(y_type, 'dirichlet')
    rho = rho(2:end-1, :);
    ngy = ngy - 2;    
    bcl = bcl(2:end-1);
    bcr = bcr(2:end-1);
end

% calculate temp - [ngy ngx] array which is used on the second step of
% algorithm; fourier image of potential u(m,n) is equal to fourier image
% of charge dencity rho(m,n) divided by temp(m,n);

% temp(m,n) = (-4)*( (sin(pi/2*m/(M+1))/dx)^2 + (sin(pi/2*n/(N+1))/dy)^2 ),
% m = 1..M, n = 1..N  -  in case of Dirichlet bc

% temp(m,n) = (-4)*( (sin(pi/2*m/(M-1))/dx)^2 + (sin(pi/2*n/(N-1))/dy)^2 ),
% m = 0..M-1, n = 0..N-1  -  in case of Neumann bc

% temp(m,n) = (-4)*( (sin(pi*m/M)/dx)^2 + (sin(pi*n/N)/dy)^2 ),
% m = 0..M-1, n = 0..N-1  - in case of periodic bc

switch x_type
    case 'periodic'
        hor_temp = sin(pi*(0:ngx-1)/ngx)/dx;
        hor_temp = hor_temp.*hor_temp;
        hor_temp = ones(ngy,1)*hor_temp;
        
    case 'dirichlet'
        hor_temp = sin(pi*(1:ngx)/2/(ngx+1))/dx;
        hor_temp = hor_temp.*hor_temp;
        hor_temp = ones(ngy,1)*hor_temp;  
        
    case 'neumann'
        hor_temp = sin(pi*(0:ngx-1)/2/(ngx-1))/dx;
        hor_temp = hor_temp.*hor_temp;
        hor_temp = ones(ngy,1)*hor_temp;          
end

switch y_type
    case 'periodic'
        vert_temp = sin(pi*(0:ngy-1)'/ngy)/dy;
        vert_temp = vert_temp.*vert_temp;
        vert_temp = vert_temp*ones(1,ngx);
        
    case 'dirichlet'
        vert_temp = sin(pi*(1:ngy)'/2/(ngy+1))/dy;
        vert_temp = vert_temp.*vert_temp;
        vert_temp = vert_temp*ones(1,ngx);    
        
    case 'neumann'
        vert_temp = sin(pi*(0:ngy-1)'/2/(ngy-1))/dy;
        vert_temp = vert_temp.*vert_temp;
        vert_temp = vert_temp*ones(1,ngx);  
end
    
temp = (-4)*(hor_temp + vert_temp);

% if (~strcmp(x_type, 'dirichlet'))&(~strcmp(y_type, 'dirichlet')) then the
% temp(1,1) is equal to 0; there will be a stage of the algorithm when we 
% divide a Fourier image of the rho by temp; to avoid a division on 0 we
% assign 1 to temp(1,1) 
if (~strcmp(x_type, 'dirichlet'))&(~strcmp(y_type, 'dirichlet'))
    temp(1,1) = 1;
end

clear hor_temp vert_temp 

% adjust bcl, bcr, bct, bcb to the purposes of the code
[t1 t2] = size(bcl);
if t1 == 1
    bcl = bcl';
end
[t1 t2] = size(bcr);
if t1 == 1
    bcr = bcr';
end
[t1 t2] = size(bct);
if t2 == 1
    bct = bct';
end
[t1 t2] = size(bcb);
if t2 == 1
    bcb = bcb';
end

clear t1 t2

switch x_type
    case 'dirichlet'
        rho(:,1)  = rho(:,1)  - bcl/dx/dx;
        rho(:,ngx) = rho(:,ngx) - bcr/dx/dx;    
        
    case 'neumann'
        rho(:,1) = rho(:,1) + 2*bcl/dx;
        rho(:,ngx) = rho(:,ngx) - 2*bcr/dx;
end

switch y_type
    case 'dirichlet'
        rho(1,:)  = rho(1,:) - bct/dy/dy;
        rho(ngy,:) = rho(ngy,:) - bcb/dy/dy;
        
    case 'neumann'
        rho(1,:)  = rho(1,:) + 2*bct/dy;
        rho(ngy,:) = rho(ngy,:) - 2*bcb/dy;        
end

% main routine

% (fourier, sine, cosine) transform of each row of fi
fi = rho;

switch x_type
    case 'periodic'
        for i = 1:ngy
            fi(i,:) = fft(fi(i,:));
        end
        
    case 'dirichlet'
        for i = 1:ngy
            fi(i,:) = m_dst(fi(i,:),0);
        end
    
    case 'neumann'    
        for i = 1:ngy
            fi(i,:) = m_dct(fi(i,:),0);
        end  
end

% (fourier, sine, cosine) transform of each column of fi
switch y_type
    case 'periodic'
        for i = 1:ngx
            fi(:,i) = fft(fi(:,i));
        end
        
    case 'dirichlet'
        for i = 1:ngx
            fi(:,i) = m_dst(fi(:,i),0);
        end
    
    case 'neumann'
        for i = 1:ngx
            fi(:,i) = m_dct(fi(:,i),0);
        end
end

fi = fi./temp;

%remember that element (1,1) of the result of division a Fourier image of
%rho by temp must be equal to 0 if 
%(~strcmp(x_type,'dirichlet'))&(~strcmp(y_type, 'dirichlet'))

if (~strcmp(x_type, 'dirichlet'))&(~strcmp(y_type, 'dirichlet'))
    fi(1,1) = 0;
end

% inverse transforms of rows and columns
switch x_type
    case 'periodic'   
        for i = 1:ngy
            fi(i,:) = ifft(fi(i,:));
        end  
        
    case 'dirichlet'
        for i = 1:ngy
            fi(i,:) = m_dst(fi(i,:),1);
        end
    
    case 'neumann'
        for i = 1:ngy
            fi(i,:) = m_dct(fi(i,:),1);
        end
end

switch y_type
    case 'periodic'
         for i = 1:ngx
            fi(:,i) = ifft(fi(:,i));
        end
        
    case 'dirichlet' 
        for i = 1:ngx
            fi(:,i) = m_dst(fi(:,i),1);
        end 
        
    case 'neumann'
        for i = 1:ngx
            fi(:,i) = m_dct(fi(:,i),1);
        end
end

fi = real(fi);

%addition of boundary values to the array with values of potential in case
%of dirichlet boundary conditions

if (strcmp(x_type,'dirichlet')&strcmp(y_type,'dirichlet'))
    bct = cat(2,up_left_value,bct,up_right_value);
    bcb = cat(2,down_left_value,bcb,down_right_value);
end

if strcmp(x_type, 'dirichlet')
        fi = cat(2,bcl,fi,bcr);
        ngx = ngx + 2;
end

if strcmp(y_type, 'dirichlet')
    fi = cat(1,bct,fi,bcb);
    ngy  = ngy + 2;
end

res = 1;
