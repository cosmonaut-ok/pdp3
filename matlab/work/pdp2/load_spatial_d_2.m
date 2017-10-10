function [resx resy] = load_particles_coords(geometry, concentration,randx,randy)

%the function assigns spatial coorditanes to particles;
%it uses 2xNbig values between 0 and 1 - randx, randy
%geometry and spatial concentration distribution - concentration

%10.03.2007
%The algorithm of particles coordinates assignment (in accordance with predefined 
%initial distribution function) involves the calculation of integral distribution 
%function of initial distribution function F(x): 
%D(x) = integrate(a,x,f(x))/integrate(a,b,f(x));
%the function load_spatial_d_2 performs this integraion numerically and
%makes the assignment of particles coordinates on the basis of two sets of
%random values randx and randy;
%in 2D geometry method of integral distribution function is decomposed on
%two subsequent 1D methods of IDF - firstly the dimensionality of initial 
%distribution funciton F(x,y) is redused by its integration along one
%dimension (here along y dimension): F(x,y) -> F'(x);
%than x coordinates of particles are calculated - x(i);
%after that the set of 1D functions F(x(i),y) is considered and y
%coordinates are calculated

N = length(randx);

x_size = geometry.x_size;
y_size = geometry.y_size;

%to perform a numerical integration the simulated volume is divided on
%small 2d bins
n_bin_x = 2048;
n_bin_y = 2048;

dx = x_size/n_bin_x;
dy = y_size/n_bin_y;

x = 0:dx:(x_size - dx);
y = 0:dy:(y_size - dy);

%the initial distribution function F(x,y) is calculated on the grid
F = concentration.handle(x,y,geometry,concentration.param);

%1D integral distribution function IDF(x) is calculated by
%integration(summation) along y dimension

IDF_1D = sum(F);

for k = 1:n_bin_x-1
    IDF_1D(k+1) = IDF_1D(k+1) + IDF_1D(k); 
end

if (IDF_1D(end) > 0)
  IDF_1D = IDF_1D/IDF_1D(end);
end
%-----------------------------------

%2D integral distribution function IDF(x,y) is calculated
IDF = F;
for k = 1:n_bin_y-1
    IDF(k+1,:) = IDF(k+1,:) + IDF(k,:);
end

for k = 1:n_bin_x
    if IDF(end,k) > 0
        IDF(:,k) = IDF(:,k)/IDF(end,k);
    end
end
%-----------------------------------

resx = zeros(1,N);
resy = zeros(1,N);

for k = 1:N
    resx(k) = b_search(randx(k), IDF_1D);
    resy(k) = b_search(randy(k), IDF(:,resx(k))); 
end
resx = resx/n_bin_x*x_size;
resy = resy/n_bin_y*y_size;
