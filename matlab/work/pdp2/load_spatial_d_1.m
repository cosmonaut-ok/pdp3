function [resx resy] = load_spatial_d_1(geometry, concentration,randx,randy)

%the function assigns spatial coorditanes to particles;
%it uses 2xNbig values between 0 and 1 - randx, randy
%geometry and spatial concentration distribution - concentration

%10.03.2007
%more specifically it performs

N = size(randx,2);

x_size = geometry.x_size;
y_size = geometry.y_size;

n_bin_x = 2048;
n_bin_y = 2048;

dx = x_size/n_bin_x;
dy = y_size/n_bin_y;

% x = 0:dx:(x_size - dx) + x0 + dx/2;
% y = 0:dy:(y_size - dy) + y0 + dy/2;
x = 0:dx:(x_size - dx);
y = 0:dy:(y_size - dy);

d = concentration.handle(x,y,geometry,concentration.param);

d1 = sum(d);

for k = 1:n_bin_x-1
    d1(k+1) = d1(k+1) + d1(k); 
end

if (d1(end) > 0)
  d1 = d1/d1(end);
end

for k = 1:n_bin_y-1
    d(k+1,:) = d(k+1,:) + d(k,:);
end

for k = 1:n_bin_x
    if d(end,k) > 0
        d(:,k) = d(:,k)/d(end,k);
    end
end

resx = zeros(1,N);
resy = zeros(1,N);

for k = 1:N

    resx(k) = b_search(randx(k), d1);
    resy(k) = b_search(randy(k), d(:,resx(k))); 
end
resx = resx/n_bin_x*x_size;
resy = resy/n_bin_y*y_size;
