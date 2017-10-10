function lambda = get_lambda(geometry, concentration, Nbig);

x_size = geometry.x_size;
y_size = geometry.y_size;

n_bin_x = 2048;
n_bin_y = 2048;

dx = x_size/n_bin_x;
dy = y_size/n_bin_y;

x = dx/2:dx:(x_size - dx/2);
y = dy/2:dy:(y_size - dy/2);


z = concentration.handle(x,y,geometry,concentration.param);

lambda = sum(sum(z))*dx*dy/Nbig;

