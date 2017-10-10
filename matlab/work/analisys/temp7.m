function res = temp7

geometry.ngx = 1024;
geometry.ngy = 1024;
geometry.x_size = 0.6;
geometry.y_size = 0.6;

bc.x_type = 'dirichlet';
bc.y_type = 'periodic';

bc.left_value = 10000000000;
bc.right_value = 0;
bc.top_value = 0;
bc.bottom_value = 0;

geometry = calc_grid_step(geometry,bc);
rho = ones(geometry.ngy,geometry.ngx);

fi = ones(geometry.ngy, geometry.ngx);

tic;
fi = field_3(rho, geometry, bc);
toc;

surf(fi)

res = 1;
return;