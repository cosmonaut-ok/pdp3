function geometry = calc_grid_step(geometry, bc)

%07/10/2006
%grid step calculation
%Division of simulated volume on the cells depends on type of boundary
%conditions;
%if bc.type is periodic then dx = x_size/ngx,
%else dx = x_size(ngx-1).

if strcmp(bc.x_type, 'periodic')
    geometry.dx = geometry.x_size/geometry.ngx;
else
    geometry.dx = geometry.x_size/(geometry.ngx - 1);
end

if strcmp(bc.y_type, 'periodic')
    geometry.dy = geometry.y_size/geometry.ngy;
else
    geometry.dy = geometry.y_size/(geometry.ngy - 1);
end