function geometry = temp5(geometry, bc)

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