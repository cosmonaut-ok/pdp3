function rho_back = load_rho_back(rho_back_struct, geometry, bc)

if strcmp(bc.x_type,'periodic')
    x_vect = geometry.dx/2:geometry.dx:geometry.x_size;
else
    x_vect = 0:geometry.dx:geometry.x_size;
end

if strcmp(bc.y_type,'periodic')
    y_vect = geometry.dy/2:geometry.dy:geometry.y_size;
else
    y_vect = 0:geometry.dy:geometry.y_size;
end

fhandle = rho_back_struct.profile_fun.handle
% rho_back = fhandle(x_vect,y_vect,geometry,rho_back_struct.profile_fun.param)*rho_back_struct.charge;
rho_back = linear_fun(x_vect,y_vect,geometry,rho_back_struct.profile_fun.param)*rho_back_struct.charge;


if ~(strcmp(bc.x_type,'periodic'))
    rho_back(:,1) = rho_back(:,1)*0.5;
    rho_back(:,end) = rho_back(:,end)*0.5;
end

if ~(strcmp(bc.y_type,'periodic'))
    rho_back(1,:) = rho_back(1,:)*0.5;
    rho_back(end,:) = rho_back(end,:)*0.5;
end
