function res = time_advance_4(p_internal_prop,time,geometry,bc,ex,ey)

%10.03.2007
%the function should be rewritten after the conception of global variables
%is implemented

%11.07.2007
%In order to make more efficient the memory usage one need to reduce the
%number of variables
%

%12.08.2007
%make ex ey fi global

global X Y VX VY F

x_size = geometry.x_size;
y_size = geometry.y_size;
ngx = geometry.ngx;
ngy = geometry.ngy;
dx = geometry.dx;
dy = geometry.dy;
x_wall_type = bc.x_type;
y_wall_type = bc.y_type;
dt = time.dt;

eq = 1.6e-19;
em = 9.1e-31;
eps0 = 8.85e-12;

c1 = eq/em*dt;

sq = dx*dy;

n_sp = length(p_internal_prop);

for k = 1:n_sp

    charge = p_internal_prop(k).charge;
    mass = p_internal_prop(k).mass;

    x_interact_type = p_internal_prop(k).x_interact_type;
    y_interact_type = p_internal_prop(k).y_interact_type;

    c2 = c1*charge/mass;

    valid_part = find(F(k).free);

    x = X(k).coord(valid_part)/dx;
    y = Y(k).coord(valid_part)/dy;
    vx = VX(k).velocity(valid_part);
    vy = VY(k).velocity(valid_part);
%     X(1).coord(242)
    %generally a particle is located in 4 cells;
    %x_index, y_index, x_i_next, y_i_next - indexes of cells, where a
    %particles is located

    if strcmp(x_wall_type, 'periodic')
        x = x + 0.5;
        x_index = floor(x);
        %x = x - floor(x);
        x = x - x_index;
    else
        x_index = floor(x) + 1; % cause the numeration of cells starts from 1
        x = x - x_index + 1;
    end

    if strcmp(y_wall_type, 'periodic')
        y = y + 0.5;
        y_index = floor(y);
        y = y - y_index;
    else
        y_index = floor(y) + 1;
        y = y - y_index + 1;
    end

    x_i_next = x_index + 1;
    y_i_next = y_index + 1;

    %one needs to take into account the periodicity of bc, if a
    %particle is on the simulation region's border
    if strcmp(x_wall_type, 'periodic')
        x_index(find((x_index==0))) = ngx;
        x_i_next(find((x_index==ngx))) = 1;
    end

    if strcmp(y_wall_type, 'periodic')
        y_index(find((y_index==0))) = ngy;
        y_i_next(find((y_index==ngy))) = 1;
    end

    ex = ex*c2;
    ey = ey*c2;

    for j = 1:length(valid_part)

        index_1 = (x_index(j) - 1)*ngy + y_index(j);
        index_2 = (x_index(j) - 1)*ngy + y_i_next(j);
        index_3 = (x_i_next(j) - 1)*ngy + y_index(j);
        index_4 = (x_i_next(j) - 1)*ngy + y_i_next(j);
        s1 = (1 - x(j))*(1 - y(j));
        s2 = (1 - x(j))*y(j);
        s3 = x(j)*(1 - y(j));
        s4 = x(j)*y(j);

        vx(j) = vx(j) + ex(index_1)*s1 + ex(index_2)*s2 +ex(index_3)*s3 + ex(index_4)*s4;
        vy(j) = vy(j) + ey(index_1)*s1 + ey(index_2)*s2 +ey(index_3)*s3 + ey(index_4)*s4;

    end
    ex = ex/c2;
    ey = ey/c2;

    x = X(k).coord(valid_part);
    y = Y(k).coord(valid_part);
    
    x = x + dt*vx;
    y = y + dt*vy;
    
    x_less_0 = find(x<0);
    switch x_interact_type
        case 'cycling'
            x(x_less_0) = x(x_less_0) + x_size;
        case 'absorption'
            F(k).free(valid_part(x_less_0)) = 0;
        case 'reflection'
            x(x_less_0) = -x(x_less_0);
            vx(x_less_0) = -vx(x_less_0);
    end

    x_more_x_size = find(x > x_size);
    switch x_interact_type
        case 'cycling'
            x(x_more_x_size) = x(x_more_x_size) - x_size;
        case 'absorption'
            F(k).free(valid_part(x_more_x_size)) = 0;
        case 'reflection'
            x(x_more_x_size) = 2*x_size - x(x_more_x_size);
            vx(x_more_x_size) = -vx(x_more_x_size);
    end

    y_less_0 = find(y < 0);
    switch y_interact_type
        case 'cycling'
            y(y_less_0) = y(y_less_0) + y_size;
        case 'absorption'
            F(k).free(valid_part(y_less_0)) = 0;
        case 'reflection'
            y(y_less_0) = -y(y_less_0);
            vy(y_less_0) = -vy(y_less_0);
    end

    y_more_y_size = find(y > y_size);
    switch y_interact_type
        case 'cycling'
            y(y_more_y_size) = y(y_more_y_size) - y_size;
        case 'absorption'
            F(k).free(valid_part(y_more_y_size)) = 0;
        case 'reflection'
            y(y_more_y_size) = 2*y_size - y(y_more_y_size);
            vy(y_more_y_size) = -vy(y_more_y_size);
    end

    X(k).coord(valid_part) = x;
    Y(k).coord(valid_part) = y;
    VX(k).velocity(valid_part) = vx;
    VY(k).velocity(valid_part) = vy;
end

res = 1;
return;
