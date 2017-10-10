function [P p_internal_prop] = time_advance_2(P,p_internal_prop,dt,geometry,bc,ex,ey)

x_size = geometry.x_size;
y_size = geometry.y_size;
ngx = geometry.ngx;
ngy = geometry.ngy;
dx = geometry.dx;
dy = geometry.dy;
x_wall_type = bc.x_type;
y_wall_type = bc.y_type;

eq = 1.6e-19;
em = 9.1e-31;
eps0 = 8.85e-12;

c1 = eq/em*dt;

sq = dx*dy;

n_sort = length(p_internal_prop);
%must be rewritten!!!
% n_sort = 1;

for k = 1:n_sort
%     N = p_struct(k).N;
    charge = p_internal_prop(k).charge;
    mass = p_internal_prop(k).mass;

    x_interact_type = p_internal_prop(k).x_interact_type;
    y_interact_type = p_internal_prop(k).y_interact_type;

    c2 = c1*charge/mass;

    valid_part = find(P(k).type);

    x = P(k).x(valid_part);
    y = P(k).y(valid_part);
    vx = P(k).vx(valid_part);
    vy = P(k).vy(valid_part);

    if strcmp(x_wall_type, 'periodic')
        x_index = floor(x/dx + 0.5);
    else
        x_index = floor(x/dx) + 1; % cause the numeration of cells starts from 1
    end

    if strcmp(y_wall_type, 'periodic')
        y_index = floor(y/dy + 0.5);
    else
        y_index = floor(y/dy) + 1;
    end

    x_i_next = x_index + 1;
    y_i_next = y_index + 1;

    if strcmp(x_wall_type, 'periodic')
        x1 = (x_i_next - 0.5)*dx - x;
        x2 = dx - x1;
    else
        x1 = x_index*dx - x;
        x2 = dx - x1;
    end
    if strcmp(y_wall_type, 'periodic')
        y1 = (y_i_next - 0.5)*dy - y;
        y2 = dy - y1;

    else
        y1 = y_index*dy - y;
        y2 = dy - y1;
    end

    s1 = x1.*y1/sq;
    s2 = x1.*y2/sq;
    s3 = x2.*y1/sq;
    s4 = x2.*y2/sq;

    if strcmp(x_wall_type, 'periodic')

        zero_ind_part = find((x_index==0));
        x_index(zero_ind_part) = ngx;

        zero_ind_part = find((x_index==ngx));
        x_i_next(zero_ind_part) = 1;

    end

    if strcmp(y_wall_type, 'periodic')
        zero_ind_part = find((y_index==0));
        y_index(zero_ind_part) = ngy;

        zero_ind_part = find((y_index==ngy));
        y_i_next(zero_ind_part) = 1;
    end

    index_vect_1 = (x_index-1)*ngy + y_index;
    index_vect_2 = (x_index-1)*ngy + y_i_next;
    index_vect_3 = (x_i_next-1)*ngy + y_index;
    index_vect_4 = (x_i_next-1)*ngy + y_i_next;

    Fx =  (ex(index_vect_1).*s1 + ex(index_vect_2).*s2 + ...
        ex(index_vect_3).*s3 + ex(index_vect_4).*s4);

    Fy =  (ey(index_vect_1).*s1 + ey(index_vect_2).*s2 + ...
        ey(index_vect_3).*s3 + ey(index_vect_4).*s4);

    vx = vx + Fx*c2;
    vy = vy + Fy*c2;


    x = x + dt*vx;
    y = y + dt*vy;

    kill = 0;
    reflect = 0;

    x_less_0 = find(x<0);
    switch x_interact_type
        case 'cycling'
            x(x_less_0) = x(x_less_0) + x_size;
        case 'absorption'
            P(k).type(valid_part(x_less_0)) = 0;
%             valid_part(x_less_0) = [];
        case 'reflection'

            x(x_less_0) = -x(x_less_0);
            vx(x_less_0) = -vx(x_less_0);
        case 'antimirror'
            x(x_less_0) = -x(x_less_0);
            vx(x_less_zero) = -vx(x_less_0);
            reflect = 1;
    end
    x_more_x_size = find(x >= x_size);
    
    switch x_interact_type
        case 'cycling'
            x(x_more_x_size) = x(x_more_x_size) - x_size;
        case 'absorption'
            P(k).type(valid_part(x_more_x_size)) = 0;
%             valide_part(x_more_x_size) = [];
        case 'reflection'
            x(x_more_x_size) = 2*x_size - x(x_more_x_size);
            equal = find(x == x_size);
            x(equal) = x(equal) - 1e-10;
            vx(x_more_x_size) = -vx(x_more_x_size);
        case 'antimirror'
            x(x_more_x_size) = 2*x_size - x(x_more_x_size);
            equal = find(x == x_size);
            x(equal) = x(equal) - 1e-10;
            vx(x_more_x_size) = -vx(x_more_x_size);
            reflect = 1;
    end

    y_less_0 = find(y < 0);
    switch y_interact_type
        case 'cycling'
            y(y_less_0) = y(y_less_0) + y_size;
        case 'absorption'
            P(k).type(valid_part(y_less_0)) = 0;
%             valid_part(y_less_0) = [];
        case 'reflection'
            y(y_less_0) = -y(y_less_0);
            vy(y_less_0) = -vy(y_less_0);
        case 'antimirror'
            y(y_less_0) = -y(y_less_0);
            vy(y_less_0) = -vy(y_less_0);
            reflect = 1*(1 - reflect);
    end

    y_more_y_size = find(y >= y_size);
    switch y_interact_type
        case 'cycling'
            y(y_more_y_size) = y(y_more_y_size) - y_size;
        case 'absorption'
            P(k).type(valid_part(y_more_y_size)) = 0;
%             valid_part(y_more_y_size) = [];
        case 'reflection'
            y(y_more_y_size) = 2*y_size - y(y_more_y_size);
            equal = find(y == y_size);
            y(equal) = y(equal) - 1e-10;

            vy(y_more_y_size) = -vy(y_more_y_size);
        case 'antimirror'
            y(y_more_y_size) = 2*y_size - y(y_more_y_size);
            equal = find(y == y_size);
            y(equal) = y(equal) - 1e-10;
            vy(y_more_y_size) = -vy(y_more_y_size);
            reflect = 1*(1 - reflect);
    end

    %     if reflect

    
    P(k).x(valid_part) = x;
    P(k).y(valid_part) = y;
    P(k).vx(valid_part) = vx;
    P(k).vy(valid_part) = vy;

%     p_internal_prop(k).n_p = length(find(P(k).type ~= 0));

end


