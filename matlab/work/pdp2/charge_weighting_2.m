function rho = charge_weighting_2(P,p_struct,geometry,bc)

% the function charge_weighting uses an array of partile's
% charges and positions to determine charge density on the grid;

% P - array of records with fields .x, .y, .vx, vy, .mass, .charge
%     length(P) = Np;
%     x, y - [m], vx, vy - [m/s], mass - [kg], charge - [Cl];
%     x, y belongs to [0,dx*ngx], [0,dy*ngy] in case of 'cyclic'
%          and to [0,dx*(ngx-1)], [0,dy*(ngy-1)] else
%     each element of the array corresponds to one particle

% ngx, ngy - number of grid cells along x and y direction;
% dx, dy - geometrical sizes of cell; [m]

% x_type, y_type - additional information about system geometry
% type =
%        'cyclic' - there are no borders
%        'mirror' - charge in the border cells is multiplied by 2
%        'antimirror' - charge in the border cells is set to 0

%03.06.2006 the question is - why x,y are in [0,dx*(ngx-1)], [0,dy*(ngy-1)]
%more logical that x,y are in [dx/2,dx*ngx-dx/2)], [dy/2,dy*ngy-dy/2]

ngx = geometry.ngx;
ngy = geometry.ngy;
dy = geometry.dy;
dx = geometry.dx;
x_wall_type = bc.x_part_wall_type;
y_wall_type = bc.y_part_wall_type;

rho = zeros(ngy,ngx);
sq = dx*dy;
eq = 1.6e-19;

n_sort = length(P);
for k = 1:n_sort
    N = p_struct(k).N;
    if N > 0
        valid_part = find(P(k).type);

        x = P(k).x(valid_part);
        y = P(k).y(valid_part);

        charge_dens = p_struct(k).charge*p_struct(k).lambda*eq/sq;

        % particle is located in 4 cells;
        % now determine indexes of cell which is
        % the closest to the origin

        if strcmp(x_wall_type, 'cyclic')
            x_index = floor(x/dx + 0.5);
        else
            x_index = floor(x/dx) + 1; % cause the numeration of cells starts from 1
        end

        if strcmp(y_wall_type, 'cyclic')
            y_index = floor(y/dy + 0.5);
        else
            y_index = floor(y/dy) + 1;
        end

        x_i_next = x_index + 1;
        y_i_next = y_index + 1;

        % s1, s2, s3, s4 - parts of particle (with comparison to one)
        % which belong tho the cells with indexes (y_index, x_index),
        % (y_index+1, x_index), (y_index, x_index+1), (y_index+1, x_index+1)
        % respectively; s1+s2+s3+s4 = 1;
        if strcmp(x_wall_type, 'cyclic')
            x1 = (x_i_next - 0.5)*dx - x;
            x2 = dx - x1;
        else
            x1 = x_index*dx - x;
            x2 = dx - x1;
        end
        if strcmp(y_wall_type, 'cyclic')
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

        if strcmp(x_wall_type, 'cyclic')

            zero_ind_part = find((x_index==0));
            x_index(zero_ind_part) = ngx;

            ngx_ind_part = find((x_index==ngx));
            x_i_next(ngx_ind_part) = 1;

        end

        if strcmp(y_wall_type, 'cyclic')
            zero_ind_part = find((y_index==0));
            y_index(zero_ind_part) = ngy;

            ngy_ind_part = find((y_index==ngy));
            y_i_next(ngy_ind_part) = 1;
        end

        index_vect_1 = (x_index-1)*ngy + y_index;
        index_vect_2 = (x_index-1)*ngy + y_i_next;
        index_vect_3 = (x_i_next-1)*ngy + y_index;
        index_vect_4 = (x_i_next-1)*ngy + y_i_next;

        for i = 1:N

            rho(index_vect_1(i)) = rho(index_vect_1(i)) + s1(i)*charge_dens;

            rho(index_vect_2(i)) = rho(index_vect_2(i)) + s2(i)*charge_dens;

            rho(index_vect_3(i)) = rho(index_vect_3(i)) + s3(i)*charge_dens;

            rho(index_vect_4(i)) = rho(index_vect_4(i)) + s4(i)*charge_dens;

        end

    end
end

% if strcmp(x_type, 'mirror')
%     rho(:,1) = 2*rho(:,1);
%     rho(:,end) = 2*rho(:,end);
% end
% if strcmp(y_type, 'mirror')
%     rho(1,:) = 2*rho(1,:);
%     rho(end,:) = 2*rho(end,:);
% end
%
% if strcmp(x_type, 'antimirror')
%     rho(:,1) = 0*rho(:,1);
%     rho(:,end) = 0*rho(:,end);
% end
% if strcmp(y_type, 'antimirror')
%     rho(1,:) = 0*rho(1,:);
%     rho(end,:) = 0*rho(end,:);
% end

