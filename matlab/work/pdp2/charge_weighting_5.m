function rho = charge_weighting_5(p_internal_prop,geometry,bc,species_to_weight)

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
%        'periodic' - there are no borders
%        'dirichlet' - charge in the border cells is multiplied by 2
%        'neumann' - charge in the border cells is set to 0

%03.06.2006 
%different boundary conditions causes different division of 
%rectangular area into cells;
%the sipmlest bc type - periodic bc:
%and an area is divided as follows:   
%   
%      ____ ____ ____ ____  
%     |    |    |    |    |   
%     |____|____|____|____|
%     |    |    |    |    |   
%     |____|____|____|____|
%     |    |    |    |    |   
%     |____|____|____|____|    - all cells are located entirely in the simulated region;  
%            
%       
%such division method is situable also for mathematical reasons:
%periodic bc needs discrete fourier transform to be used - it means
%that the whole 2D space need to be filled with such rectangular blocks
% 
%another case is dirichlet and neumann types of bc:
%the pattern of simulated area division is another:
%   ____ ____ ____ ____ ____
%  |####|####|####|####|####|
%  |##__|____|____|____|__##|
%  |##  |    |    |    |  ##| 
%  |##__|____|____|____|__##|
%  |##  |    |    |    |  ##| 
%  |##__|____|____|____|__##|
%  |##  |    |    |    |  ##| - boundary cells lay partly in simulated area
%  |####|####|####|####|####|
%
%  ngx and ngy are greater by 1 than ngx and ngy for periodic bc.
% Such division is nessesary from mathematical reason. Poisson equation with
% dirichlet and neuman bc is solved by dicrete sine and cosine transforms
% respectively. OK, there is previous division type is also possible
% 
%
%  The procedure of charge weighting is also different for periodic and
%  non-periodic cases
%   

%10.03.2007
%New conception of global X,Y,VX,VY variables has been implementing. It
%requires some parts of code to be rewrited.

%10.07.2007
%The function charge_weighting is very memory consuming. One need to reduce
%the number of used variables.
%

ngx = geometry.ngx;
ngy = geometry.ngy;
dy = geometry.dy;
dx = geometry.dx;
x_wall_type = bc.x_type;
y_wall_type = bc.y_type;

rho = zeros(ngy,ngx);

sq = dx*dy;
eq = 1.6e-19;

global X Y VX VY F

% n_sp = length(p_internal_prop);
for k = species_to_weight
     %find particles which take part in the simulation;
     %there are empty places in X,Y,VX,VY array that marked
     %in F array with zero value
     valid_part = find(F(k).free);
     n_p = length(valid_part);
     if n_p > 0
        x = X(k).coord(valid_part)/dx;
        y = Y(k).coord(valid_part)/dy;
        
        %charge density, which corresponds to one particle in one cell  -
        %Cl/m^-3
        charge_dens = p_internal_prop(k).charge*p_internal_prop(k).lambda*eq/sq;

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
        % s1, s2, s3, s4 - parts of particle (with comparison to one)
        % which belong tho the cells with indexes (y_index, x_index),
        % (y_index+1, x_index), (y_index, x_index+1), (y_index+1, x_index+1)
        % respectively; s1+s2+s3+s4 = 1;
        %
        %  -----------------------------------------------------
        %  ||  s2(x_index, y_i_next) || s4(x_i_next, y_i_next) ||
        %  -----------------------------------------------------
        %  ||  s1(x_index, y_index)  || s3(x_i_next, y_index)  ||
        %  -----------------------------------------------------

        x = x*charge_dens;
        clear valid_part
        %-------------------------    
        
%         rho_temp = zeros(ngy,ngx);
        for i = 1:n_p
            % + s1;
            index_1 = (x_index(i)-1)*ngy + y_index(i);
            ff = find(index_1 < 0);
            if ~isempty(ff)
                index_1(ff)
                i
            end
            
            rho(index_1) = rho(index_1) + (charge_dens - x(i))*(1-y(i));
            
            % + s2;
            index_2 = (x_index(i)-1)*ngy + y_i_next(i);
            rho(index_2) = rho(index_2) + (charge_dens - x(i))*y(i);
            
            % + s3;
            index_3 = (x_i_next(i)-1)*ngy + y_index(i);
            rho(index_3) = rho(index_3) + x(i)*(1 - y(i));
            
            % + s4;
            index_4 = (x_i_next(i)-1)*ngy + y_i_next(i);
            rho(index_4) = rho(index_4) + x(i)*y(i);
        end
        
     end

end