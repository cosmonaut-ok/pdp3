function rho = charge_weighting_6(p_internal_prop,geometry,bc,species_to_weight)

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
%Superscalar efficiency of MATLAB appears only for simple operations. Each
%particles iteration method could be more efficient of charge weighting;


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

for k = species_to_weight
     %find particles which take part in the simulation;
     %there are empty places in X,Y,VX,VY array that marked
     %in F array with zero value
     valid_part = find(F(k).free);
     n_p = length(valid_part);
     x = X(k).coord(valid_part)/dx;
     y = Y(k).coord(valid_part)/dx;     

     %charge density, which corresponds to one particle in one cell - Cl/m^-3
     charge_dens = p_internal_prop(k).charge*p_internal_prop(k).lambda*eq/sq;     
     
     for j = 1:n_p
     
        %generally a particle is located in 4 cells;
        %now determine indexes of cell which is
        %the closest to the origin

        if strcmp(x_wall_type, 'periodic')
            xj = x(j) + 0.5;
            x_index = floor(xj);
            x_i_next = x_index + 1;
            
            if x_index == 0
                x_index = ngx;
            end
            if x_index == ngx;
                x_i_next = 1;
            end
            
        else
            xj = x(j);
            x_index = floor(xj) + 1; % cause the numeration of cells starts from 1
            x_i_next = x_index + 1;
        end

        if strcmp(y_wall_type, 'periodic')
            yj = y(j) + 0.5;
            y_index = floor(yj);
            y_i_next = y_index + 1;
            
            if y_index == 0
                y_index = ngy;
            end
            if y_index == ngy;
                y_i_next = 1;
            end
            
        else
            yj = y(j);
            y_index = floor(yj) + 1; % cause the numeration of cells starts from 1
            y_i_next = y_index + 1;
        end        

% s1, s2, s3, s4 - parts of particle (with comparison to one)
        % which belong tho the cells with indexes (y_index, x_index),
        % (y_index+1, x_index), (y_index, x_index+1), (y_index+1, x_index+1)
        % respectively; s1+s2+s3+s4 = 1;        
                
        if strcmp(x_wall_type, 'periodic')
            s3 = xj - floor(xj);
            s4 = xj - floor(xj);
            s1 = 1 - s3;
            s2 = 1 - s3;
        else
            s3 = xj - floor(xj);
            s4 = xj - floor(xj);
            s1 = 1 - s3;
            s2 = 1 - s3;
        end
        
        if strcmp(y_wall_type, 'periodic')
            s2 = s2.*(yj - floor(yj));
            s4 = s4.*(yj - floor(yj));
            s1 = s1 - s2;
            s3 = s3 - s4;
        else
            s2 = s2.*(yj - floor(yj));
            s4 = s4.*(yj - floor(yj));
            s1 = s1 - s2;
            s3 = s3 - s4;
        end        
    
        rho(x_index, y_index) = rho(x_index, y_index) + s1;
        rho(x_index, y_index+1) = rho(x_index, y_index+1) + s2;
        rho(x_index+1, y_index) = rho(x_index+1, y_index) + s3;
        rho(x_index+1, y_index+1) = rho(x_index+1, y_index+1) + s4;
        
        rho = rho*charge_dens;
        
    end
end


