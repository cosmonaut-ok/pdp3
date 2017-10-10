% function [P p_table Np] = load_particles(N,type,TeV,x0,xend,y0,yend,lambda)
function P = load_particles_1(p_struct, geometry)

% each record in _struct corresponds to the one sort of particles;
% it has the follows fields:
% N - number of particles of the specified sort;
% Nmax - maximal number of particles of the specified sort;
% mass - mass of the elementary particle (in electron mass units);
% charge - charge of the elem. particle (in elect. charge units *(-1));
% lambda - linear parameter - number of elem. part. in big part [number/m];
% geometry - [x0 y0 xSize ySize] - determine the space region of particle
%            distribution;
% f_handle - determine the function of particle distribution;
% sp_distr_type - determine the type of particle distribution in space;
% TeV - temperature of specified sort of particles



% This function determines the initial distribution
% of big particles in x- and v-spaces.

% Input:
%    N - number of big particles (in present version
%        Ne = Ni = N)
%    type - type of start. The variable may take on
%           two values - 'quiet' and 'random'

% Output:
%    P - structure array which contains all information
%        about particles in the system; P has the following
%        fields:
%        .type - type of big particle(electron, ion, neutral,
%                  beam);
%        .x - x coordinate of particle;
%        .y - y coordinate of particle;
%        .vx - velocity of particle along x axis;
%        .vy - velocity of particle along y axis;
%        .charge - charge of particle in elementary charge units;
%        .mass - mass of particle in electron mass units;

qe = 1.6e-19;
me = 9.1e-31;
mp2me = 1.836e3;
kb = 1.38e-23;
eVconst = 1.16e4;

n_sort = length(p_struct);


two_pi = 2*pi;

    x_size = geometry.x_size;
    y_size = geometry.y_size;

for k = 1:n_sort

    n_p = p_struct(k).init_load_param.n_p;
    n_p_max = p_struct(k).init_load_param.n_p_max;
    load_type = p_struct(k).init_load_param.load_type;

    T = p_struct(k).init_load_param.velocity_fun.param(1)*eVconst;
    sp_fun = p_struct(k).init_load_param.spatial_fun;

    if n_p_max > 0
        zer_vect = zeros(1,n_p_max);
        P(k).x = zer_vect;
        P(k).y = zer_vect;
        P(k).vx = zer_vect;
        P(k).vy = zer_vect;
        P(k).type = zer_vect;
        if n_p > 0
            switch load_type

                case 'random'
                    
                    switch func2str(sp_fun.handle)
                        case 'linear_fun'
%                             [resx resy] = load_linear_d(geometry, sp_fun,rand(1,n_p),rand(1,n_p));
                            [resx resy] = load_linear_d(geometry, sp_fun,revers(1:n_p,3),revers(1:n_p,5));
                        case 'homogene_fun'
%                             resx = rand(1,n_p)*geometry.x_size;
%                             resy = rand(1,n_p)*geometry.y_size;
                            
                            resx = revers(1:n_p,3)*geometry.x_size;
                            resy = revers(1:n_p,5)*geometry.y_size;
                            
                        otherwise
                            [resx resy] = load_spatial_d_1(geometry, sp_fun,rand(1,n_p),rand(1,n_p));
                    end

                case 'uniform_in_space'
                    [resx resy] = load_spatial_d_1(geometry, sp_fun,revers(1:n_p,3),revers(1:n_p,5));

                case 'quiet'

                otherwise

            end
            
            charge = p_struct(k).internal_prop.charge;
            mass = p_struct(k).internal_prop.mass;
            vconst = (2*kb*T/me)^0.5*mass^(-0.5);
            
            Rv_vect = rand(1,n_p);
            Rtheta_vect = rand(1,n_p);
            v = vconst*(-log(Rv_vect)).^0.5;
            theta = two_pi*Rtheta_vect;
            vx = v.*cos(theta);
            vy = v.*sin(theta);
            
            P(k).x(1:n_p) = resx;
            P(k).y(1:n_p) = resy;
            P(k).vx(1:n_p) = vx;
            P(k).vy(1:n_p) = vy;
            P(k).type(1:n_p) = k*ones(1,n_p);
            end
        end
    end