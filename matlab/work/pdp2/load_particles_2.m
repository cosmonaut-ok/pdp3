% function [P p_table Np] = load_particles(N,type,TeV,x0,xend,y0,yend,lambda)
function res = load_particles_2(p_struct, geometry)

% each record in p_struct corresponds to the one sort of particles;
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
%10.03.2007
% The nessesity to rewrite this function has appeared as the conception of
% the main variable P had changed. From now on the variable P is divided on
% 4 variables X,Y,VX,VY which are global througout the several function

global X Y VX VY F

qe = 1.6e-19;
me = 9.1e-31;
mp2me = 1.836e3;
kb = 1.38e-23;
eVconst = 1.16e4;

n_sp = length(p_struct);


two_pi = 2*pi;

    x_size = geometry.x_size;
    y_size = geometry.y_size;

for k = 1:n_sp

    n_p = p_struct(k).init_load_param.n_p;
    load_type = p_struct(k).init_load_param.load_type;

    T = p_struct(k).init_load_param.velocity_fun.param(1)*eVconst;
    sp_fun = p_struct(k).init_load_param.spatial_fun;
    if n_p > 0
        switch load_type
            case 'random'
%                 seed_x = rand(1, n_p);
%                 seed_y = rand(1, n_p);
                seed_x = revers(1:n_p,3);
                seed_y = revers(1:n_p,5);                
            case 'uniform_in_space'
                seed_x = revers(1:n_p,3);
                seed_y = revers(1:n_p,5);
        end
        
        [resx resy] = load_particles_coords(geometry, sp_fun, seed_x, seed_y);
            
        charge = p_struct(k).internal_prop.charge;
        mass = p_struct(k).internal_prop.mass;
        vconst = (2*kb*T/me)^0.5*mass^(-0.5);
            
        Rv_vect = rand(1,n_p);
        Rtheta_vect = rand(1,n_p);
        v = vconst*(-log(Rv_vect)).^0.5;
        theta = two_pi*Rtheta_vect;
        vx = v.*cos(theta);
        vy = v.*sin(theta);
            
        X(k).coord = resx;
        Y(k).coord = resy;
        VX(k).velocity = vx;
        VY(k).velocity = vy;
        F(k).free(1:n_p) = ones(1,n_p);
    end

end

res = 1;
return;