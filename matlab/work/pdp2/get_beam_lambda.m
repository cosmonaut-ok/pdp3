function lambda = get_beam_lambda(inject_param, geometry, time);

%the function calculates linear density of beam particles

%density (or mean density in case of modulation) of the beam in m^-3; one
%of the injection parameters
density = inject_param.density;

%width of beam in normalized (to Ly) units;
beam_width = inject_param.beam_width;

%absolute velocity of the beam (m/s);

vx = inject_param.vx;
vy = inject_param.vy;

v = (vx^2 + vy^2)^0.5;

%number of injected particles per time step
n_injected = inject_param.n_injected;

%value time step (s)
dt = time.dt;

%width of the system (m)
y_size = geometry.y_size;

lambda = round(density*beam_width*y_size*v*dt/n_injected);