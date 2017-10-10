function P = beam_injection_1(P,inject_param,geometry,time)

% The important feature of "PDP2, the plasma simulator" is the possibility of
% deep and detailed investigation of beam-plasma interaction.

% The beam_injection_1 function performs an injection of beam particles into
% the system. Beams could be injected from the left electrode at different 
% positions, widthes and velocities. They could be also modulated by
% density or by velocity.

% The variable inject_param is an array of structures which provides all
% necessary information for particle's injection in the system. The length of
% inject_param array is equal to the number of beam species in the system.
% 
% inject_param structures consist of such fields:
% .n_injected(k) - number of particles of k-th specie which are injected
%                  into the simulated volume per one time step of simulation;
% .vx(k), .vy(k) - initial velocity of injected particle of k-th specie; 
% .beam_width(k) - beam width in normalized (to the length of the system in
%                  the y direction) units;
% .beam_pos(k)   - the position of the center of a beam in normalized (to
%                  the length of the system in the y direction) units;
% .mod_type(k)   - type of modulation: 'density' - density modulation, 
%                  'velocity' - velocity modulation;
% .mod_depth(k)  - depth of modulation, i.e. the relative alteration of the
%                  beam amplitude;
% .mod_frqn(k)   - frequency of modulation;

%10.03.2007
%function is modificated because of global variables introduction in the
%program

global X,Y,VX,VY,F

beam_n_sp = length(inject_param);
total_n_sp = length(X);

for k = 1:beam_n_sp
    cur_sp = total_n_sp - beam_n_sp + 1;
    n_injected = inject_param(k).n_injected;
    vx = inject_param(k).vx;
    vy = inject_param(k).vy;
    beam_pos = inject_param(k).beam_pos;
    beam_width = inject_param(k).beam_width;
    y_size = geometry.y_size;
    mod_type = inject_param(k).mod_type;
    mod_depth = inject_param(k).mod_depth;
    mod_frqn = inject_param(k).mod_frqn;
    
    switch mod_type
        case 'density'
            n_injected = round(n_injected*(1 + mod_depth*sin(mod_frqn*time.t_current)));
        case 'velocity'
            vx = vx*(1 + mod_depth*sin(mod_frqn*time.t_current));
            vy = vy*(1 + mod_depth*sin(mod_frqn*time.t_current));
    end

    if vx*time.dt > y_size
     disp('Warning: the beam is too intensive');
     return;
    end
    
    free_cells = find(F(cur_sp).free == 0);
    
    if (length(free_cells) < n_injected)
        disp('Error: number of particles exceens its maximum');
        return;
    end
    
    free_cells = free_cells(1:n_injected);
    F(cur_sp).free(free_cells) = 1;
    X(cur_sp).coord(free_cells) = vx*time.dt*rand(1,n_injected);
    Y(cur_sp).coord(free_cells) = (beam_pos + beam_width*(rand(1,n_injected)-0.5))*y_size;
    VX(cur_sp).velocity(free_cells) = vx;
    VY(cur_sp).velocity(free_cells) = vy;
end