function varargout = pdp2_4(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pdp2_4_OpeningFcn, ...
                   'gui_OutputFcn',  @pdp2_4_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before pdp2_2 is made visible.
function pdp2_4_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% Choose default command line output for pdp2_2
handles.output = hObject;

const.em = 9.1e-31;
const.eq = 1.6e-19;
const.eps0 = 8.85e-12;
const.kb = 1.38e-23;

handles.const = const;

load('init_param.mat', 'param');
handles.param = param;

handles.visual = initialize_visual(handles);


% handles.CurFig = get(0, 'CurrentFigure');
handles.state = 'Stopped';
handles.logfile = '';
set(handles.text1,'String','Stopped', 'ForegroundColor', [0 0 0]);
set(handles.Pause_MenuItem, 'Enable', 'off');
set(handles.Resume_MenuItem, 'Enable', 'off');
set(handles.Stop_MenuItem, 'Enable', 'off');
set(handles.save_sys_state_item, 'Enable', 'off');

% Update handles structure
guidata(hObject, handles);

visualization_field_update(hObject, guidata(hObject));

% UIWAIT makes pdp2_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = pdp2_4_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function File_MenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
%------------------FILE MENU COMMANDS---------------------------------

function results_saving_item_Callback(hObject, eventdata, handles)
beam_struct = handles.param.beam_struct;
p_struct = handles.param.p_struct;

if (length(beam_struct) > 0)
    
    for i = 1:length(p_struct)
        p_internal_prop(i) = p_struct(i).internal_prop;
    end
    for i = 1:length(beam_struct)
        beam_internal_prop(i) = beam_struct(i).internal_prop;
    end
    name_list = {p_internal_prop.name, beam_internal_prop.name};
else
    for i = 1:length(p_struct)
        p_internal_prop(i) = p_struct(i).internal_prop;
    end
    name_list = {p_internal_prop.name};  
end

handles.param.results_saving_struct = results_saving_param(name_list, handles.param.results_saving_struct);
guidata(hObject, handles)

% --------------------------------------------------------------------
function save_sim_param_item_Callback(hObject, eventdata, handles)

[file, path] = uiputfile('*.mat', 'Save simulation parameters as');
if isstr(path)
    param = handles.param;
    save([path file], 'param')
end

% --------------------------------------------------------------------
function load_sim_param_item_Callback(hObject, eventdata, handles)

[file, path] = uigetfile('*.mat', 'Load simulation parameters');
if isstr(path)
    load([path file], 'param')
    handles.param = param;
    guidata(hObject, handles);
end

set(handles.Run_MenuItem, 'Enable', 'on');
set(handles.Pause_MenuItem, 'Enable', 'off');
set(handles.Resume_MenuItem, 'Enable', 'off');
set(handles.Stop_MenuItem, 'Enable', 'off');
set(handles.save_sys_state_item, 'Enable', 'off');
set(handles.text1, 'String', 'Stopped', 'ForegroundColor', [0 0 0]);
handles.state = 'Stopped';

guidata(hObject, handles);

visualization_field_update(hObject, guidata(hObject));

% --------------------------------------------------------------------
function save_sys_state_item_Callback(hObject, eventdata, handles)
global X Y VX VY F

[file, path] = uiputfile('*.mat', 'Save current system state as');
if isstr(path)
    param = handles.param;
    logfile = handles.logfile;
    save([path file],'X','Y','VX','VY','F','param','logfile')
end

% --------------------------------------------------------------------
function load_sys_state_item_Callback(hObject, eventdata, handles)
global X Y VX VY F

[file, path] = uigetfile('*.mat', 'Load simulation parameters');
if isstr(path)
    load([path file],'X','Y','VX','VY','F','param','logfile')
    handles.param = param;
    handles.logfile = logfile;
    guidata(hObject, handles);
end

set(handles.Run_MenuItem, 'Enable', 'off');
set(handles.Pause_MenuItem, 'Enable', 'off');
set(handles.Resume_MenuItem, 'Enable', 'on');
set(handles.Stop_MenuItem, 'Enable', 'on');
set(handles.save_sys_state_item, 'Enable', 'on');
set(handles.text1, 'String', 'Paused', 'ForegroundColor', [0 0 1]);
handles.state = 'Paused';

visualization_field_update(hObject, guidata(hObject));

guidata(hObject, handles);

% --------------------------------------------------------------------
function Exit_MenuItem_Callback(hObject, eventdata, handles)

clear X Y VX VY F ex ey fi

delete(handles.figure1)

%-----------END OF FILE MENU COMMANDS---------------------------------
%---------------------------------------------------------------------

function Operation_MenuItem_Callback(hObject, eventdata, handles)
pause(0.01);

% ----------OPERATION MENU COMMANDS-----------------------------------
% --------------------------------------------------------------------
function Run_MenuItem_Callback(hObject, eventdata, handles)

%----------------------
set(handles.Run_MenuItem, 'Enable', 'off');
set(handles.Pause_MenuItem, 'Enable', 'on');
set(handles.Resume_MenuItem, 'Enable', 'off');
set(handles.Stop_MenuItem, 'Enable', 'on');
handles.state = 'Running';
set(handles.text1, 'String', 'Loading...', 'ForegroundColor', [0 0 0]);
drawnow;
%information about program state is transfered trough the text1 field

geometry = handles.param.geometry;
bc = handles.param.bc;
time = handles.param.time;
p_struct = handles.param.p_struct;
beam_struct = handles.param.beam_struct;

beam_n_sp = length(beam_struct);
n_sp = length(p_struct);

geometry = calc_grid_step(geometry, bc);
handles.param.geometry = geometry;

k = 0;

%compose a structure that contains information about particles properties
%(e.g. name, mass, charge, type of interaction with electrodes)
for i = 1:n_sp
    p_internal_prop(i) = p_struct(i).internal_prop;
end
if beam_n_sp > 0
 p_internal_prop(end+1:end+beam_n_sp) = beam_struct.internal_prop;
end
%-------------------------------

%the main variables of simulation program are created; these variables (X,
%Y, VX, VY) contain information about coordinates and velocities of 
%particles in each moment of time; each variables is structure of (n_sp + beam_n_sp)
%length and this structure has the only field - X.coord, Y.coord,
%VX.velocity, VY.velocity, F.free; such data organization allows us to store in one
%variable (X,Y,VX,VY) arrays of different lengthes

global X Y VX VY F
for i = 1:n_sp
    X(i).coord = ones(1, p_struct(i).init_load_param.n_p_max,'double');
    Y(i).coord = ones(1, p_struct(i).init_load_param.n_p_max,'double');
    VX(i).velocity = ones(1, p_struct(i).init_load_param.n_p_max,'double');
    VY(i).velocity = ones(1, p_struct(i).init_load_param.n_p_max,'double');
    F(i).free = logical(zeros(1, p_struct(i).init_load_param.n_p_max));
end

for i = 1:beam_n_sp
    X(i+n_sp).coord = ones(1, beam_struct(i).inject_param.n_p_max,'double');
    Y(i+n_sp).coord = ones(1, beam_struct(i).inject_param.n_p_max,'double');
    VX(i+n_sp).velocity = ones(1, beam_struct(i).inject_param.n_p_max,'double');
    VY(i+n_sp).velocity = ones(1, beam_struct(i).inject_param.n_p_max,'double');
    F(i+n_sp).free = logical(zeros(1, beam_struct(i).inject_param.n_p_max));
end

res = load_particles_2(p_struct(1:n_sp), geometry);

time.t_current = time.t_start;
handles.param.time = time;

set(handles.save_sys_state_item, 'Enable', 'on');
set(handles.Current_Time, 'String', num2str(handles.param.time.t_current, '%2.1e'));

set(handles.text1, 'String', 'Running', 'ForegroundColor', [0 1 0]);


guidata(hObject, handles);


%---------------LOG FILE------------------------------------------------
logfile = handles.logfile;

logfile = [logfile, 'Started at ' ,datestr(now), '.\n'];
logfile = [logfile, '*********PARAMETERS**********\n'];
logfile = [logfile, '----------Geometry----------\n',...
    'x_size = ', num2str(handles.param.geometry.x_size), ' m,\n',...
    'y_size = ', num2str(handles.param.geometry.y_size), ' m,\n',...
    'Nx = ', num2str(handles.param.geometry.ngx), ',\n',...
    'Ny = ', num2str(handles.param.geometry.ngy), '.\n'];
logfile = [logfile, '\n'];
logfile = [logfile, '--------Boundary conditions-------\n',...
    'bc type on top/bottom electrodes: ', handles.param.bc.y_type, ',\n',...
    'bc type on left/right electrodes: ', handles.param.bc.x_type, ',\n',...
    'value on top electrode: ', num2str(handles.param.bc.top_value), ',\n',...
    'value on bottom electrode: ', num2str(handles.param.bc.bottom_value), ',\n',...   
    'value on left electrode: ', num2str(handles.param.bc.left_value), ',\n',...
    'value on right electrode: ', num2str(handles.param.bc.bottom_value), '.\n'];   
logfile = strcat(logfile, '\n');
logfile = strcat(logfile, '----------Particles----------\n');
for j = 1:n_sp
    logfile = [logfile, '   Specie ', num2str(j),':\n', ...
        '      name: ', handles.param.p_struct(j).internal_prop.name, ',\n', ...
        '      charge = ', num2str(handles.param.p_struct(j).internal_prop.charge), ' qe,\n', ...
        '      mass = ', num2str(handles.param.p_struct(j).internal_prop.mass), ' me,\n', ...
        '      lambda = ', num2str(handles.param.p_struct(j).internal_prop.lambda, '%2.1e'), ' 1/m,\n', ...
        '      left/right electrodes interaction type: ', handles.param.p_struct(j).internal_prop.x_interact_type, ',\n', ...
        '      top/bottom electrodes interaction type: ', handles.param.p_struct(j).internal_prop.y_interact_type, ',\n', ...
        '      number of particles = ', num2str(handles.param.p_struct(j).init_load_param.n_p, '%2.1e'), ',\n',...
        '      max number of particles = ' , num2str(handles.param.p_struct(j).init_load_param.n_p_max, '%2.1e'), ',\n',...
        '      function of spatial distribution: ' , func2str(handles.param.p_struct(j).init_load_param.spatial_fun.handle), ',\n',...
        '      parameters of spatial distribution: ', num2str(handles.param.p_struct(j).init_load_param.spatial_fun.param, '%2.1e %2.1e'), ' 1/m^3,\n',...
        '      temperature = ', num2str(handles.param.p_struct(j).init_load_param.velocity_fun.param), ' eV.\n'];
    rho = charge_weighting_5(p_internal_prop, geometry, bc, [j]);
    logfile = [logfile, '     Mean density of specie ', num2str(j), ' = ', num2str(sum(sum(rho))/geometry.ngx/geometry.ngy/handles.const.eq,'%4.3e'), '.\n'];
end
logfile = [logfile, '\n'];
logfile = [logfile, '----------Beam----------\n'];
if beam_n_sp > 0
    for j = 1:beam_n_sp
        logfile = [logfile, '   Beam specie ', num2str(j),':\n', ...
        '      name: ', handles.param.beam_struct(j).internal_prop.name, ',\n', ...
        '      charge = ', num2str(handles.param.beam_struct(j).internal_prop.charge), ' qe,\n', ...
        '      mass = ', num2str(handles.param.beam_struct(j).internal_prop.mass), ' me,\n', ...
        '      lambda = ', num2str(handles.param.beam_struct(j).internal_prop.lambda, '%2.1e'), ' 1/m,\n', ...
        '      left/right electrodes interaction type: ', handles.param.beam_struct(j).internal_prop.x_interact_type, ',\n', ...
        '      top/bottom electrodes interaction type: ', handles.param.beam_struct(j).internal_prop.y_interact_type, ',\n', ...
        '      number of injected particles = ', num2str(handles.param.beam_struct(j).inject_param.n_injected, '%2.1e'), ',\n',...
        '      max number of particles = ' , num2str(handles.param.beam_struct(j).inject_param.n_p_max, '%2.1e'), ',\n',...
        '      vx = ' , num2str(handles.param.beam_struct(j).inject_param.vx, '%2.1e'), ' m/s,\n',...
        '      vy = ' , num2str(handles.param.beam_struct(j).inject_param.vy, '%2.1e'), ' m/s,\n',...
        '      density = ' , num2str(handles.param.beam_struct(j).inject_param.density, '%2.1e'), ' 1/m^3,\n',...
        '      beam width = ' , num2str(handles.param.beam_struct(j).inject_param.beam_width), ' Ly,\n',...
        '      beam position = ' , num2str(handles.param.beam_struct(j).inject_param.beam_pos), ' Ly,\n',...
        '      modulation type: ' , handles.param.beam_struct(j).inject_param.mod_type, ',\n',...
        '      modulation frequency = ' , num2str(handles.param.beam_struct(j).inject_param.mod_frqn, '%2.1e'), ' 1/s,\n',...
        '      modulation depth = ' , num2str(handles.param.beam_struct(j).inject_param.mod_depth), '.\n'];
    
    end
else
    logfile = [logfile, 'No beam in the system. \n'];
        
end
logfile = [logfile, '\n'];
logfile = [logfile, '----------Rho back----------\n'];
if handles.param.rho_back_struct.enabled
    logfile = [logfile, ...
    'name: ', handles.param.rho_back_struct.name, ',\n',...
    'charge = ', num2str(handles.param.rho_back_struct.charge), ' qe,\n',...
    'function of spatial distribution: ' , func2str(handles.param.rho_back_struct.profile_fun.handle), ',\n',...
    'parameters of spatial distribution: ', num2str(handles.param.rho_back_struct.profile_fun.param, '%2.1e %2.1e'), ' 1/m^3.\n'];
else
    logfile = [logfile, 'No background charge in the system. \n'];
    
end
logfile = [logfile, '\n'];

logfile = [logfile, '-----------Time-------------\n'];
logfile = [logfile,...
    't_start = ', num2str(handles.param.time.t_start, '%2.1e'), ' s,\n'...
    't_end: ', num2str(handles.param.time.t_end, '%2.1e'), ' s,\n',...
    'dt: ', num2str(handles.param.time.dt, '%2.1e'), ' s.\n'];
logfile = [logfile, '\n'];
logfile = [logfile, '*******END OF PARAMETERS*****\n'];
    
handles.logfile = logfile;

guidata(hObject, handles);

%-------------- END OF LOG FILE--------------------------------------------

% handles = main_routine(handles);
 res = main_routine(hObject, guidata(hObject));


% --------------------------------------------------------------------
function Pause_MenuItem_Callback(hObject, eventdata, handles)
set(handles.text1,'String','Paused', 'ForegroundColor', [0 0 1]);
handles.state = 'Paused';
set(handles.Resume_MenuItem, 'Enable', 'on');
set(handles.Pause_MenuItem, 'Enable', 'off');
    
guidata(hObject,handles);

% --------------------------------------------------------------------
function Resume_MenuItem_Callback(hObject, eventdata, handles)

handles.logfile = [handles.logfile, 'Resumed at ' ,datestr(now), '.\n'];
set(handles.text1,'String','Running', 'ForegroundColor', [0 1 0]);
handles.state = 'Running';
set(handles.Resume_MenuItem, 'Enable', 'off');
set(handles.Pause_MenuItem, 'Enable', 'on');

guidata(hObject,handles);

res = main_routine(hObject, guidata(hObject));


% --------------------------------------------------------------------
function Stop_MenuItem_Callback(hObject, eventdata, handles)

handles.logfile = [handles.logfile, 'Stopped at ' ,datestr(now), '.\n'];

handles.state = 'Stopped';
set(handles.text1, 'String', 'Stopped', 'ForegroundColor', [0 0 0]);
guidata(hObject, handles);

% ----------END OF OPERATION MENU COMMANDS----------------------------
% --------------------------------------------------------------------

function Parameters_MenuItem_Callback(hObject, eventdata, handles)

% ----------PARAMETERS' DIALOGUES-------------------------------------
% --------------------------------------------------------------------
function Geometry_MenuItem_Callback(hObject, eventdata, handles)

handles.param.geometry = geometry_param(handles.param.geometry);
guidata(hObject,handles);

visualization_field_update(hObject, guidata(hObject));

% --------------------------------------------------------------------
function BC_MenuItem_Callback(hObject, eventdata, handles)
handles.param.bc = boundary_conditions(handles.param.bc);
guidata(hObject, handles);

% --------------------------------------------------------------------
function Particles_MenuItem_Callback(hObject, eventdata, handles)
[handles.param.p_struct handles.n_sp] = particles_param(handles.param.p_struct, handles.param.geometry, handles.param.time);

guidata(hObject,handles);

visualization_field_update(hObject, guidata(hObject));

% --------------------------------------------------------------------
function Beam_MenuItem_Callback(hObject, eventdata, handles)
handles.param.beam_struct = beam_param(handles.param.beam_struct, handles.param.geometry, handles.param.time);

guidata(hObject, handles);

visualization_field_update(hObject, guidata(hObject));

% --------------------------------------------------------------------
function Rho_back_MenuItem_Callback(hObject, eventdata, handles)
handles.param.rho_back_struct = rho_back_param(handles.param.rho_back_struct);
guidata(hObject, handles);

% --------------------------------------------------------------------
function Time_MenuItem_Callback(hObject, eventdata, handles)
handles.param.time = time_param(handles.param.time);

%change of time step leads to the change of linear density of beam
%particles. Therefore one needs to correct beam parameters
for j = 1:length(handles.param.beam_struct)
    handles.param.beam_struct(j).internal_prop.lambda = get_beam_lambda...
        (handles.param.beam_struct(j).inject_param, handles.param.geometry, handles.param.time);
end

guidata(hObject,handles);

% ----------END OF PARAMETERS' DIALOGUES------------------------------

% --------------------------------------------------------------------
function About_MenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------

%-------------------------------------------------------
function res = main_routine(hObject, handles)

p_struct = handles.param.p_struct;
beam_struct = handles.param.beam_struct;
rho_back_struct = handles.param.rho_back_struct;
time = handles.param.time;
geometry = handles.param.geometry;
bc = handles.param.bc;
results_saving_struct = handles.param.results_saving_struct;
const = handles.const;

beam_n_sp = length(beam_struct);
n_sp = length(p_struct);

global X Y VX VY F fi ex ey
if rho_back_struct.enabled
    rho_back = load_rho_back(rho_back_struct, geometry, bc)*const.eq;
end

%compose a structure that contains information about particles properties
%(e.g. name, mass, charge, type of interaction with electrodes)
for i = 1:n_sp
    p_internal_prop(i) = p_struct(i).internal_prop;
end
if beam_n_sp > 0
 p_internal_prop(end+1:end+beam_n_sp) = beam_struct.internal_prop;
end
%-------------------------------
tic
tt = 0;
res = 0;
while (floor(time.t_current/time.dt) <= floor(time.t_end/time.dt))

handles = guidata(hObject);
%     switch get(handles.text1, 'String')
      switch handles.state  
        case 'Running'
            if beam_n_sp > 0
                res = beam_injection_2(beam_struct.inject_param, geometry, time);
            end
            rho = charge_weighting_5(p_internal_prop,geometry,bc,1:(n_sp+beam_n_sp));

            if rho_back_struct.enabled
                rho = rho + rho_back;
            end

            res = field_5(rho, geometry, bc);
            res = e_from_fi_2(geometry, bc);

            res = visualization(hObject, guidata(hObject), p_internal_prop);
            
            res = time_advance_5(p_internal_prop, time, geometry, bc);
            
            res = results_saving(handles, p_internal_prop);

            time.t_current = time.t_current + time.dt;
            handles.param.time = time;


            set(handles.Current_Time, 'String', num2str(time.t_current/1e-9, '%5.1f'));
            set(handles.Duration, 'String', num2str(toc - tt, '%2.1f'));
            guidata(hObject,handles);
            drawnow;
            tt = toc;
        case 'Stopped'
            handles.logfile = [handles.logfile, 'Stopped at ' ,datestr(now), '.\n'];
            break;
        case 'Paused'
            handles.logfile = [handles.logfile, 'Paused at ' ,datestr(now), '.\n'];
            break;
    end
end

%-------------------------------------------------------

if (time.t_current > time.t_end)
    handles.logfile = [handles.logfile, 'Ended at ' ,datestr(now), '.\n'];
end
handles.logfile = [handles.logfile, 'Time elapsed  ' ,num2str(tt), '.\n'];
handles.logfile = [handles.logfile, 'Iterations executed ' ,num2str(floor((time.t_current-time.t_start)/time.dt)), '.\n'];
handles.logfile = [handles.logfile, 'Current simulation time ' ,num2str(time.t_current,'%2.1e'), '.\n'];
handles.param.time = time;

% ----------SAVE LOG FILE ------------------------------
if (time.t_current > time.t_end)|strcmp(handles.state,'Stopped')
    path = handles.param.results_saving_struct.path;
    if ~strcmp(path(end),'\')
        path = strcat(path, '\');
    end
    f = fopen([path, 'logfile.txt'], 'w');
    fprintf(f, handles.logfile);
    fclose(f);

    handles.logfile = '';

    %---temporary
%     set(handles.im, 'CData', []);
%     set(handles.axes1, 'Color', get(handles.uipanel1, 'BackgroundColor'));
%     drawnow;
   
%     res = visualization(hObject, guidata(hObject), p_internal_prop);
    
    handles.state = 'Stopped';
    set(handles.text1, 'String', 'Stopped', 'ForegroundColor', [0 0 0]);
    set(handles.Run_MenuItem, 'Enable', 'on');
    set(handles.Pause_MenuItem, 'Enable', 'off');
    set(handles.Resume_MenuItem, 'Enable', 'off');
    set(handles.Stop_MenuItem, 'Enable', 'off');
end

% -----------END OF SAVE LOG FILE ----------------------

guidata(hObject, handles);

% -------------------------------------------------------
function res = results_saving(handles, p_internal_prop)

global X Y VX VY F ex ey fi

results_saving_struct = handles.param.results_saving_struct;
time = handles.param.time;
p_struct = handles.param.p_struct;
beam_struct = handles.param.beam_struct;
geometry = handles.param.geometry;
bc = handles.param.bc;

n_sp = length(p_struct);
beam_n_sp = length(beam_struct);

if results_saving_struct.enabled

    t_current = time.t_current;
    t_start = results_saving_struct.t_start;
    t_end = results_saving_struct.t_end;
    dt = time.dt;
    k = floor((t_current - t_start)/dt);
    n_dt = results_saving_struct.n_dt;
    path = results_saving_struct.path;
    if ~strcmp(path(end),'\')
        path = strcat(path, '\');
    end
    saving_list = results_saving_struct.saving_list;

    if (t_current >= t_start)&(t_current <= t_end)
        if (mod(k,n_dt) == 0)
            for i = 1:length(saving_list)
                switch(saving_list{i})
                    case 'Potential'
                        save(strcat(path,'fi_',num2str(t_current),'.dat'),'fi')
                    case 'RHO, total charge density'
                        rho = charge_weighting_5(p_internal_prop,geometry,bc,1:(n_sp+beam_n_sp));
                        save(strcat(path,'rho_total_',num2str(t_current),'.dat'),'rho')
                    case 'El. field, X component'
                        save(strcat(path,'ex_',num2str(t_current),'.dat'),'ex')
                    case 'El. field, Y component'
                        save(strcat(path,'ey_',num2str(t_current),'.dat'),'ey')
                    case 'El. field, ABSOLUTE val.'
                        e_abs = (ex.^2 + ey.^2).^0.5;
                        save(strcat(path,'e_abs_',num2str(t_current),'.dat'),'e_abs')
                    case 'Distribution function'
                        save(strcat(path,'Df_',num2str(t_current),'.dat'),'X','Y','VX','VY','F')
                    otherwise
                        name_sp = saving_list{i};
                        name_sp = name_sp(find(name_sp==',')+2:end);
                        for i = 1:n_sp
                            if strcmp(name_sp, p_struct(i).internal_prop.name)
                                rho_sp = charge_weighting_5(p_internal_prop,geometry,bc,[i]);
                                save(strcat(path,'rho_',name_sp,'_',num2str(t_current),'.dat'),'rho_sp')
                            end
                        end
                        for i = n_sp+1 : n_sp+beam_n_sp
                            if strcmp(name_sp, beam_struct(i-n_sp).internal_prop.name)
                                rho_sp = charge_weighting_5(p_internal_prop,geometry,bc,[i]);
                                save(strcat(path,'rho_',name_sp,'_',num2str(t_current),'.dat'),'rho_sp')
                            end
                        end
                end
            end
        end
    end
end

res = 1;

% -----------visualizATION--------------------------------------------

function min_value_1_Callback(hObject, eventdata, handles)

tt = str2num(get(handles.min_value_1, 'String'));

if ~isempty(tt)
    k = get(handles.visualize_popup_1, 'Value');
    handles.visual.lim_struct(k).min_value_1 = str2num(get(handles.min_value_1, 'String'));
    guidata(hObject, handles);
else
    set(handles.min_value_1, 'String', num2str(handles.visual.min_value_1, '%2.1e'));
end

function min_value_1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function max_value_1_Callback(hObject, eventdata, handles)

tt = str2num(get(handles.max_value_1, 'String'));

if ~isempty(tt)
    k = get(handles.visualize_popup_1, 'Value');
    handles.visual.lim_struct(k).max_value_1 = str2num(get(handles.max_value_1, 'String'));
    guidata(hObject, handles);
else
    set(handles.max_value_1, 'String', num2str(handles.visual.max_value_1, '%2.1e'));    
end

% handles.visual.max_value_1 = str2num(get(handles.max_value_1, 'String'));
% guidata(hObject, handles);


function max_value_1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function visualize_popup_1_Callback(hObject, eventdata, handles)

list = get(handles.visualize_popup_1,'String');

k = get(handles.visualize_popup_1,'Value');
handles.visual.visualized_value_1 = list{k};

switch handles.visual.lim_struct(k).mode_1
    case 'auto'
        set(handles.limits_popup_1, 'Value', 1);
        set(handles.min_value_1, 'Enable', 'off');
        set(handles.max_value_1, 'Enable', 'off');
    case 'manual'
        set(handles.limits_popup_1, 'Value', 2);
        set(handles.min_value_1, 'Enable', 'on');
        set(handles.max_value_1, 'Enable', 'on');
end

set(handles.min_value_1, 'String', num2str(handles.visual.lim_struct(k).min_value_1));
set(handles.max_value_1, 'String', num2str(handles.visual.lim_struct(k).max_value_1));

guidata(hObject, handles);

function visualize_popup_1_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function limits_popup_1_Callback(hObject, eventdata, handles)
list = get(handles.limits_popup_1,'String');

k = get(handles.visualize_popup_1, 'Value');
handles.visual.lim_struct(k).mode_1 = list{get(handles.limits_popup_1,'Value')};

switch handles.visual.lim_struct(k).mode_1
    case 'auto'
        set(handles.min_value_1, 'Enable', 'off');
        set(handles.max_value_1, 'Enable', 'off');
    case 'manual'
        set(handles.min_value_1, 'Enable', 'on');
        set(handles.max_value_1, 'Enable', 'on');
end

guidata(hObject, handles);


function limits_popup_1_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in colormap_popup_1.
function colormap_popup_1_Callback(hObject, eventdata, handles)
list = get(handles.colormap_popup_1,'String');

handles.visual.colormap_1 = list{get(handles.colormap_popup_1,'Value')};

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function colormap_popup_1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu13.
function popupmenu13_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu13 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu13


% --- Executes during object creation, after setting all properties.
function popupmenu13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu14.
function popupmenu14_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu14 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu14


% --- Executes during object creation, after setting all properties.
function popupmenu14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu11.
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu12.
function popupmenu12_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu12 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu12


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu16.
function popupmenu16_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu16 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu16


% --- Executes during object creation, after setting all properties.
function popupmenu16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu15.
function popupmenu15_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu15 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu15


% --- Executes during object creation, after setting all properties.
function popupmenu15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --------------------------------------------------------------------
function visualization_mode_SelectionChangeFcn(hObject, eventdata, handles)
Tag = get(get(handles.visualization_mode,'SelectedObject'),'Tag');

switch Tag

    case 'one_axes'
        set(handles.uipanel1, 'Visible', 'on');
        set(handles.uipanel8, 'Visible', 'off');
        set(handles.uipanel9, 'Visible', 'off');
        set(handles.uipanel10, 'Visible', 'off');
        
        set(handles.uipanel5, 'Visible', 'on');
        set(handles.uipanel14, 'Visible', 'off');
        set(handles.uipanel15, 'Visible', 'off');
        set(handles.uipanel16, 'Visible', 'off');
        

    case 'three_axes'
        set(handles.uipanel1, 'Visible', 'off');
        set(handles.uipanel8, 'Visible', 'on');
        set(handles.uipanel9, 'Visible', 'on');
        set(handles.uipanel10, 'Visible', 'on');
        
        set(handles.uipanel5, 'Visible', 'of');
        set(handles.uipanel14, 'Visible', 'on');
        set(handles.uipanel15, 'Visible', 'on');
        set(handles.uipanel16, 'Visible', 'on');        

end
guidata(hObject, handles);

function visualization_field_update(hObject, handles)

% handles = guidata(hObject);
geometry = handles.param.geometry;

Tag = get(get(handles.visualization_mode,'SelectedObject'),'Tag');

switch Tag
    case 'one_axes'
        %--------axes position update--------------------------------
        set(gcf,'Units', 'pixels');
        a = get(gcf,'Position');
        width = a(3);
        height = a(4);
        set(gcf,'Units', 'characters');
        figure_ratio = width/height;

        a = get(handles.uipanel1,'Position');
        width1 = a(3);
        height1 = a(4);
        figure_ratio = figure_ratio*width1/height1;

        system_ratio = geometry.x_size/geometry.y_size;

        if system_ratio >= figure_ratio
            set(handles.axes1, 'Position', [0.1 (1 - 0.8*figure_ratio/system_ratio)/2 0.8 0.8*figure_ratio/system_ratio])
        else
            set(handles.axes1, 'Position', [(1 - 0.8*system_ratio/figure_ratio)/2 0.1 0.8*system_ratio/figure_ratio 0.8])
        end
        set(handles.axes1, 'xlim', [0 geometry.x_size], 'ylim', [0 geometry.y_size], 'Color', get(handles.uipanel1, 'BackgroundColor'));

        %------list of visualize_popup update------------------------
        name_specie_list = [];

        for i = 1:length(handles.param.p_struct)
            name_specie_list{end + 1} = ['Rho, ', handles.param.p_struct(i).internal_prop.name];
        end;
        for i = 1:length(handles.param.beam_struct)
            name_specie_list{end + 1} = ['Rho, ', handles.param.beam_struct(i).internal_prop.name];
        end;
                
        name_specie_list{end+1} = 'Rho, total';
        set(handles.visualize_popup_1, 'String', [handles.visual.init_list; name_specie_list']);
        
        
        for i = length(handles.visual.lim_struct)+1:length([handles.visual.init_list; name_specie_list'])
            handles.visual.lim_struct(i).mode_1 = 'auto';
            handles.visual.lim_struct(i).min_value_1 = 0;
            handles.visual.lim_struct(i).max_value_1 = 1;

        end
        
        
        %------limits update------------------------
        k = get(handles.visualize_popup_1, 'Value');
        switch handles.visual.lim_struct(k).mode_1
            case 'auto'
                set(handles.min_value_1, 'Enable', 'off');
                set(handles.max_value_1, 'Enable', 'off');
            case 'manual'
                set(handles.min_value_1, 'Enable', 'on');
                set(handles.max_value_1, 'Enable', 'on');
        end
    case 'three_axes'

end

guidata(hObject, handles);


function res = visualization(hObject, handles, p_internal_prop)

global X Y VX VY F fi ex ey

geometry = handles.param.geometry;
bc = handles.param.bc;
rho_back_struct = handles.param.rho_back_struct;
visual = handles.visual;

switch visual.visualized_value_1
    case 'none'
        var = [];
    case 'Potential'
        var = fi;
    case 'El. field, abs'
        var = (ex.^2 + ey.^2).^0.5;
    case 'El. field, x'
        var = ex;
    case 'El. field, y'
        var = ey;
    case 'El. field, vector'
        var = [];
    case 'Rho, total'
        var = charge_weighting_5(p_internal_prop,geometry,bc,1:length(p_internal_prop));
        if rho_back_struct.enabled
            rho_back = load_rho_back(rho_back_struct, geometry, bc)*handles.const.eq;
            var = var + rho_back;
        end
        
    otherwise
        
        specie = find(strcmp({p_internal_prop(:).name}, visual.visualized_value_1(6:end)));
        
        var = charge_weighting_5(p_internal_prop,geometry,bc, [specie]);

        
        if specie <= length(handles.param.p_struct)
            
            back_struct.profile_fun.handle = handles.param.p_struct(specie).init_load_param.spatial_fun.handle;
            back_struct.profile_fun.param = handles.param.p_struct(specie).init_load_param.spatial_fun.param;
            back_struct.charge = handles.param.p_struct(specie).internal_prop.charge;
            
            mean_rho = load_rho_back(back_struct, geometry, bc)*handles.const.eq;
            var = var - mean_rho;
        end
      
end

%index of visualized variable
k = get(handles.visualize_popup_1, 'Value');

if ~isempty(var)
    handles.im = image([0 geometry.x_size],[0 geometry.y_size], ...
        var,'Parent',handles.axes1,'cDataMapping','scaled');

    if strcmp(visual.lim_struct(k).mode_1, 'auto')
        mmin = min(min(var));
        mmax = max(max(var));
        set(handles.axes1, 'CLim', [mmin mmax]);
        set(handles.min_value_1, 'String', num2str(mmin, '%2.1e'));
        set(handles.max_value_1, 'String', num2str(mmax, '%2.1e'));
        handles.visual.lim_struct(k).min_value_1 = mmin;
        handles.visual.lim_struct(k).max_value_1 = mmax;
    else
        set(handles.axes1, 'CLim', [visual.lim_struct(k).min_value_1 visual.lim_struct(k).max_value_1]);
    end

    colormap(visual.colormap_1);
end

res = 1;
guidata(hObject, handles);


function visual = initialize_visual(handles)

visual.init_list = get(handles.visualize_popup_1, 'String');

visual.lim_struct = [];

for i = 1: length(visual.init_list)
    visual.lim_struct(i).mode_1 = 'auto';
    visual.lim_struct(i).min_value_1 = 0;
    visual.lim_struct(i).max_value_1 = 1;
    
end

list = get(handles.visualize_popup_1,'String');
visual.visualized_value_1 = list{get(handles.visualize_popup_1, 'Value')}; 

% list = get(handles.limits_popup_1,'String');
% visual.limits_mode_1 = list{get(handles.limits_popup_1,'Value')};

list = get(handles.colormap_popup_1,'String');
visual.colormap_1 = list{get(handles.colormap_popup_1,'Value')};

% visual.max_value_1 = str2num(get(handles.max_value_1, 'String'));
% visual.min_value_1 = str2num(get(handles.min_value_1, 'String'));
