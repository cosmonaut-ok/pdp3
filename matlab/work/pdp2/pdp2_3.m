function varargout = pdp2_3(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pdp2_2_OpeningFcn, ...
                   'gui_OutputFcn',  @pdp2_2_OutputFcn, ...
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
function pdp2_2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% Choose default command line output for pdp2_2
handles.output = hObject;

global X Y VX VY F

const.em = 9.1e-31;
const.eq = 1.6e-19;
const.eps0 = 8.85e-12;
const.kb = 1.38e-23;

handles.const = const;

load('init_param.mat', 'param');
handles.param = param;

visualisation_field_update(hObject, handles);

% handles.CurFig = get(0, 'CurrentFigure');
handles.state = 'Stopped';
set(handles.text1,'String','Stopped');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pdp2_2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = pdp2_2_OutputFcn(hObject, eventdata, handles) 
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


handles.param.results_saving_struct = results_saving(name_list, handles.param.results_saving_struct);
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
visualisation_field_update(hObject, handles);
guidata(hObject, handles);

% --------------------------------------------------------------------
function save_sys_state_item_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function load_sys_state_item_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Exit_MenuItem_Callback(hObject, eventdata, handles)

delete(handles.figure1)

%-----------END OF FILE MENU COMMANDS---------------------------------
%---------------------------------------------------------------------


function Operation_MenuItem_Callback(hObject, eventdata, handles)
pause(0.01);


% ----------OPERATION MENU COMMANDS-----------------------------------
% --------------------------------------------------------------------
function Run_MenuItem_Callback(hObject, eventdata, handles)

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

geometry = calc_grid_step(geometry, bc);

k = 0;

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

%the main variables of simulation program are created; these variables (X,
%Y, VX, VY) contain information about coordinates and velocities of 
%particles in each moment of time; each variables is structure of (n_sp + beam_n_sp)
%length and this structure has the only field - X.coord, Y.coord,
%VX.velocity, VY.velocity, F.free; such data organization allows us to store in one
%variable (X,Y,VX,VY) arrays of different lengthes


if strcmp(handles.state,'Stopped')
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
        F(i+n_sp).free = logical(zeros(1, p_struct(i).init_load_param.n_p_max));
    end
    
    res = load_particles_2(p_struct(1:n_sp), geometry);
    
    time.t_current = time.t_start;

end

% % 
% load('-mat','D:\test\Particles1e-010.dat');
% time.t_start = 1e-010;
% k = 1;

const.eq = 1.6e-19;
rho_back = load_rho_back(rho_back_struct, geometry, bc)*const.eq;

handles.state = 'Running';
set(handles.text1,'String','Running');


%--------should be writed in log file---------------------------
disp('rho back')
sum(sum(rho_back))

for j = 1:n_sp

rho_el = charge_weighting_5(p_internal_prop, geometry, bc, [j]);
disp(strcat('specie ', num2str(j)))
sum(sum(rho_el))
end

%---------------------------------------------------------------


for t_cur = time.t_current:time.dt:time.t_end
%     time.t_current = t_cur;
    tic
    switch handles.state
        case 'Running'

            if beam_n_sp > 0
                res = beam_injection_2(beam_struct.inject_param, geometry, time);
            end

            rho = charge_weighting_5(p_internal_prop,geometry,bc,1:(n_sp+beam_n_sp));

            if rho_back_struct.enabled
                rho = rho + rho_back;
            end

            fi = field_3(rho, geometry, bc);

            [ex ey] = e_from_fi(fi, geometry, bc);

            rho_beam = charge_weighting_5(p_internal_prop,geometry,bc,[2]);
            %             rho_el = rho_el + rho_back;
            im = image([0 geometry.x_size],[0 geometry.y_size], ...
                rho_beam,'Parent',handles.axes1,'cDataMapping','scaled');

            colormap('gray');
            drawnow;
            
            res = time_advance_4(p_internal_prop, time, geometry, bc, ex, ey);
            
            %---------SAVING----------------

            if results_saving_struct.enabled

                t_current = time.t_current;
                t_start = results_saving_struct.t_start;
                t_end = results_saving_struct.t_end;
                dt = time.dt;
                n_dt = results_saving_struct.n_dt;
                path = results_saving_struct.path;
                if ~strcmp(path(end),'\')
                    path = strcat(path, '\');
                end
                saving_list = results_saving_struct.saving_list;

                if (t_current >= t_start)&(t_current <= t_end)

                    %     if (mod(floor((t_current - t_start)/dt),n_dt) == 0)
                    if (mod(k,n_dt) == 0)
                        for i = 1:length(saving_list)
                            switch(saving_list{i})
                                case 'Potential'
                                    save(strcat(path,'fi_',num2str(t_current),'.dat'),'fi')
                                case 'RHO, total charge density'
                                    rho_ions = charge_weighting_4(p_internal_prop,geometry,bc,[3]);
                                    rho = rho - rho_beam - rho_ions;
                                    save(strcat(path,'rho_electrons_',num2str(t_current),'.dat'),'rho')
                                case 'El. field, X component'
                                    save(strcat(path,'ex_',num2str(t_current),'.dat'),'ex')
                                case 'El. field, Y component'
                                    save(strcat(path,'ey_',num2str(t_current),'.dat'),'ey')
                                case 'El. field, ABSOLUTE val.'
                                    e_abs = (ex.^2 + ey.^2).^0.5;
                                    save(strcat(path,'e_abs_',num2str(t_current),'.dat'),'e_abs')
                                otherwise
                                    name_sp = saving_list{i};
                                    name_sp = name_sp(find(name_sp==',')+2:end);
                                    for i = 1:n_sp
                                        if strcmp(name_sp, p_struct(i).internal_prop.name)
                                            rho_sp = charge_weighting_4(p_internal_prop,geometry,bc,[i]);
                                            save(strcat(path,'rho_',name_sp,'_',num2str(t_current),'.dat'),'rho_sp')
                                        end
                                    end
                                    for i = n_sp+1 : n_sp+beam_n_sp

                                        if strcmp(name_sp, beam_struct(i-n_sp).internal_prop.name)
                                            rho_sp = charge_weighting_4(p_internal_prop,geometry,bc,[i]);
                                            save(strcat(path,'rho_',name_sp,'_',num2str(t_current),'.dat'),'rho_sp')
                                        end
                                    end
                            end
                        end
                    end
                end
                if mod(k,1000) == 1
                    save(strcat(path,'Particles',num2str(time.t_current),'.dat'),'X','Y','VX','VY','F')
                end                

            end

            %------------END----OF----SAVING---------------

            time.t_current = time.t_current + time.dt;
            k = k + 1;


        case 'Stopped'
            break;
        case 'Paused'
            break;
    end
    toc
end


handles.time = time;
guidata(hObject, handles);

% --------------------------------------------------------------------
function Step_MenuItem_Callback(hObject, eventdata, handles)
handles.rho = charge_weighting(handles.P, handles.Np, handles.ngx, handles.ngy, handles.dx, handles.dy);
[handles.ex handles.ey] = poifft(handles.rho, 1, 1, 1, 1, handles.bcl, handles.bcr, handles.bcu, handles.bcd, handles.dx, handles.dy);
handles.P = time_advance(handles.P,handles.Np,handles.dt,handles.dx,handles.dy,handles.ex,handles.ey);
handles.tCurrent = handles.tCurrent + handles.dt;

imshow(abs(handles.rho))

% imshow(handles.rho)
guidata(hObject, handles);

% --------------------------------------------------------------------
function Pause_MenuItem_Callback(hObject, eventdata, handles)
set(handles.text1,'String','Paused');
handles.state = 'Paused';
guidata(handles.figure1,handles);
% --------------------------------------------------------------------
function Stop_MenuItem_Callback(hObject, eventdata, handles)

handles.state = 'Stopped';
set(handles.text1, 'String', 'Stopped');
guidata(gcf, handles);

% ----------END OF OPERATION MENU COMMANDS----------------------------
% --------------------------------------------------------------------

function Parameters_MenuItem_Callback(hObject, eventdata, handles)

% ----------PARAMETERS' DIALOGUES-------------------------------------
% --------------------------------------------------------------------
function Geometry_MenuItem_Callback(hObject, eventdata, handles)

handles.param.geometry = geometry_param(handles.param.geometry);
guidata(hObject,handles);

visualisation_field_update(hObject, handles);

% --------------------------------------------------------------------
function BC_MenuItem_Callback(hObject, eventdata, handles)
handles.param.bc = boundary_conditions(handles.param.bc);
guidata(hObject, handles);

% --------------------------------------------------------------------
function Particles_MenuItem_Callback(hObject, eventdata, handles)
[handles.param.p_struct handles.n_sp] = particles_param(handles.param.p_struct, handles.param.geometry);
guidata(hObject,handles);

% --------------------------------------------------------------------
function Beam_MenuItem_Callback(hObject, eventdata, handles)
handles.param.beam_struct = beam_param(handles.param.beam_struct, handles.param.geometry, handles.param.time);
guidata(hObject, handles);

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

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles.visualized_value = get(handles.popupmenu1,'Value');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --------------------------------------------------------------------
function About_MenuItem_Callback(hObject, eventdata, handles)

function visualisation_field_update(hObject, handles)

geometry = handles.param.geometry;

set(handles.axes1, 'xlim', [0 geometry.x_size]);
set(handles.axes1, 'ylim', [0 geometry.y_size]);
if (geometry.x_size > geometry.y_size)
    set(handles.axes1, 'Position', [0.05 (1 - 0.9*geometry.y_size/geometry.x_size)/2 ...
                                    0.9 0.9*geometry.y_size/geometry.x_size]);
else
    set(handles.axes1, 'Position', [(1 - 0.9*geometry.x_size/geometry.y_size)/2 0.05 ...
                                    0.9*geometry.x_size/geometry.y_size 0.9]);
end
guidata(hObject, handles);


