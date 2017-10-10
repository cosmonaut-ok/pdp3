function varargout = pdp2_2(varargin)

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

const.em = 9.1e-31;
const.eq = 1.6e-19;
const.eps0 = 8.85e-12;
const.kb = 1.38e-23;

handles.const = const;

% load('init_param.mat', 'param')
% 
% handles.param = param;

geometry.x_size = 0.6;
geometry.y_size = 0.1;
geometry.ngx = 1024;
geometry.ngy = 256;

handles.geometry = geometry;

time.dt = 1e-10;
time.t_start = 0.0;
time.t_end = 1e-6;
time.t_current = 0.0;

handles.time = time;

bc.x_type = 'dirichlet';
bc.y_type = 'periodic';

bc.top_value = 0;
bc.bottom_value = 0;
bc.left_value = 0;
bc.right_value = 0;

handles.bc = bc;

%particles parameters

handles.n_sp = 1;
handles.beam_n_sp = 0;


% particles

p_struct(1) = struct('internal_prop',    {struct('name',            {'electrons'},    ...
                                                 'mass',            {1},              ...
                                                 'charge',          {-1},             ...
                                                 'lambda',          {[]},             ...
                                                 'x_interact_type', {'reflection'},   ...
                                                 'y_interact_type', {'cycling'}                                                        )},...
                     'init_load_param',  {struct('spatial_fun',     {struct('handle',{@linear_fun}, 'param',{0.64*[0.5e16 1.5e16]})}, ... 
                                                 'velocity_fun',    {struct('handle',{[]}, 'param',{[0.5]})},                        ...
                                                 'load_type',       {'uniform_in_space'}, ...
                                                 'n_p',             {1e5},      ...
                                                 'n_p_max',         {1e5}                                                              )});
                                                                  
sp_fun = p_struct(1).init_load_param.spatial_fun;
N = p_struct(1).init_load_param.n_p;
p_struct(1).internal_prop.lambda = get_lambda(geometry, sp_fun, N);

rho_back_struct(1) = struct('name', {'ions'},  ...
                            'charge', {1},     ...
                            'profile_fun',  {struct('handle',{@linear_fun}, 'param',{[2.0e14 4.4e14]})},...
                            'enabled', {1});
                                            
n_middle = sp_fun.handle(geometry.x_size/2,geometry.y_size/2,geometry,sp_fun.param);    
TeV = p_struct(1).init_load_param.velocity_fun.param;
Te = TeV*1.16e4;

wp = (n_middle*const.eq^2/const.em/const.eps0)^0.5;
ld = (const.kb*Te*const.eps0/const.eq^2/n_middle)^0.5;
vTe = (const.kb*Te/const.em)^0.5;
                            
beam_struct(1) = struct('internal_prop',  {struct('name',            {'light_electrons'},    ...
                                                 'mass',            {1},                    ...
                                                 'charge',          {-1},                   ...
                                                 'lambda',          {[]},                   ...
                                                 'x_interact_type', {'absorption'},          ...
                                                 'y_interact_type', {'cycling'}                                                   )},...
                        'inject_param',  {struct('vx',             {3e7},                   ... 
                                                 'vy',             {0},                    ...
                                                 'density',         {1e12},                 ...
                                                 'n_injected',      {500},                  ...
                                                 'n_p_max',         {5e5},                  ...
                                                 'mod_type',      {'density'},              ...
                                                 'mod_depth',     {1},                      ...
                                                 'mod_frqn',      {4.7e9},                    ...
                                                 'beam_width',    {1},                    ...
                                                 'beam_pos',      {0.5}                                                         )});
                                             
                                                 
results_saving_struct = struct('t_start', {0}, 't_end', {1e-7}, 'n_dt', {5}, 'path', {'d:\'}, 'enabled', {1},...
                               'saving_list', {{'RHO, total charge density'}});
                           
handles.beam_n_sp = 1;
handles.n_sp = 0;

density = beam_struct(1).inject_param.density;
beam_width = beam_struct(1).inject_param.beam_width;
velocity = sqrt(beam_struct(1).inject_param.vx^2 + beam_struct(1).inject_param.vy^2);

n_injected = beam_struct(1).inject_param.n_injected;
dt = handles.time.dt;
                                             
beam_struct(1).internal_prop.lambda = round(density*beam_width*geometry.y_size*velocity*dt/n_injected);



% l = beam_struct(1).internal_prop.lambda
                 
handles.p_struct = p_struct;
handles.beam_struct = beam_struct;
handles.rho_back_struct = rho_back_struct;
handles.results_saving_struct = results_saving_struct;
handles.visualized_value = 1;


set(handles.axes1,'xlim',[0 handles.geometry.x_size], 'ylim', [0 handles.geometry.y_size], 'color', [0 0 0]);
if (handles.geometry.x_size > handles.geometry.y_size)
    set(handles.axes1, 'Position', [0.05 (1 - 0.9*handles.geometry.y_size/handles.geometry.x_size)/2 ...
                                    0.9 0.9*handles.geometry.y_size/handles.geometry.x_size]);
else
    set(handles.axes1, 'Position', [(1 - 0.9*handles.geometry.x_size/handles.geometry.y_size)/2 0.05 ...
                                    0.9*handles.geometry.x_size/handles.geometry.y_size 0.9]);
end

handles.CurFig = get(0, 'CurrentFigure');
handles.Status = 'Stopped';
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
function Operation_MenuItem_Callback(hObject, eventdata, handles)
pause(0.01);

% --------------------------------------------------------------------
function About_MenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Open_MenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Save_MenuItem_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function Exit_MenuItem_Callback(hObject, eventdata, handles)

delete(handles.figure1)


% --------------------------------------------------------------------
function Run_MenuItem_Callback(hObject, eventdata, handles)

p_struct = handles.p_struct;
beam_struct = handles.beam_struct;
rho_back_struct = handles.rho_back_struct;
time = handles.time;
geometry = handles.geometry;
bc = handles.bc;
results_saving_struct = handles.results_saving_struct;
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
    F(i+n_sp).free = logical(zeros(1, p_struct(i).init_load_param.n_p_max));
end

res = load_particles_2(p_struct(1:n_sp), geometry);
  
% % 
% load('-mat','D:\test\Particles1e-010.dat');
% time.t_start = 1e-010;
% k = 1;


const.eq = 1.6e-19;
rho_back = load_rho_back(rho_back_struct, geometry, bc)*const.eq;


handles.Status = 'Running';
set(handles.text1,'String','Running');


disp('rho back')
sum(sum(rho_back))

for j = 1:n_sp

rho_el = charge_weighting_5(p_internal_prop, geometry, bc, [j]);
disp(strcat('specie ', num2str(j)))
sum(sum(rho_el))
end


for t_cur = time.t_start:time.dt:time.t_end
    time.t_current = t_cur;
    tic
    switch get(handles.text1,'String')
        case 'Running'

            if beam_n_sp > 0
                res = beam_injection_2(beam_struct.inject_param, geometry, time);
            end

            rho = charge_weighting_5(p_internal_prop,geometry,bc,1:(n_sp+beam_n_sp));

            if handles.rho_back_struct.enabled
                rho = rho + rho_back;
            end

            fi = field_3(rho, geometry, bc);

            [ex ey] = e_from_fi(fi, geometry, bc);

            rho_beam = charge_weighting_5(p_internal_prop,geometry,bc,[n_sp+1:n_sp+beam_n_sp]);
            %             rho_el = rho_el + rho_back;
            im = image([0 handles.geometry.x_size],[0 handles.geometry.y_size], ...
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

            % actions for this condition must be rewrited
            %         case 'Paused'
            %             pause(0.01);
            %             continue;
        case 'Stopped'
            %             set(handles.text1,'String','Stopped');
            %               disp(handles.Status);
            break;

    end
    toc
end
%     set(handles.text1,'String','Breaked');


    guidata(hObject, handles);
% start(handles.t);
% while (handles.Running)&(handles.tCurrent < handles.tTotal)
%     handles.rho = charge_weighting(handles.P, handles.Np, handles.ngx, handles.ngy, handles.dx, handles.dy);
%     [handles.ex handles.ey] = poifft(handles.rho, 1, 1, 1, 1, handles.bcl, handles.bcr, handles.bcu, handles.bcd, handles.dx, handles.dy);
%     [handles.P handles.Np] = time_advance(handles.P,handles.Np,handles.dt,handles.xSize,handles.ySize, handles.ngx, handles.ngy,handles.ex,handles.ey);
%     handles.tCurrent = handles.tCurrent + handles.dt;
%     imshow(handles.rho/max(max(abs(handles.rho))), [-1 1])

%     pause(.1)
% end


% --------------------------------------------------------------------
function Step_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Step_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.rho = charge_weighting(handles.P, handles.Np, handles.ngx, handles.ngy, handles.dx, handles.dy);
[handles.ex handles.ey] = poifft(handles.rho, 1, 1, 1, 1, handles.bcl, handles.bcr, handles.bcu, handles.bcd, handles.dx, handles.dy);
handles.P = time_advance(handles.P,handles.Np,handles.dt,handles.dx,handles.dy,handles.ex,handles.ey);
handles.tCurrent = handles.tCurrent + handles.dt;

imshow(abs(handles.rho))

% imshow(handles.rho)
guidata(hObject, handles);


% --------------------------------------------------------------------
function Stop_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Stop_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% stop(handles.t);
% delete(handles.t);
% uiwait(gcf);

% waitfor(handles.text1, 'String');
handles.Status = 'Stopped';
% guidata(handles.figure1, handles);
 set(handles.text1, 'String', 'Stopped');
% 
% drawnow;
% % handles.Stauts= 'Stopped';
% % stop(handles.t);
% 
 guidata(gcf, handles);


%  uiresume(handles.figure1);


% --------------------------------------------------------------------
function Geometry_MenuItem_Callback(hObject, eventdata, handles)

handles.geometry = geometry_param(handles.geometry);
guidata(hObject,handles);
set(handles.axes1, 'xlim', [0 handles.geometry.x_size]);
set(handles.axes1, 'ylim', [0 handles.geometry.y_size]);
if (handles.geometry.x_size > handles.geometry.y_size)
    set(handles.axes1, 'Position', [0.05 (1 - 0.9*handles.geometry.y_size/handles.geometry.x_size)/2 ...
                                    0.9 0.9*handles.geometry.y_size/handles.geometry.x_size]);
else
    set(handles.axes1, 'Position', [(1 - 0.9*handles.geometry.x_size/handles.geometry.y_size)/2 0.05 ...
                                    0.9*handles.geometry.x_size/handles.geometry.y_size 0.9]);
end
guidata(hObject, handles);
%[handles.xSize handles.ySize handles.ngx handles.ngy]   = a;
% --------------------------------------------------------------------
function Time_MenuItem_Callback(hObject, eventdata, handles)
handles.time = time_param(handles.time);
guidata(hObject,handles);


% --------------------------------------------------------------------
function Particles_MenuItem_Callback(hObject, eventdata, handles)
[handles.p_struct handles.n_sp] = particles_param(handles.p_struct, handles.geometry);
guidata(hObject,handles);


% --------------------------------------------------------------------
function Parameters_MenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to Parameters_MenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function BC_MenuItem_Callback(hObject, eventdata, handles)
handles.bc = boundary_conditions(handles.bc);
guidata(hObject, handles);

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

handles.visualized_value = get(handles.popupmenu1,'Value');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
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
function Pause_MenuItem_Callback(hObject, eventdata, handles)
set(handles.text1,'String','Paused');
guidata(handles.figure1,handles);


% --------------------------------------------------------------------
function Beam_MenuItem_Callback(hObject, eventdata, handles)
handles.beam_struct = beam_param(handles.beam_struct, handles.geometry, handles.time);
guidata(hObject, handles);


% --------------------------------------------------------------------
function Rho_back_MenuItem_Callback(hObject, eventdata, handles)
handles.rho_back_struct = rho_back_param(handles.rho_back_struct);
guidata(hObject, handles);


% --------------------------------------------------------------------
function save_sim_param_item_Callback(hObject, eventdata, handles)

[file, path] = uiputfile('*.mat', 'Save simulation parameters as');
if isstr(path)
    param.time = handles.time;
    param.p_struct = handles.p_struct;
    param.beam_struct = handles.beam_struct;
    param.geometry = handles.geometry;
    param.bc = handles.bc;
    param.rho_back_struct = handles.rho_back_struct;
    param.results_saving_struct = handles.results_saving_struct;
    save([path file], 'param')
end


% --------------------------------------------------------------------
function load_sim_param_item_Callback(hObject, eventdata, handles)
[file, path] = uigetfile('*.mat', 'Load simulation parameters');
if isstr(path)
    load([path file], 'param')
    handles.time = param.time;
    handles.p_struct = param.p_struct;
    handles.beam_struct = param.beam_struct;
    handles.geometry = param.geometry;
    handles.bc = param.bc;
    handles.rho_back_struct = param.rho_back_struct;
    handles.results_saving_struct = param.results_saving_struct;
    guidata(hObject, handles);
end


% --------------------------------------------------------------------
function save_sys_state_item_Callback(hObject, eventdata, handles)
% hObject    handle to save_sys_state_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function load_sys_state_item_Callback(hObject, eventdata, handles)
% hObject    handle to load_sys_state_item (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --------------------------------------------------------------------
function results_saving_item_Callback(hObject, eventdata, handles)



if (length(handles.beam_struct) > 0)
    
    for i = 1:length(handles.p_struct)
        p_internal_prop(i) = handles.p_struct(i).internal_prop;
    end
    for i = 1:length(handles.beam_struct)
        beam_internal_prop(i) = handles.beam_struct(i).internal_prop;
    end
    name_list = {p_internal_prop.name, beam_internal_prop.name};
else
    for i = 1:length(handles.p_struct)
        p_internal_prop(i) = handles.p_struct(i).internal_prop;
    end
    name_list = {p_internal_prop.name};  
end

handles.results_saving_struct = results_saving(name_list, handles.results_saving_struct);
guidata(hObject, handles)




