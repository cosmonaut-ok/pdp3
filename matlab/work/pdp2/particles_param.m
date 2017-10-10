function varargout = particles_param(varargin)

% Last Modified by GUIDE v2.5 30-Aug-2007 11:51:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @particles_param_OpeningFcn, ...
                   'gui_OutputFcn',  @particles_param_OutputFcn, ...
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


% --- Executes just before particles_param is made visible.
function particles_param_OpeningFcn(hObject, eventdata, handles, varargin)

handles.p_struct = varargin{1};
handles.geometry = varargin{2};
handles.time = varargin{3};
handles.backup = varargin{1};
handles.n_sp = length(handles.p_struct);
handles.backup_n_sp = length(handles.p_struct);

set(handles.n_species,'String',num2str(length(handles.p_struct)));

     
empty = struct('internal_prop',          {struct('name',            {'none'},    ...
                                                 'mass',            {1},              ...
                                                 'charge',          {-1},             ...
                                                 'lambda',          {[]},             ...
                                                 'x_interact_type', {'reflection'},   ...
                                                 'y_interact_type', {'cycling'}                                                        )},...
                     'init_load_param',  {struct('spatial_fun',     {struct('handle',{@homogene_fun}, 'param',{[1e16 1e16]})}, ... 
                                                 'velocity_fun',    {struct('handle',{[]}, 'param',{[0.5]})},                        ...
                                                 'load_type',       {'uniform_in_space'}, ...
                                                 'n_p',             {4e5},      ...
                                                 'n_p_max',         {5e5}                                                              )});
handles.empty = empty;

for i = 1:length(handles.p_struct)
    temp{i} = num2str(i);
end

set(handles.popupmenu1, 'String', temp);

set(handles.popupmenu1, 'Value', 1);

popupmenu1_Callback(handles.popupmenu1, eventdata, handles);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes particles_param wait for user response (see UIRESUME)
uiwait(gcf);


% --- Outputs from this function are returned to the command line.
function varargout = particles_param_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.p_struct;
varargout{2} = handles.n_sp;
delete(gcf);



function n_species_Callback(hObject, eventdata, handles)
len = length(handles.p_struct);
len_new = str2num(get(hObject,'String'));
if (isfloat(len_new)&(len_new > 0))
    if (len_new > len)
        for j = len+1:len_new
            handles.p_struct(j) = handles.empty;
        end
    else
        handles.p_struct = handles.p_struct(1:len_new);
    end
else
    set(hObject,'String',num2str(len));
end

for i = 1:length(handles.p_struct)
    temp{i} = num2str(i);
end

set(handles.popupmenu1, 'String', temp);
set(handles.popupmenu1, 'Value', 1);

handles.n_sp = str2double(get(hObject,'String'));

guidata(hObject,handles);
    

% --- Executes during object creation, after setting all properties.
function n_species_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function name_Callback(hObject, eventdata, handles)

handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.name = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function charge_Callback(hObject, eventdata, handles)

handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.charge = str2num(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function charge_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function mass_Callback(hObject, eventdata, handles)

handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.mass = str2num(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function mass_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function TeV_Callback(hObject, eventdata, handles)

handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.velocity_fun.param = str2num(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function TeV_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function N_Callback(hObject, eventdata, handles)

handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.n_p = str2num(get(hObject,'String'));

fun = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun;
n_p = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.n_p;
set(handles.text20,'String',num2str(get_lambda(handles.geometry, fun, n_p),'%3.2e'));

handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.lambda = get_lambda(handles.geometry, fun, n_p);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Nmax_Callback(hObject, eventdata, handles)

handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.n_p_max = str2num(get(hObject,'String'));
guidata(hObject, handles);

function Nmax_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)


%internal properties
selected_specie = handles.p_struct(get(handles.popupmenu1,'Value'));

set(handles.name,'String',num2str(selected_specie.internal_prop.name));
set(handles.charge,'String',num2str(selected_specie.internal_prop.charge));
set(handles.mass, 'String', num2str(selected_specie.internal_prop.mass));

switch selected_specie.internal_prop.x_interact_type
    case 'cycling'
        set(handles.int_type_x, 'Value', 1);
    case 'reflection'
        set(handles.int_type_x, 'Value', 2);
    case 'absorption'
        set(handles.int_type_x, 'Value', 3);
end

switch selected_specie.internal_prop.y_interact_type
    case 'cycling'
        set(handles.int_type_y, 'Value', 1);
    case 'reflection'
        set(handles.int_type_y, 'Value', 2);
    case 'absorption'
        set(handles.int_type_y, 'Value', 3);
end

set(handles.TeV, 'String', num2str(selected_specie.init_load_param.velocity_fun.param(1)));
set(handles.N, 'String', num2str(selected_specie.init_load_param.n_p));
set(handles.Nmax, 'String', num2str(selected_specie.init_load_param.n_p_max));


switch func2str(selected_specie.init_load_param.spatial_fun.handle)
    case 'homogene_fun'
        set(handles.sp_distribution, 'Value', 1);
        set(handles.n0, 'String', num2str(selected_specie.init_load_param.spatial_fun.param(1),'%3.2e'));
        set(handles.n1, 'String', num2str(selected_specie.init_load_param.spatial_fun.param(1),'%3.2e'));
        set(handles.n1, 'Enable', 'off');
        
    case 'linear_fun'
        set(handles.sp_distribution, 'Value', 2);
        set(handles.n0, 'String', num2str(selected_specie.init_load_param.spatial_fun.param(1),'%3.2e'));
        set(handles.n1, 'String', num2str(selected_specie.init_load_param.spatial_fun.param(2),'%3.2e'));
%     otherwise
%         set(handles.sp_distribution, 'Value', 3);  

end

fun = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun;
n_p = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.n_p;
set(handles.text20,'String',num2str(get_lambda(handles.geometry, fun, n_p),'%3.2e'));

handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.lambda = get_lambda(handles.geometry, fun, n_p);


guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in int_type_x.
function int_type_x_Callback(hObject, eventdata, handles)

switch get(hObject, 'Value')
    case 1
        handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.x_interact_type = 'cycling';
    case 2
        handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.x_interact_type = 'reflection';
    case 3
        handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.x_interact_type = 'absorption';
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function int_type_x_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in int_type_y.
function int_type_y_Callback(hObject, eventdata, handles)

switch get(hObject, 'Value')
    case 1
        handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.y_interact_type = 'cycling';
    case 2
        handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.y_interact_type = 'reflection';
    case 3
        handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.y_interact_type = 'absorption';
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function int_type_y_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)

function popupmenu4_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in sp_distribution.
function sp_distribution_Callback(hObject, eventdata, handles)

switch get(hObject,'Value')
    case 1
        handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.handle = @homogene_fun;
        handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.param(2) = ...
            handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.param(1);

        set(handles.n1, 'String',num2str(handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.param(1),'%3.2e'));
        set(handles.n1, 'Enable', 'off');
    case 2
        handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.handle = @linear_fun;
%         handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.param = [0 0];
%         set(handles.n0, 'String','0');
%         set(handles.n1, 'String','0');
        set(handles.n1, 'Enable', 'on');
    case 3
        handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.handle = @circle_fun;
        handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.param = [];
        set(handles.n0, 'Enable', 'off');
        set(handles.n1, 'Enable', 'off');
        
end
fun = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun;
n_p = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.n_p;
set(handles.text20,'String',num2str(get_lambda(handles.geometry, fun, n_p),'%3.2e'));

handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.lambda = get_lambda(handles.geometry, fun, n_p);

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function sp_distribution_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
uiresume(gcf);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
handles.p_struct = handles.backup;
handles.n_sp = handles.backup_n_sp;
guidata(hObject,handles);
uiresume(gcf);


function n0_Callback(hObject, eventdata, handles)

handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.param(1) = str2num(get(hObject,'String'));

fun = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun;
n_p = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.n_p;
set(handles.text20,'String',num2str(get_lambda(handles.geometry, fun, n_p),'%3.2e'));

handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.lambda = get_lambda(handles.geometry, fun, n_p);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function n0_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function n1_Callback(hObject, eventdata, handles)

handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun.param(2) = str2num(get(hObject,'String'));
fun = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.spatial_fun;
n_p = handles.p_struct(get(handles.popupmenu1,'Value')).init_load_param.n_p;
set(handles.text20,'String',num2str(get_lambda(handles.geometry, fun, n_p),'%3.2e'));

handles.p_struct(get(handles.popupmenu1,'Value')).internal_prop.lambda = get_lambda(handles.geometry, fun, n_p);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function n1_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end




% --- Executes on button press in Control_param_pb.
function Control_param_pb_Callback(hObject, eventdata, handles)
% hObject    handle to Control_param_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = guidata(hObject);
res = control_param(handles.p_struct(get(handles.popupmenu1,'Value')), handles.geometry, handles.time);







