function varargout = beam_param(varargin)

% Last Modified by GUIDE v2.5 01-Mar-2007 23:46:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @beam_param_OpeningFcn, ...
    'gui_OutputFcn',  @beam_param_OutputFcn, ...
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


% --- Executes just before beam_param is made visible.
function beam_param_OpeningFcn(hObject, eventdata, handles, varargin)

handles.beam_struct = varargin{1};
handles.geometry = varargin{2};
handles.time = varargin{3};

handles.backup = varargin{1};

set(handles.n_species,'String',num2str(length(handles.beam_struct)));

load('init_param.mat', 'param');
handles.empty = param.beam_struct;

if length(handles.beam_struct) > 0

    for i = 1:length(handles.beam_struct)
        temp{i} = num2str(i);
    end

    set(handles.specie_select_menu, 'String', temp);
    set(handles.specie_select_menu, 'Value', 1);

else
    set(handles.specie_select_menu, 'String', '0');
    set(handles.specie_select_menu, 'Value', 1);
end
specie_selection_event(hObject, handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes beam_param wait for user response (see UIRESUME)
uiwait(gcf);

% --- Outputs from this function are returned to the command line.
function varargout = beam_param_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.beam_struct;
delete(gcf);


function n_species_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function specie_select_menu_Callback(hObject, eventdata, handles)
specie_selection_event(hObject, handles);

function specie_select_menu_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%---------------------------------------------------------------
%--------------EDIT BEAM PARAMETERS-----------------------------
function name_Callback(hObject, eventdata, handles)

handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.name = get(hObject,'String');
guidata(hObject, handles);

function name_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------CHARGE------------------------------------------
function charge_Callback(hObject, eventdata, handles)
handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.charge = str2double(get(hObject,'String'));
guidata(hObject, handles);

function charge_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------------MASS---------------------------------------------
function mass_Callback(hObject, eventdata, handles)
handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.mass = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function mass_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-----------------X INTERACTION TYPE------------------------------
function int_type_x_Callback(hObject, eventdata, handles)
switch get(hObject, 'Value')
    case 1
        handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.x_interact_type = 'cycling';
    case 2
        handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.x_interact_type = 'reflection';
    case 3
        handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.x_interact_type = 'absorption';
end
guidata(hObject, handles);

function int_type_x_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-----------------Y INTERACTION TYPE------------------------------
function int_type_y_Callback(hObject, eventdata, handles)
switch get(hObject, 'Value')
    case 1
        handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.y_interact_type = 'cycling';
    case 2
        handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.y_interact_type = 'reflection';
    case 3
        handles.beam_struct(get(handles.specie_select_menu,'Value')).internal_prop.y_interact_type = 'absorption';
end

guidata(hObject, handles);

function int_type_y_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%------------------N INJECTED------------------------------------
function n_injected_Callback(hObject, eventdata, handles)
specie_num = get(handles.specie_select_menu,'Value');
handles.beam_struct(specie_num).inject_param.n_injected = str2double(get(hObject,'String'));

handles.beam_struct(specie_num).internal_prop.lambda = get_beam_lambda(handles.beam_struct(specie_num).inject_param, handles.geometry, handles.time);
set(handles.lambda,'String',num2str(handles.beam_struct(specie_num).internal_prop.lambda,'%3.2e'));
guidata(hObject, handles);

function n_injected_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-----------------N_P_MAX---------------------------------------
function n_p_max_Callback(hObject, eventdata, handles)
handles.beam_struct(get(handles.specie_select_menu,'Value')).inject_param.n_p_max = str2double(get(hObject,'String'));
guidata(hObject, handles);

function n_p_max_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%---------------------VX----------------------------------------
function vx_Callback(hObject, eventdata, handles)
specie_num = get(handles.specie_select_menu,'Value');
handles.beam_struct(specie_num).inject_param.vx = str2double(get(hObject,'String'));

handles.beam_struct(specie_num).internal_prop.lambda = get_beam_lambda(handles.beam_struct(specie_num).inject_param, handles.geometry, handles.time);
set(handles.lambda,'String',num2str(handles.beam_struct(specie_num).internal_prop.lambda,'%3.2e'));
guidata(hObject, handles);

function vx_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%--------------------VY-----------------------------------------
function vy_Callback(hObject, eventdata, handles)
specie_num = get(handles.specie_select_menu,'Value');
handles.beam_struct(specie_num).inject_param.vy = str2double(get(hObject,'String'));

handles.beam_struct(specie_num).internal_prop.lambda = get_beam_lambda(handles.beam_struct(specie_num).inject_param, handles.geometry, handles.time);
set(handles.lambda,'String',num2str(handles.beam_struct(specie_num).internal_prop.lambda,'%3.2e'));
guidata(hObject, handles);

function vy_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%---------------------DENSITY-----------------------------------
function density_Callback(hObject, eventdata, handles)
specie_num = get(handles.specie_select_menu,'Value');
handles.beam_struct(specie_num).inject_param.density = str2double(get(hObject,'String'));

handles.beam_struct(specie_num).internal_prop.lambda = get_beam_lambda(handles.beam_struct(specie_num).inject_param, handles.geometry, handles.time);
set(handles.lambda,'String',num2str(handles.beam_struct(specie_num).internal_prop.lambda,'%3.2e'));
guidata(hObject, handles);

function density_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------BEAM WIDTH-------------------------------------------
function beam_width_Callback(hObject, eventdata, handles)
specie_num = get(handles.specie_select_menu,'Value');
handles.beam_struct(specie_num).inject_param.beam_width = str2double(get(hObject,'String'));

handles.beam_struct(specie_num).internal_prop.lambda = get_beam_lambda(handles.beam_struct(specie_num).inject_param, handles.geometry, handles.time);
set(handles.lambda,'String',num2str(handles.beam_struct(specie_num).internal_prop.lambda,'%3.2e'));
guidata(hObject, handles);

function beam_width_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%----------BEAM POSITION---------------------------------------
function beam_pos_Callback(hObject, eventdata, handles)
handles.beam_struct(get(handles.specie_select_menu,'Value')).inject_param.beam_pos = str2double(get(hObject,'String'));
guidata(hObject, handles);

function beam_pos_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%------------MODULATION TYPE-----------------------------------
function mod_type_Callback(hObject, eventdata, handles)
switch get(hObject, 'Value')
    case 1
        handles.beam_struct(get(handles.specie_select_menu,'Value')).inject_param.mod_type = 'none';
        set(handles.mod_depth, 'Enable', 'off');
        set(handles.mod_frqn,  'Enable', 'off');
    case 2
        handles.beam_struct(get(handles.specie_select_menu,'Value')).inject_param.mod_type = 'density';
        set(handles.mod_depth, 'Enable', 'on');
        set(handles.mod_frqn,  'Enable', 'on');
    case 3
        handles.beam_struct(get(handles.specie_select_menu,'Value')).inject_param.mod_type = 'velocity';
        set(handles.mod_depth, 'Enable', 'on');
        set(handles.mod_frqn,  'Enable', 'on');
end

function mod_type_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%---------------------MODULATION FREQUENCY----------------------
function mod_frqn_Callback(hObject, eventdata, handles)
handles.beam_struct(get(handles.specie_select_menu,'Value')).inject_param.mod_frqn = str2double(get(hObject,'String'));
guidata(hObject, handles);

function mod_frqn_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%--------------------MODULATION DEPTH---------------------------
function mod_depth_Callback(hObject, eventdata, handles)
handles.beam_struct(get(handles.specie_select_menu,'Value')).inject_param.mod_depth = str2double(get(hObject,'String'));
guidata(hObject, handles);

function mod_depth_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%------------------INJECTION ELECTRODE (disabled)--------------
function popupmenu7_Callback(hObject, eventdata, handles)

function popupmenu7_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-----------------BEAM ELECTRONS TEMPERATURE (disabled) -------
function TeV_Callback(hObject, eventdata, handles)
handles.beam_struct(get(handles.specie_select_menu,'Value')).TeV = str2double(get(hObject,'String'));
guidata(hObject, handles);

function TeV_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%-----------VELOCITY DISTRIBUTION (disabled)-------------------
function popupmenu4_Callback(hObject, eventdata, handles)

function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
%-----------------END OF EDIT BEAM PARAMETERS-----------------
%-------------------------------------------------------------

% -----------------pushbutton1 (OK)---------------------------
function pushbutton1_Callback(hObject, eventdata, handles)
uiresume(gcf);

%------------------pushbutton2 (CANCEL)-----------------------
function pushbutton2_Callback(hObject, eventdata, handles)
handles.beam_struct = handles.backup;
guidata(hObject,handles);
uiresume(gcf);

% ------------------- add_specie (ADD)----------------------
function add_specie_Callback(hObject, eventdata, handles)
len = length(handles.beam_struct);

%update the "n_species" field
set(handles.n_species, 'String', num2str(len+1));

%create/update "beam_struct" structure
if len == 0
    %this if .. else ... end statement is nesessary because empty
    %array is a double array ( a=[]) and it's impossible to assign
    %structure to it (b = struct(); a != b)
    handles.beam_struct = handles.empty;
else
    handles.beam_struct(end+1) = handles.empty;
end

%update "specie_select_menu" - species selection menu
for i = 1:len+1
    temp{i} = num2str(i);
end
set(handles.specie_select_menu, 'String', temp);
set(handles.specie_select_menu, 'Value', 1);

% handles.beam_n_sp = len + 1;
guidata(hObject,handles);
specie_selection_event(hObject, handles);


%---------------- del_specie (DEL) ------------------------------.
function del_specie_Callback(hObject, eventdata, handles)
len = length(handles.beam_struct);
if len > 1
    %update the "n_species" field
    set(handles.n_species, 'String', num2str(len-1));

    %update "beam_struct" structure
    handles.beam_struct = handles.beam_struct(1:end-1);

    %update "specie_select_menu" - species selection menu
    for i = 1:len-1
        temp{i} = num2str(i);
    end
    set(handles.specie_select_menu, 'String', temp);
    set(handles.specie_select_menu, 'Value', 1);

%     handles.beam_n_sp = len - 1;
elseif len == 1
    %update the "n_species" field
    set(handles.n_species, 'String', num2str(len-1));

    %update "beam_struct" structure
    handles.beam_struct = [];

    %update "specie_select_menu"
    set(handles.specie_select_menu, 'String', '0');
    set(handles.specie_select_menu, 'Value', 1);

    %make the fields of "beam_param" dialogue empty and uneditable

%     handles.beam_n_sp = len - 1;
else
    return;
end

guidata(hObject,handles);
specie_selection_event(hObject, handles);

% --- Executes on selection change in specie_select_menu.
function specie_selection_event(hObject, handles)
%internal properties

if length(handles.beam_struct) > 0
    
    set(handles.specie_select_menu, 'Enable', 'on');
    set(handles.name, 'Enable', 'on');
    set(handles.charge, 'Enable', 'on');
    set(handles.mass, 'Enable', 'on');
    set(handles.n_injected, 'Enable', 'on');
    set(handles.n_p_max, 'Enable', 'on');
    set(handles.beam_width, 'Enable', 'on');
    set(handles.beam_pos, 'Enable', 'on');
    
    set(handles.vx, 'Enable', 'on');
    set(handles.vy, 'Enable', 'on');
    set(handles.density, 'Enable', 'on');
    set(handles.mod_type, 'Enable', 'on');
    set(handles.mod_frqn, 'Enable', 'on');
    set(handles.mod_depth, 'Enable', 'on');
    set(handles.int_type_y, 'Enable', 'on');
    set(handles.int_type_x, 'Enable', 'on');
  
    selected_specie = handles.beam_struct(get(handles.specie_select_menu,'Value'));
    
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

    % set(handles.TeV, 'String', num2str(selected_specie.init_load_param.velocity_fun.param(1)));
    set(handles.n_injected, 'String', num2str(selected_specie.inject_param.n_injected, '%3.2e'));
    set(handles.n_p_max, 'String', num2str(selected_specie.inject_param.n_p_max, '%3.2e'));

    set(handles.beam_width, 'String', num2str(selected_specie.inject_param.beam_width));
    set(handles.beam_pos, 'String', num2str(selected_specie.inject_param.beam_pos));

    set(handles.vx, 'String', num2str(selected_specie.inject_param.vx, '%3.2e'));
    set(handles.vy, 'String', num2str(selected_specie.inject_param.vy, '%3.2e'));
    set(handles.density, 'String', num2str(selected_specie.inject_param.density, '%3.2e'));

    switch selected_specie.inject_param.mod_type
        case 'none'
            set(handles.mod_type, 'Value', 1);
            set(handles.mod_depth, 'Enable', 'off');
            set(handles.mod_frqn, 'Enable', 'off');
        case 'density'
            set(handles.mod_type, 'Value', 2);
            set(handles.mod_depth, 'Enable', 'on');
            set(handles.mod_frqn, 'Enable', 'on');
        case 'velocity'
            set(handles.mod_type, 'Value', 3);
            set(handles.mod_depth, 'Enable', 'on');
            set(handles.mod_frqn, 'Enable', 'on');
    end
    set(handles.mod_depth, 'String', num2str(selected_specie.inject_param.mod_depth));
    set(handles.mod_frqn, 'String', num2str(selected_specie.inject_param.mod_frqn, '%3.2e'));
    set(handles.lambda,'String',num2str(selected_specie.internal_prop.lambda,'%3.2e'));
else
    set(handles.specie_select_menu, 'Enable', 'off');
    set(handles.name, 'Enable', 'off', 'String', '');
    set(handles.charge, 'Enable', 'off', 'String', '');
    set(handles.mass, 'Enable', 'off', 'String', '');
    set(handles.lambda, 'String', '');
    set(handles.n_injected, 'Enable', 'off', 'String', '');
    set(handles.n_p_max, 'Enable', 'off', 'String', '');
    set(handles.beam_width, 'Enable', 'off', 'String', '');
    set(handles.beam_pos, 'Enable', 'off', 'String', '');
    set(handles.vx, 'Enable', 'off', 'String', '');
    set(handles.vy, 'Enable', 'off', 'String', '');
    set(handles.density, 'Enable', 'off', 'String', '');
    set(handles.mod_type, 'Enable', 'off');
    set(handles.mod_frqn, 'Enable', 'off', 'String', '');
    set(handles.mod_depth, 'Enable', 'off', 'String', '');
    set(handles.int_type_y, 'Enable', 'off');
    set(handles.int_type_x, 'Enable', 'off');
end

guidata(hObject, handles);

