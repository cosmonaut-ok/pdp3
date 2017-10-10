function varargout = results_saving_param(varargin)

%results_saving is a dialogue for determining results saving parameters.
%During a simulation process appears a lot of intermediate results -
%electrical field, species concentration, potential etc.
%As a rule the main subject of interest is a dynamic of these physical
%values. 
%The heart of results_saving dialogue is data structure:
% results_saving_struct. This structure has the following fields:
% .enabled - indicates the state of saving: values - 0/1;
% .path - defines the path to directory where the results will be saved;
% .saving_list - cell array of strings which contains names of variables to
% be saved;
% .t_start, .t_end, .n_dt - time parameters which defines time interval within which
%saving is perfomed and time interval between two subsequent savings

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @results_saving_OpeningFcn, ...
                   'gui_OutputFcn',  @results_saving_OutputFcn, ...
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




% --- Executes just before results_saving is made visible.
function results_saving_OpeningFcn(hObject, eventdata, handles, varargin)

% whole_list_str = {{'rho'}, {'el_field_abs'}, {'el_field_x'}, {'el_field_y'}, {'potential'}};

whole_list_str = {'RHO, total charge density', 'El. field, ABSOLUTE val.',...
                  'El. field, X component', 'El. field, Y component', 'Potential', 'Distribution function'};
name_list = varargin{1};
for i = 1:length(name_list)
% the variable temp is used to prevent an error occuring when
% length(name_list) = 1
    name_list{i} = ['Charge density, ', name_list{i}];
end


whole_list_str = {whole_list_str{:}, name_list{:}};
set(handles.whole_list, 'String', whole_list_str);

handles.results_saving_struct = varargin{2};
if (length(handles.results_saving_struct.saving_list) == 0)
    handles.results_saving_struct.saving_list = {''};
end
    

set(handles.t_start, 'String', num2str(handles.results_saving_struct.t_start, '%3.2e'));
set(handles.t_end, 'String', num2str(handles.results_saving_struct.t_end, '%3.2e'));
set(handles.n_dt, 'String', num2str(handles.results_saving_struct.n_dt));
set(handles.path, 'String', handles.results_saving_struct.path);
set(handles.saving_list, 'String', handles.results_saving_struct.saving_list);
set(handles.enabled, 'Value', handles.results_saving_struct.enabled);

if (handles.results_saving_struct.enabled == 0)
    set(handles.path, 'Enable', 'off');
    set(handles.pb_path_select, 'Enable', 'off');
    set(handles.whole_list, 'Enable', 'off');
    set(handles.saving_list, 'Enable', 'off');
    set(handles.pb_add, 'Enable', 'off');
    set(handles.pb_delete, 'Enable', 'off');
    set(handles.t_start, 'Enable', 'off');
    set(handles.t_end, 'Enable', 'off');
    set(handles.n_dt, 'Enable', 'off');
end

handles.backup = varargin{2};

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes results_saving wait for user response (see UIRESUME)
uiwait(gcf);


% --- Outputs from this function are returned to the command line.
function varargout = results_saving_OutputFcn(hObject, eventdata, handles) 
%the length of saving_list must coincide with the number of saved variables
if strcmp(handles.results_saving_struct.saving_list{1}, '')
    handles.results_saving_struct.saving_list = {};
    handles.results_saving_struct.enabled = 0;
end
varargout{1} = handles.results_saving_struct;
delete(gcf);


% --- Executes on selection change in whole_list.
function whole_list_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function whole_list_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in saving_list.
function saving_list_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function saving_list_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pb_add.
function pb_add_Callback(hObject, eventdata, handles)
selected = get(handles.whole_list, 'Value');
whole_list = get(handles.whole_list, 'String');
saving_list = get(handles.saving_list, 'String');
saving_list_len = length(saving_list);


if strcmp(saving_list{1},'')
    saving_list{1} = whole_list{selected};
else
    add_condition = 1;
    for i = 1:saving_list_len
        if strcmp(saving_list{i}, whole_list{selected})
            add_condition = 0;
            break;
        end
    end
    if add_condition
        saving_list{end+1} = whole_list{selected};
    end
end
set(handles.saving_list, 'String', saving_list);     
handles.results_saving_struct.saving_list = saving_list;

guidata(hObject, handles);

        
% --- Executes on button press in pb_delete.
function pb_delete_Callback(hObject, eventdata, handles)

selected = get(handles.saving_list, 'Value');
whole_list = get(handles.whole_list, 'String');
saving_list = get(handles.saving_list, 'String');

if strcmp(saving_list{selected},'')

else
    saving_list_len = length(saving_list);
    for i = selected:saving_list_len-1
        saving_list{i} = saving_list{i+1};
    end
    if (saving_list_len > 1)
        saving_list = {saving_list{1:end-1}};
    else
        saving_list = {''};
    end
    set(handles.saving_list, 'String', saving_list);
    if (selected == saving_list_len)&(selected > 1)
        set(handles.saving_list, 'Value', selected - 1);
    end
end
handles.results_saving_struct.saving_list = saving_list;

guidata(hObject, handles);


function t_start_Callback(hObject, eventdata, handles)
handles.results_saving_struct.t_start = str2num(get(handles.t_start, 'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function t_start_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function t_end_Callback(hObject, eventdata, handles)
handles.results_saving_struct.t_end = str2num(get(handles.t_end, 'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function t_end_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function n_dt_Callback(hObject, eventdata, handles)
handles.results_saving_struct.n_dt = str2num(get(handles.n_dt, 'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function n_dt_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function path_Callback(hObject, eventdata, handles)
if exist(get(handles.path, 'String'))
    handles.results_saving_struct.path = get(handles.path, 'String');
else
    set(handles.path, 'String', handles.results_saving_struct.path);
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in pb_ok.
function pb_ok_Callback(hObject, eventdata, handles)
uiresume(gcf);


% --- Executes on button press in pb_cansel.
function pb_cansel_Callback(hObject, eventdata, handles)
handles.results_saving_struct = handles.backup;
guidata(hObject, handles);
uiresume(gcf);


% --- Executes on button press in enabled.
function enabled_Callback(hObject, eventdata, handles)
handles.results_saving_struct.enabled = get(handles.enabled, 'Value');
if get(handles.enabled, 'Value')
    set(handles.path, 'Enable', 'on');
    set(handles.pb_path_select, 'Enable', 'on');
    set(handles.whole_list, 'Enable', 'on');
    set(handles.saving_list, 'Enable', 'on');
    set(handles.pb_add, 'Enable', 'on');
    set(handles.pb_delete, 'Enable', 'on');
    set(handles.t_start, 'Enable', 'on');
    set(handles.t_end, 'Enable', 'on');
    set(handles.n_dt, 'Enable', 'on');
else
    set(handles.path, 'Enable', 'off');
    set(handles.pb_path_select, 'Enable', 'off');
    set(handles.whole_list, 'Enable', 'off');
    set(handles.saving_list, 'Enable', 'off');
    set(handles.pb_add, 'Enable', 'off');
    set(handles.pb_delete, 'Enable', 'off');
    set(handles.t_start, 'Enable', 'off');
    set(handles.t_end, 'Enable', 'off');
    set(handles.n_dt, 'Enable', 'off');
end
guidata(hObject, handles);



% --- Executes on button press in pb_path_select.
function pb_path_select_Callback(hObject, eventdata, handles)
path = uigetdir(handles.results_saving_struct.path, 'Select results saving directory');
if (path ~= 0)
    handles.results_saving_struct.path = path;
    set(handles.path, 'String', path);
end
guidata(hObject, handles);


