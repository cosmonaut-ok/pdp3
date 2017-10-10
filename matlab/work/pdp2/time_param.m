function varargout = time_param(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @time_param_OpeningFcn, ...
                   'gui_OutputFcn',  @time_param_OutputFcn, ...
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


% --- Executes just before time_param is made visible.
function time_param_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for time_param
handles.time = varargin{1};
handles.backup = varargin{1};

set(handles.edit1,'String',handles.time.t_start*1e9);
set(handles.edit2,'String',handles.time.t_end*1e9);
set(handles.edit3,'String',handles.time.dt*1e9);



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes time_param wait for user response (see UIRESUME)
uiwait(gcf);


% --- Outputs from this function are returned to the command line.
function varargout = time_param_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.time;
delete(gcf);

function edit1_Callback(hObject, eventdata, handles)

handles.time.t_start = str2double(get(hObject,'String'))*1e-9;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit2_Callback(hObject, eventdata, handles)

handles.time.t_end = str2double(get(hObject,'String'))*1e-9;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit3_Callback(hObject, eventdata, handles)

handles.time.dt = str2double(get(hObject,'String'))*1e-9;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)

guidata(hObject, handles);
uiresume(gcf);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)

handles.geometry = handles.backup;
guidata(hObject, handles);
uiresume(gcf);


