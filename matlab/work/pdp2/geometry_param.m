function varargout = geometry_param(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @geometry_param_OpeningFcn, ...
                   'gui_OutputFcn',  @geometry_param_OutputFcn, ...
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


% --- Executes just before geometry_param is made visible.
function geometry_param_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for geometry_param
handles.geometry = varargin{1};
handles.backup = varargin{1};

set(handles.edit1,'String',handles.geometry.x_size*100);
set(handles.edit2,'String',handles.geometry.y_size*100);
set(handles.edit3,'String',handles.geometry.ngx);
set(handles.edit4,'String',handles.geometry.ngy);


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes geometry_param wait for user response (see UIRESUME)
uiwait(gcf);


% --- Outputs from this function are returned to the command line.
function varargout = geometry_param_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.geometry;
delete(gcf);

function edit1_Callback(hObject, eventdata, handles)

handles.geometry.x_size = str2double(get(hObject,'String'))/100;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit2_Callback(hObject, eventdata, handles)

handles.geometry.y_size = str2double(get(hObject,'String'))/100;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit3_Callback(hObject, eventdata, handles)

handles.geometry.ngx = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit4_Callback(hObject, eventdata, handles)

handles.geometry.ngy = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)

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


