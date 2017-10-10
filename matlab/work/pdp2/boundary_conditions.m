function varargout = boundary_conditions(varargin)


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @boundary_conditions_OpeningFcn, ...
                   'gui_OutputFcn',  @boundary_conditions_OutputFcn, ...
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


% --- Executes just before boundary_conditions is made visible.
function boundary_conditions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.

handles.bc = varargin{1};
handles.backup = varargin{1};

switch handles.bc.x_type
    case 'periodic'
        set(handles.x_el_type,'SelectedObject',handles.x_per);
    case 'dirichlet'
        set(handles.x_el_type,'SelectedObject',handles.x_dir);
    case 'neumann'
        set(handles.x_el_type,'SelectedObject',handles.x_neu);
end

switch handles.bc.y_type
    case 'periodic'
        set(handles.y_el_type,'SelectedObject',handles.y_per);
    case 'dirichlet'
        set(handles.y_el_type,'SelectedObject',handles.y_dir);
    case 'neumann'
        set(handles.y_el_type,'SelectedObject',handles.y_neu);
end

set(handles.left_el_value,'String',handles.bc.left_value);
set(handles.right_el_value,'String',handles.bc.right_value);
set(handles.top_el_value,'String',handles.bc.top_value);
set(handles.bottom_el_value,'String',handles.bc.bottom_value);

x_el_type_SelectionChangeFcn(handles.x_el_type, eventdata, handles);
y_el_type_SelectionChangeFcn(handles.y_el_type, eventdata, handles);

guidata(hObject, handles);
% UIWAIT makes boundary_conditions wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = boundary_conditions_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.bc;
delete(gcf);


% --- Executes when EditBox left_el_value is changed
function left_el_value_Callback(hObject, eventdata, handles)
handles.bc.left_value = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function left_el_value_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes when EditBox right_el_value is changed
function right_el_value_Callback(hObject, eventdata, handles)

handles.bc.right_value = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function right_el_value_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes when EditBox top_el_value is changed
function top_el_value_Callback(hObject, eventdata, handles)

handles.bc.top_value = str2double(get(hObject,'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function top_el_value_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes when EditBox bottom_el_value is changed
function bottom_el_value_Callback(hObject, eventdata, handles)

handles.bc.bottom_value = str2double(get(hObject,'String'));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function bottom_el_value_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in OK_button.
function OK_button_Callback(hObject, eventdata, handles)
guidata(hObject, handles);
uiresume(gcf);

% --- Executes on button press in Cancel_button.
function Cancel_button_Callback(hObject, eventdata, handles)

handles.bc = handles.backup;
guidata(hObject, handles);
uiresume(gcf);

% --------------------------------------------------------------------
function y_el_type_SelectionChangeFcn(hObject, eventdata, handles)

Tag = get(get(handles.y_el_type,'SelectedObject'),'Tag');

switch Tag
    
    case 'y_per'
        set(handles.top_el_value,'Enable','off');
        set(handles.bottom_el_value,'Enable','off');
        set(handles.top_el_unit,'String','');
        set(handles.bottom_el_unit,'String',''); 
        handles.bc.y_type = 'periodic';
        
    case 'y_dir'
        set(handles.top_el_value,'Enable','on');
        set(handles.bottom_el_value,'Enable','on');
        set(handles.top_el_unit,'String','V');
        set(handles.bottom_el_unit,'String','V');   
        handles.bc.y_type = 'dirichlet';
            
    case 'y_neu'
        set(handles.top_el_value,'Enable','on');
        set(handles.bottom_el_value,'Enable','on');
        set(handles.top_el_unit,'String','V/m');
        set(handles.bottom_el_unit,'String','V/m');  
        handles.bc.y_type = 'neumann';
        
end
guidata(hObject, handles);


% --------------------------------------------------------------------
function x_el_type_SelectionChangeFcn(hObject, eventdata, handles)

Tag = get(get(handles.x_el_type,'SelectedObject'),'Tag');

switch Tag
    
    case 'x_per'
        set(handles.left_el_value,'Enable','off');
        set(handles.right_el_value,'Enable','off');
        set(handles.left_el_unit,'String','');
        set(handles.right_el_unit,'String',''); 
        handles.bc.x_type = 'periodic';
        
    case 'x_dir'
        set(handles.left_el_value,'Enable','on');
        set(handles.right_el_value,'Enable','on');
        set(handles.left_el_unit,'String','V');
        set(handles.right_el_unit,'String','V'); 
        handles.bc.x_type = 'dirichlet';
            
    case 'x_neu'
        set(handles.left_el_value,'Enable','on');
        set(handles.right_el_value,'Enable','on');
        set(handles.left_el_unit,'String','V/m');
        set(handles.right_el_unit,'String','V/m');  
        handles.bc.x_type = 'neumann';
        
end
guidata(hObject, handles);
