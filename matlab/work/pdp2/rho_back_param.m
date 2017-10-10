function varargout = rho_back_param(varargin)
% RHO_BACK_PARAM M-file for rho_back_param.fig
%      RHO_BACK_PARAM, by itself, creates a new RHO_BACK_PARAM or raises the existing
%      singleton*.
%
%      H = RHO_BACK_PARAM returns the handle to a new RHO_BACK_PARAM or the handle to
%      the existing singleton*.
%
%      RHO_BACK_PARAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RHO_BACK_PARAM.M with the given input arguments.
%
%      RHO_BACK_PARAM('Property','Value',...) creates a new RHO_BACK_PARAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rho_back_param_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rho_back_param_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help rho_back_param

% Last Modified by GUIDE v2.5 23-Jul-2006 14:22:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rho_back_param_OpeningFcn, ...
                   'gui_OutputFcn',  @rho_back_param_OutputFcn, ...
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


% --- Executes just before rho_back_param is made visible.
function rho_back_param_OpeningFcn(hObject, eventdata, handles, varargin)

handles.backup = varargin{1};
handles.rho_back_struct = varargin{1};

set(handles.name, 'String', handles.rho_back_struct.name);
set(handles.charge, 'String', num2str(handles.rho_back_struct.charge));

profile_array = get(handles.profile_fun, 'String');

cur_profile_value = [];

%Concentration profiles are represented in the dialogues with their common
%names: homogene, linear etc.
%Hanle of function consists of common name of profile (linear) plus
%addition (_fun): linear_fun, homogene_fun

handle_str = func2str(handles.rho_back_struct.profile_fun.handle);

for i = 1:length(profile_array)
    if strcmp(profile_array{i}, handle_str(1:end-4))
        cur_profile_value(end+1) = i;
    end
end

switch length(cur_profile_value)
    case 0
        disp('Warning: there is no option for concentration profile defined in rho_back_struct. Concentration profile is taken as homogeneous')
        set(handles.profile_fun, 'Value', 1);
    case 1
        set(handles.profile_fun, 'Value', cur_profile_value);
    otherwise
        disp('Warning: list of options for concentration profile is incorrect!')
        set(handles.profile_fun, 'Value', cur_profile_value(1));
end

switch get(handles.profile_fun, 'Value')
    case 1
        set(handles.n0, 'String', num2str(handles.rho_back_struct.profile_fun.param(1),'%3.2e'));
        set(handles.n1, 'String', num2str(handles.rho_back_struct.profile_fun.param(1),'%3.2e'));
        set(handles.n1, 'Enable', 'off');
    case 2
        set(handles.n0, 'String', num2str(handles.rho_back_struct.profile_fun.param(1),'%3.2e'));
        set(handles.n1, 'String', num2str(handles.rho_back_struct.profile_fun.param(2),'%3.2e'));
end

if (handles.rho_back_struct.enabled == 0)
    set(handles.name, 'Enable', 'off');
    set(handles.charge, 'Enable', 'off');
    set(handles.profile_fun, 'Enable', 'off');
    set(handles.n0, 'Enable', 'off');
    set(handles.n1, 'Enable', 'off');
end

set(handles.enabled, 'Value', handles.rho_back_struct.enabled);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rho_back_param wait for user response (see UIRESUME)
uiwait(gcf);


% --- Outputs from this function are returned to the command line.
function varargout = rho_back_param_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.rho_back_struct;
delete(gcf);


function name_Callback(hObject, eventdata, handles)

handles.rho_back_struct.name = get(handles.name, 'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function name_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function charge_Callback(hObject, eventdata, handles)

handles.rho_back_struct.charge = str2num(get(handles.charge, 'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function charge_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in profile_fun.
function profile_fun_Callback(hObject, eventdata, handles)

profile_array = get(handles.profile_fun, 'String');
handles.rho_back_struct.profile_fun.handle = str2func(strcat(profile_array(get(handles.profile_fun, 'Value')),'_fun'));

switch get(handles.profile_fun, 'Value')
    case 1
        set(handles.n0, 'Enable', 'on');
        set(handles.n1, 'String', num2str(handles.rho_back_struct.profile_fun.param(1),'%3.2e'));
        set(handles.n1, 'Enable', 'off');
    case 2
        set(handles.n0, 'Enable', 'on');
        set(handles.n1, 'Enable', 'on');
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function profile_fun_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function n0_Callback(hObject, eventdata, handles)

handles.rho_back_struct.profile_fun.param(1) = str2num(get(handles.n0, 'String'));
if (get(handles.profile_fun, 'Value')==1)
    set(handles.n1, 'String', get(handles.n0,'String'));
    handles.rho_back_struct.profile_fun.param(2) = handles.rho_back_struct.profile_fun.param(1);
end
    
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function n0_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function n1_Callback(hObject, eventdata, handles)
handles.rho_back_struct.profile_fun.param(2) = str2num(get(handles.n1, 'String'));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function n1_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
uiresume(gcf);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
handles.rho_back_struct = handles.backup;
guidata(hObject, handles);
uiresume(gcf);







% --- Executes on button press in enabled.
function enabled_Callback(hObject, eventdata, handles)
handles.rho_back_struct.enabled = get(handles.enabled, 'Value');
if get(handles.enabled, 'Value')
    set(handles.name, 'Enable', 'on');
    set(handles.charge, 'Enable', 'on');
    set(handles.profile_fun, 'Enable', 'on');
    set(handles.n0, 'Enable', 'on');
    set(handles.n1, 'Enable', 'on');    
    
else
    set(handles.name, 'Enable', 'off');
    set(handles.charge, 'Enable', 'off');
    set(handles.profile_fun, 'Enable', 'off');
    set(handles.n0, 'Enable', 'off');
    set(handles.n1, 'Enable', 'off');    
    
end

guidata(hObject, handles);


