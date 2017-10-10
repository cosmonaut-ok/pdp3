function varargout = control_param(varargin)
% CONTROL_PARAM M-file for control_param.fig
%      CONTROL_PARAM, by itself, creates a new CONTROL_PARAM or raises the existing
%      singleton*.
%
%      H = CONTROL_PARAM returns the handle to a new CONTROL_PARAM or the handle to
%      the existing singleton*.
%
%      CONTROL_PARAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTROL_PARAM.M with the given input arguments.
%
%      CONTROL_PARAM('Property','Value',...) creates a new CONTROL_PARAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before control_param_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to control_param_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help control_param

% Last Modified by GUIDE v2.5 30-Aug-2007 13:16:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @control_param_OpeningFcn, ...
                   'gui_OutputFcn',  @control_param_OutputFcn, ...
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


% --- Executes just before control_param is made visible.
function control_param_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to control_param (see VARARGIN)

% Choose default command line output for control_param
% handles.output = hObject;

p_struct = varargin{1};
geometry = varargin{2};
time = varargin{3};

kb = 1.38e-23;
eq = 1.6e-19;
me = 9.1e-31;
eps0 = 8.85e-12;

f_handle = p_struct.init_load_param.spatial_fun.handle;
param = p_struct.init_load_param.spatial_fun.param;



n_x0 = f_handle(0,(geometry.y_size)/2,geometry,param);
n_xL = f_handle(geometry.x_size,(geometry.y_size)/2,geometry,param);
n_middle = f_handle((geometry.x_size)/2,(geometry.y_size)/2,geometry,param);


T = 11600*p_struct.init_load_param.velocity_fun.param;
lambda = p_struct.internal_prop.lambda;



wp_x0 = (n_x0*eq^2/eps0/me)^0.5;
wp_xL = (n_xL*eq^2/eps0/me)^0.5;

Ld_x0 = (kb*T*eps0/n_x0/eq^2)^0.5;
Ld_xL = (kb*T*eps0/n_xL/eq^2)^0.5;

Nd = kb*T*eps0*pi/eq^2/lambda;

vT = (kb*T/me/p_struct.internal_prop.mass)^0.5;

set(handles.edit1, 'String', num2str(1/wp_x0*1e9, '%2.1e'));
set(handles.edit3, 'String', num2str(1/wp_xL*1e9, '%2.1e'));

set(handles.edit2, 'String', num2str(Ld_x0*1e3, '%2.1e'));
set(handles.edit4, 'String', num2str(Ld_xL*1e3, '%2.1e'));

set(handles.edit6, 'String', num2str(geometry.x_size/geometry.ngx*1e3, '%2.1e'));
set(handles.edit7, 'String', num2str(geometry.y_size/geometry.ngy*1e3, '%2.1e'));
set(handles.edit8, 'String', num2str(time.dt*1e9, '%2.1e'));
set(handles.edit9, 'String', num2str(vT, '%2.1e'));

if min(1/wp_x0, 1/wp_xL) > 10*time.dt
    set(handles.text26, 'String', 'OK', 'ForegroundColor', 'green');
else
    set(handles.text26, 'String', 'Not OK', 'ForegroundColor', 'red');
end
    

if min(Ld_x0, Ld_xL) > 0.3*max(geometry.x_size/geometry.ngx, geometry.y_size/geometry.ngy)
    set(handles.text27, 'String', 'OK', 'ForegroundColor', 'green');
else
    set(handles.text27, 'String', 'Not OK', 'ForegroundColor', 'red');
end

if 3*vT < min(geometry.x_size/geometry.ngx/time.dt, geometry.y_size/geometry.ngy/time.dt)
    set(handles.text28, 'String', 'OK', 'ForegroundColor', 'green');
else
    set(handles.text28, 'String', 'Not OK', 'ForegroundColor', 'red');
end

if Nd > 10
    set(handles.text29, 'String', 'OK', 'ForegroundColor', 'green');
else
    set(handles.text29, 'String', 'Not OK', 'ForegroundColor', 'red');
end

uiwait(gcf);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes control_param wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = control_param_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = 1;
delete(gcf);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
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



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in OK_pb.
function OK_pb_Callback(hObject, eventdata, handles)
% hObject    handle to OK_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcf);


