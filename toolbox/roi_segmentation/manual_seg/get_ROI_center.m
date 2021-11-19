function varargout = get_ROI_center(varargin)
% GET_ROI_CENTER MATLAB code for get_ROI_center.fig
%      GET_ROI_CENTER, by itself, creates a new GET_ROI_CENTER or raises the existing
%      singleton*.
%
%      H = GET_ROI_CENTER returns the handle to a new GET_ROI_CENTER or the handle to
%      the existing singleton*.
%
%      GET_ROI_CENTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GET_ROI_CENTER.M with the given input arguments.
%
%      GET_ROI_CENTER('Property','Value',...) creates a new GET_ROI_CENTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before get_ROI_center_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to get_ROI_center_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help get_ROI_center

% Last Modified by GUIDE v2.5 24-May-2021 15:33:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @get_ROI_center_OpeningFcn, ...
                   'gui_OutputFcn',  @get_ROI_center_OutputFcn, ...
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


% --- Executes just before get_ROI_center is made visible.
function get_ROI_center_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to get_ROI_center (see VARARGIN)

% create global variable
global p

p = [];
% get all files from current directory
fprintf('Reading files from local directory\n')
files2load = rdir('*metadata.mat');
p.files2load = {files2load.name}';

% populate menu
set(handles.popupmenu1, 'String', p.files2load);
set(handles.popupmenu1, 'Value', 1);

% set default value as file to load
p.currentfile = get(handles.popupmenu1, 'Value');

handles.axes1.XTick = [];
handles.axes1.YTick = [];

% Choose default command line output for get_ROI_center
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes get_ROI_center wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = get_ROI_center_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global p

fprintf('Drawing ROI center and radious\n')

[XX, YY] = getpts(handles.axes1);

p.roi_center(p.idx, :) = [XX(1) YY(1), p.z];
p.roi_idx(p.idx, :) = str2double(get(handles.edit3, 'String'));
p.roi_border(p.idx, :) = [XX(2) YY(2)];
p.roi_radious(p.idx, :) = sqrt((p.roi_border(p.idx, 1) - p.roi_center(p.idx, 1))^2 ...
    + (p.roi_border(p.idx, 2) - p.roi_center(p.idx, 2))^2);
        
plot_plane(handles)

p.idx = p.idx + 1;

handles.text6.String = ['ROI n (out of ', num2str(numel(p.roi_idx)), ')'];
set(handles.edit3, 'String', num2str(p.idx))

% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global p

plot_plane(handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global p

p.cha = round(get(handles.slider2, 'Value'));

% reset min max to default
set(handles.edit1, 'String', num2str(p.min_max(p.cha, 1)))
set(handles.edit2, 'String', num2str(p.min_max(p.cha, 2)))

plot_plane(handles)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global p
p.currentfile = get(hObject,'Value');

disp(['Loading ' p.files2load{p.currentfile}])

% load file
load(p.files2load{p.currentfile})
p.wDat = wDat;

% roi related variables
p.roi_center = [];
p.roi_idx = [];
p.roi_border = [];
p.roi_radious = [];
p.idx = 1;

p.cha_matrix = [];
p.min_max = [];

% get channels
default_cha = {'RedChaMean', 'GreenChaMean', 'GreenChaDfof', 'GreenChaSNR', 'lc3D'};
k = 1;

for i = 1:numel(default_cha)
    
    if isfield(p.wDat, default_cha{i})
        eval(['p.cha_matrix(:, :, :, k) = p.wDat.' default_cha{i} ';']);
        eval(['p.min_max(k, :) = [min(p.wDat.' default_cha{i} ...
            '(:)) max(p.wDat.' default_cha{i} '(:))];']);
        p.chaname{k} = default_cha{i};
        k = k + 1;
    end
    
end

% set plane number
set(handles.slider1, 'Value', 1)
set(handles.slider1, 'Min', 1)
set(handles.slider1, 'Max', p.wDat.vSize(3))
numSteps = p.wDat.vSize(3);
set(handles.slider1, 'SliderStep', [1/(numSteps-1) , 1/(numSteps-1) ]);

% set channel number
set(handles.slider2, 'Value', 1)
set(handles.slider2, 'Min', 1)
set(handles.slider2, 'Max', size(p.cha_matrix, 4))
numSteps = size(p.cha_matrix, 4);
set(handles.slider2, 'SliderStep', [1/(numSteps-1) , 1/(numSteps-1) ]);

p.z = get(handles.slider1, 'Value');
p.cha = get(handles.slider2, 'Value');

% display min max
set(handles.edit1, 'String', num2str(p.min_max(p.cha, 1)))
set(handles.edit2, 'String', num2str(p.min_max(p.cha, 2)))

% plot image
plot_plane(handles)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global p
plot_plane(handles)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global p
plot_plane(handles)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
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
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global p
% update wDat

fprintf(['Updating File ', p.files2load{p.currentfile}])

load(p.files2load{p.currentfile}, 'wDat')

% pass ROI variables
wDat.manROI.roi_center = p.roi_center;
wDat.manROI.roi_idx = p.roi_idx;
wDat.manROI.roi_border = p.roi_border;
wDat.manROI.roi_radious = p.roi_radious;
wDat.ROI_center_matrix = zeros(wDat.vSize);

for i = 1:numel(p.roi_idx)
   wDat.ROI_center_matrix(round(p.roi_center(i, 2)), ...
       round(p.roi_center(i, 1)), p.roi_center(i, 3)) = 1;
end

save(p.files2load{p.currentfile}, 'wDat', '-append')

fprintf(' ... Done\n')
