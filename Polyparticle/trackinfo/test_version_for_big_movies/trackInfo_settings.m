function varargout = trackInfo_settings(varargin)
% TRACKINFO_SETTINGS M-file for trackInfo_settings.fig
%      TRACKINFO_SETTINGS, by itself, creates a new TRACKINFO_SETTINGS or raises the existing
%      singleton*.
%
%      H = TRACKINFO_SETTINGS returns the handle to a new TRACKINFO_SETTINGS or the handle to
%      the existing singleton*.
%
%      TRACKINFO_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKINFO_SETTINGS.M with the given input arguments.
%
%      TRACKINFO_SETTINGS('Property','Value',...) creates a new TRACKINFO_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackInfo_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackInfo_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackInfo_settings

% Last Modified by GUIDE v2.5 17-Aug-2011 10:21:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trackInfo_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @trackInfo_settings_OutputFcn, ...
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


% --- Executes just before trackInfo_settings is made visible.
function trackInfo_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trackInfo_settings (see VARARGIN)

% Choose default command line output for trackInfo_settings
handles.output = hObject;

% get the main_gui handle (access to the gui)
mainGUIhandle = trackInfo;       
% get the data from the gui (all handles inside gui_main)
mainGUIdata  = guidata(mainGUIhandle);

set(handles.checkbox_Stitch,'Value',mainGUIdata.stitchOpt.perform)
set(handles.edit_stitch_tThresh,'String',num2str(mainGUIdata.stitchOpt.tThresh))
set(handles.edit_stitchPixThresh,'String',num2str(mainGUIdata.stitchOpt.pixThresh))
set(handles.edit_dispThresh,'String',num2str(mainGUIdata.track_analysisOpt.displacementThresh))
set(handles.edit_runhalt_lpix,'String',num2str(mainGUIdata.runhaltOpt.Lpix))
set(handles.edit_runhalt_tscale,'String',num2str(mainGUIdata.runhaltOpt.haltTscale))
set(handles.edit_contour_lpix,'String',num2str(mainGUIdata.positionAveragingOpt.Lpix))
set(handles.editMinNoFrames,'String',num2str(mainGUIdata.minNoFrames))
set(handles.editFrameJump,'String',num2str(mainGUIdata.FrameJump))
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trackInfo_settings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trackInfo_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_Stitch.
function checkbox_Stitch_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Stitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_Stitch



function edit_stitch_tThresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stitch_tThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stitch_tThresh as text
%        str2double(get(hObject,'String')) returns contents of edit_stitch_tThresh as a double


% --- Executes during object creation, after setting all properties.
function edit_stitch_tThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stitch_tThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stitchPixThresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stitchPixThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%gets the contents of the edit text field
%and then converts it into a numeric value
%if not a number then input will be empty
input = str2num(get(hObject,'String'));
 
%checks to see if input is empty. if so, default input1_editText to zero
if (isempty(input))
     set(hObject,'String','0')
end
guidata(hObject, handles);
% Hints: get(hObject,'String') returns contents of edit_stitchPixThresh as text
%        str2double(get(hObject,'String')) returns contents of edit_stitchPixThresh as a double


% --- Executes during object creation, after setting all properties.
function edit_stitchPixThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stitchPixThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dispThresh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dispThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dispThresh as text
%        str2double(get(hObject,'String')) returns contents of edit_dispThresh as a double


% --- Executes during object creation, after setting all properties.
function edit_dispThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dispThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_runhalt_lpix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_runhalt_lpix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_runhalt_lpix as text
%        str2double(get(hObject,'String')) returns contents of edit_runhalt_lpix as a double


% --- Executes during object creation, after setting all properties.
function edit_runhalt_lpix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_runhalt_lpix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_runhalt_tscale_Callback(hObject, eventdata, handles)
% hObject    handle to edit_runhalt_tscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_runhalt_tscale as text
%        str2double(get(hObject,'String')) returns contents of edit_runhalt_tscale as a double


% --- Executes during object creation, after setting all properties.
function edit_runhalt_tscale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_runhalt_tscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_contour_lpix_Callback(hObject, eventdata, handles)
% hObject    handle to edit_contour_lpix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_contour_lpix as text
%        str2double(get(hObject,'String')) returns contents of edit_contour_lpix as a double


% --- Executes during object creation, after setting all properties.
function edit_contour_lpix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_contour_lpix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_saveSettings.
function pushbutton_saveSettings_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_saveSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get the main_gui handle (access to the gui)
mainGUIhandle = trackInfo;       
% get the data from the gui (all handles inside gui_main)
mainGUIdata  = guidata(mainGUIhandle);
 
% change main gui strings
mainGUIdata.stitchOpt.perform = get(handles.checkbox_Stitch, 'Value');
mainGUIdata.stitchOpt.tThresh = str2double(get(handles.edit_stitch_tThresh, 'String'));
mainGUIdata.stitchOpt.pixThresh = str2double(get(handles.edit_stitchPixThresh, 'String'));
mainGUIdata.track_analysisOpt.displacementThresh = str2double(get(handles.edit_dispThresh, 'String'));
mainGUIdata.runhaltOpt.Lpix = str2double(get(handles.edit_runhalt_lpix, 'String'));
mainGUIdata.runhaltOpt.haltTscale = str2double(get(handles.edit_runhalt_tscale, 'String'));
mainGUIdata.positionAveragingOpt.Lpix = str2double(get(handles.edit_contour_lpix,'String'));
mainGUIdata.minNoFrames = str2double(get(handles.editMinNoFrames,'String'));
mainGUIdata.FrameJump = str2double(get(handles.editFrameJump,'String'));
% save changed data back into main_gui
%this line updates the data of the Main Gui
guidata(mainGUIhandle, mainGUIdata);
 
% close this gui 
close(trackInfo_settings);



function editMinNoFrames_Callback(hObject, eventdata, handles)
% hObject    handle to editMinNoFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinNoFrames as text
%        str2double(get(hObject,'String')) returns contents of editMinNoFrames as a double


% --- Executes during object creation, after setting all properties.
function editMinNoFrames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinNoFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editFrameJump_Callback(hObject, eventdata, handles)
% hObject    handle to editFrameJump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrameJump as text
%        str2double(get(hObject,'String')) returns contents of editFrameJump as a double


% --- Executes during object creation, after setting all properties.
function editFrameJump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrameJump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
