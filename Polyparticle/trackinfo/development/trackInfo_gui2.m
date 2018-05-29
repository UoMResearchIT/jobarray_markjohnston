function varargout = trackInfo_gui2(varargin)
% STILL IN DEVELOPMENT!
% Author: David A. Kenwright May 2011
% For use with the results of PolyParticleTracker (developed by Salman S. Rogers)
% Displays the results of the tracking algorithm overlaid on the movie
% frames (currently TIFF only), with slider to move between frames

%TRACKINFO_GUI2 M-file for trackInfo_gui2.fig
%      TRACKINFO_GUI2, by itself, creates a new TRACKINFO_GUI2 or raises the existing
%      singleton*.
%
%      H = TRACKINFO_GUI2 returns the handle to a new TRACKINFO_GUI2 or the handle to
%      the existing singleton*.
%
%      TRACKINFO_GUI2('Property','Value',...) creates a new TRACKINFO_GUI2 using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to trackInfo_gui2_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      TRACKINFO_GUI2('CALLBACK') and TRACKINFO_GUI2('CALLBACK',hObject,...) call the
%      local function named CALLBACK in TRACKINFO_GUI2.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackInfo_gui2

% Last Modified by GUIDE v2.5 05-May-2011 11:12:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @trackInfo_gui2_OpeningFcn, ...
    'gui_OutputFcn',  @trackInfo_gui2_OutputFcn, ...
    'gui_LayoutFcn',  [], ...
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


% --- Executes just before trackInfo_gui2 is made visible.
function trackInfo_gui2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for trackInfo_gui2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trackInfo_gui2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trackInfo_gui2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function moviePathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to moviePathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of moviePathEdit as text
%        str2double(get(hObject,'String')) returns contents of moviePathEdit as a double


% --- Executes during object creation, after setting all properties.
function moviePathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to moviePathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in moviePathButton.
function moviePathButton_Callback(hObject, eventdata, handles)
textString1 = get(handles.moviePathEdit,'String');
textString2 = get(handles.prefix,'String');
images = dir([textString1 '\' textString2 '*.tif']);
set(handles.text1,'String',images(1).name)
image = imread([textString1 '\' images(1).name]);

%Set the image for contourSub
set(gcf,'CurrentAxes',handles.contourSub)
set(handles.contourSub,'NextPlot','add')
imshow(single(image-min(min(image)))/single(max(max(image))-min(min(image)))*256,gray(256));
set(handles.contourSub,'NextPlot','add')

%Set the image for displaying the tracks
set(handles.axes1,'NextPlot','add')
set(gcf,'CurrentAxes',handles.axes1)
imshow(single(image-min(min(image)))/single(max(max(image))-min(min(image)))*256,gray(256));
set(handles.axes1,'NextPlot','add')

set(handles.tracksSlider,'Min',1)
set(handles.tracksSlider,'Max',length(images));
set(handles.tracksSlider,'Sliderstep',[1/length(images) 0.2]);
set(handles.particleSlider,'Min',1)
set(handles.particleSlider,'Max',length(images));
set(handles.particleSlider,'Sliderstep',[1/length(images) 0.2]);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over moviePathButton.
function moviePathButton_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to moviePathButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function dataPathEdit_Callback(hObject, eventdata, handles)
% hObject    handle to dataPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dataPathEdit as text
%        str2double(get(hObject,'String')) returns contents of dataPathEdit as a double


% --- Executes during object creation, after setting all properties.
function dataPathEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataPathEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function tracksSlider_Callback(hObject, eventdata, handles)
textString1 = get(handles.moviePathEdit,'String');
textString2 = get(handles.prefix,'String');
images = dir([textString1 '\' textString2 '*.tif']);
frame_no = round(get(hObject,'Value')); %returns position of slider
image = imread([textString1 '\' textString2 num2str(frame_no) '.tif']);
% set(gcf,'CurrentAxes',handles.axes1)
imshow(single(image-min(min(image)))/single(max(max(image))-min(min(image)))*256,gray(256));
set(handles.text1,'String',images(frame_no).name)

% textString2 = get(handles.dataPathEdit,'String');
% load(textString2);
% hold on
% plot(handles.axes1,eval([textString3 '{i}(:,1)']),eval([textString3 '{i}(:,2)']),'r')
% for i = 1:length(eval(textString3))
% plot(handles.axes1,eval([textString3 '{i}(frame_no,1)']),eval([textString3 '{i}(frame_no,2)']),'r')
% end

% hObject    handle to tracksSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function tracksSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracksSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function dataPathButton_Callback(hObject, eventdata, handles)

global lifetime meanSpeed totalDisplacement tr tracksmooth position

textString2 = get(handles.dataPathEdit,'String');
% textString3 = get(handles.edit3,'String');
% load(textString2,textString3);
load(textString2)
%Just in case the stitched tracks haven't been created
if exist('tr','var') == 0
    tThresh=14;
    pixThresh=5;
    tr=stitchTracksTogether(tracks, tThresh, pixThresh);
end

if exist('tracksmooth','var') == 0 || exist('contour','var') == 0 || exist('position','var') == 0
    tracksmooth=cell(size(tr));
    contour=cell(size(tr));
    position=cell(size(tr));
    Lpix=4;
    for n=1:length(tr)
        [position{n},contour{n},tracksmooth{n}]=positionAveraging(tr{n},Lpix);
    end
end

[meanSpeed,totalDisplacement,lifetime,indices] = track_analysis_dak(tr);

hold on

for i = indices
    [nRuns,~,nHalts,~]=runhalt(position{i});
    eval(['track_line.n' num2str(i) '= plot(handles.axes1,  tr{i}(:,1),tr{i}(:,2));']),...
        set(eval(['track_line.n' num2str(i)]),'color','r','DisplayName',['Particle:' num2str(i)],...
        'Tag',['nRuns:' num2str(nRuns)...
        ' nHalts:' num2str(nHalts)...
        ]);
end

set(gcf,'CurrentAxes',handles.axes1)

dcm_obj = datacursormode(handles.figure1);
set(dcm_obj,'DisplayStyle','window',...%%'datatip',...
    'SnapToDataVertex','on','Enable','on')
set(dcm_obj,'UpdateFcn',@get_data)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function txt = get_data(~,event_obj)
global lifetime meanSpeed totalDisplacement tr n
handles = guihandles(trackInfo_gui2);

tar = get(event_obj,'Target');
nam = get(tar,'DisplayName');
tag = get(tar,'Tag');

n = str2double(nam(10:end));

if ~isempty(n)
set(handles.particleSlider,'Value',tr{n}(1,3))
set(handles.particleSlider,'Min',tr{n}(1,3))
set(handles.particleSlider,'Max',tr{n}(end,3));
set(handles.particleSlider,'Sliderstep',[1/(length(tr{n})-1) 1/(length(tr{n})-1)]);

set(handles.min_frame,'String',num2str(tr{n}(1,3)))
set(handles.max_frame,'String',num2str(tr{n}(end,3)))

txt = {nam,...
    tag,...
    ['Length: ' num2str(round((totalDisplacement(str2double(nam(10:end)))/11.72)*100)/100) 'microns'],...
    ['Lifetime: ' num2str(round((lifetime(str2double(nam(10:end)))/28)*100)/100) 's'],...
    ['avSpeed: ' num2str(round(meanSpeed(str2double(nam(10:end)))*(28/11.72)*100)/100) 'microns/s']...
    ...['Data1: ',''],...
    ... ['Data2: ',''
    };
end


function [meanSpeed,totalDisplacement,lifetime,indices] = track_analysis_dak(stitched_tracks)

tr = stitched_tracks;

Lpix=4;
contour=cell(size(tr));
position=cell(size(tr));
tracksmooth=cell(size(tr));

for n=1:length(tr)
    [position{n},contour{n},tracksmooth{n}]=positionAveraging(tr{n},Lpix);
end

% 1. select particles with contour length greater than totalDisplacementThresh pixels
totalDisplacementThresh=23.44; %23.44=2um, 11.72=1um or Set value,cascade

finalPosition=zeros(size(position));
totalDisplacement=zeros(size(position));
lifetime=calclifetime(tr);
for n=1:length(tr)
    finalPosition(n)=position{n}(end);
    totalDisplacement(n)=sqrt(sum(diff(tr{n}([1 end],1:2)).^2));
end
meanSpeed=finalPosition./lifetime;

temp=sortrows([totalDisplacement' (1:length(position))']);
indices=temp(temp(:,1)>totalDisplacementThresh,2)';


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1

% --- Executes during object creation, after setting all properties.
function contourSub_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on slider movement.
function particleSlider_Callback(hObject, eventdata, handles)
global tr n position
textString1 = get(handles.moviePathEdit,'String');
textString2 = get(handles.prefix,'String');
% images = dir([textString1 '\*.tif']);
frame_no = round(get(hObject,'Value')); %returns position of slider

set(handles.curr_frame,'String',num2str(frame_no))

% image = imread([textString1 '\' images(frame_no).name]);
image = imread([textString1 '\' textString2 num2str(frame_no) '.tif']);
% set(gcf,'CurrentAxes',handles.axes1)
set(handles.contourSub,'NextPlot','replace')
imshow(single(image-min(min(image)))/single(max(max(image))-min(min(image)))*256,gray(256),'Parent',handles.contourSub);
axis(handles.contourSub,[min(tr{n}(:,1))-10 max(tr{n}(:,1))+10 ...
    min(tr{n}(:,2))-10 max(tr{n}(:,2))+10])
set(handles.contourSub,'NextPlot','add')
%Plot the position for that frame
plot(handles.contourSub,tr{n}((tr{n}(:,3) == frame_no),1),tr{n}((tr{n}(:,3) == frame_no),2),'go','markersize',20);
% %Mark the start and end positions of each run
plot(handles.contourSub,tr{n}(1,1),tr{n}(1,2),'yo','markersize',10,'linewidth',2);
text(double(tr{n}(1,1)+3), double(tr{n}(1,2)+3), ...
    'Start','fontsize',14,'color','w','Parent',handles.contourSub);
plot(handles.contourSub,tr{n}(end,1),tr{n}(end,2),'y+','markersize',10,'linewidth',2);
text(double(tr{n}(end,1)+3),double(tr{n}(end,2)+3),...
    'End','fontsize',14,'color','w','Parent',handles.contourSub);

%Plot the runs and rests of the track
[nRuns,runToggle,nHalts,haltToggle]=runhalt(position{n});
for N=1:nHalts
    indices2=find(haltToggle==N);
    plot(handles.contourSub,tr{n}( tr{n}(indices2,3)<frame_no,1),tr{n}(tr{n}(indices2,3)<frame_no,2),'g','linewidth',1);
end
for N=1:nRuns
    indices2=find(runToggle==N);
    plot(handles.contourSub,tr{n}( tr{n}(indices2,3)<frame_no,1),tr{n}(tr{n}(indices2,3)<frame_no,2),'m','linewidth',1);
end

% %Set the axis of the track subplot and insert length scale marker
rectangle('Position',[min(tr{n}(:,1))-9,min(tr{n}(:,2))-9,11.72,1]...
    ,'facecolor','w','Parent',handles.contourSub)
text(double(min(tr{n}(:,1))-6),double(min(tr{n}(:,2))-6)...
    ,'1{\mu}m','fontsize',14,'color','w','Parent',handles.contourSub);

% hObject    handle to particleSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function particleSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to particleSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function prefix_Callback(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefix as text
%        str2double(get(hObject,'String')) returns contents of prefix as a double


% --- Executes during object creation, after setting all properties.
function prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
