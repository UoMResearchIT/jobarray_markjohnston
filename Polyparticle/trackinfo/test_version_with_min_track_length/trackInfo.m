function varargout = trackInfo(varargin)
% Further help information to go here
% Created by Dr D A Kenwright May 2011
% For use with PolyParticleTracker, developed by SS Rogers.
%
% This GUI is for evaluating the tracking results of PolyParticleTracker
%
% Type trackInfo to begin
%
% Movies can be either a single avi file, or in the form of tif files with
% a common
% prefix, followed by the frame number. It should be the same as used for
% generating the tracking data with PolyParticleTracker.
%
% FOR AVI FILES:
% E.g. in the movie prefix box, type:
% C:\Movies\my_movie.avi
% and then hit 'Load movie'.
% The first frame of the movie should appear below.
%
% FOR TIF FILES
% Movies should (currently) be in the form of tif files with a common
% prefix, followed by the frame number.
%
% E.g. If you have a folder 'C:\Movies\' containing tif files control_t1,
% control_t2, control_t3, control_t4...
% then in the movie prefix box, type:
% C:\Movies\control_t
% and then hit 'Load movie'.
% The first frame of the movie should appear below.
%
%
% To load the tracks, type the location of the MAT file with the output from
% PolyParticleTracker (the 'tracks' variable) and hit 'Load tracks'.
%
% Individual tracks can now be selected by clicking on the trajectories in
% the movie window. Once selected, their zoomed trajectory is plotted in
% the central region. The trajectory can be traced throguh time with the
% slider above.
%

%
% trackInfo M-file for trackInfo.fig
%      trackInfo, by itself, creates a new trackInfo or raises the existing
%      singleton*.
%
%      H = trackInfo returns the handle to a new trackInfo or the handle to
%      the existing singleton*.
%
%      trackInfo('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in trackInfo.M with the given input arguments.
%
%      trackInfo('Property','Value',...) creates a new trackInfo or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trackInfo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trackInfo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trackInfo

% Last Modified by GUIDE v2.5 15-Jun-2011 15:25:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @trackInfo_OpeningFcn, ...
    'gui_OutputFcn',  @trackInfo_OutputFcn, ...
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


% --- Executes just before trackInfo is made visible.
function trackInfo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trackInfo (see VARARGIN)

% Choose default command line output for trackInfo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trackInfo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trackInfo_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in moviebutton.
function moviebutton_Callback(hObject, eventdata, handles)
% hObject    handle to moviebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

movie = get(handles.moviePathEdit,'String');
if strcmp(movie(end-3:end),'.avi') == 1
    handles.image = read(mmreader(movie), 1);
    imshow(handles.image,'Parent',handles.movieAxes);
else
    handles.image = imread([movie '1.tif']);
    imshow(single(handles.image-min(min(handles.image)))/single(max(max(handles.image))-min(min(handles.image)))*256,gray(256),'Parent',handles.movieAxes);
end
%
set(handles.MovieName,'String',[movie '1'])
% set(handles.movieAxes,'Visible','On')

set(handles.MovieName,'Visible','On')
set(handles.movieAxes,'NextPlot','Add')
%Set some default values which can be changed with the Settings button
handles.stitchOpt.perform = 1; %Stitch tracks by default
handles.stitchOpt.tThresh = 0.5; %Time over which to stitch tracks - Set at half a second
handles.stitchOpt.pixThresh = 0.5; %Distance over which to stitch tracks - Set at half a micron
handles.track_analysisOpt.displacementThresh = 1;
handles.runhaltOpt.Lpix = 0.1; %Distance over which minimum run occurs - set at 0.1 microns
handles.runhaltOpt.haltTscale = 0.5; %Set at half a second
handles.positionAveragingOpt.Lpix = 0.5; %Set at half a micron
handles.minNoFrames = 1; %
set(handles.pushbuttonSettings,'Visible','On')
%Keep all the handle changes
guidata(hObject, handles);

function tracksedit_Callback(hObject, eventdata, handles)
% hObject    handle to tracksedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tracksedit as text
%        str2double(get(hObject,'String')) returns contents of tracksedit as a double


% --- Executes during object creation, after setting all properties.
function tracksedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tracksedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tracksbutton.
function tracksbutton_Callback(hObject, eventdata, handles)
% hObject    handle to tracksbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.plotPanel,'Visible','On')

track_string = get(handles.tracksedit,'String');
load(track_string);

    %Stitch tracks if box is ticked
if handles.stitchOpt.perform == 1
    tThresh=handles.stitchOpt.tThresh*str2double(get(handles.framerateedit,'String'));
    pixThresh=handles.stitchOpt.pixThresh*str2double(get(handles.pixelsizeedit,'String'));
    handles.tr=stitchTracksTogether(tracks, tThresh, pixThresh);
else handles.tr = tracks;
end

%Remove those than last less than min number of frames
lengthOK = [];
for j = 1:length(handles.tr)
    if length(handles.tr{j}) >= handles.minNoFrames
        lengthOK = [lengthOK j];
    end
end
tempcell = {};
for k = 1:length(lengthOK)
    tempcell{k} = handles.tr{lengthOK(k)};
end
handles.tr = tempcell;
clear tempcell;

% handles.tr = tracks;
% if exist('tracksmooth','var') == 0 || exist('contour','var') == 0 || exist('position','var') == 0
handles.tracksmooth=cell(size(handles.tr));
handles.contour=cell(size(handles.tr));
handles.position=cell(size(handles.tr));
Lpix=handles.positionAveragingOpt.Lpix*str2double(get(handles.pixelsizeedit,'String'));
for n=1:length(handles.tr)
    [handles.position{n},handles.contour{n},handles.tracksmooth{n}]=positionAveraging(handles.tr{n},Lpix);
end
% else handles.position = position;
%     handles.tracksmooth = tracksmooth;
% end
displacementThresh = handles.track_analysisOpt.displacementThresh*str2double(get(handles.pixelsizeedit,'String'));
[handles.meanSpeed,handles.totalDisplacement,handles.lifetime,handles.indices] = track_analysis_dak3(handles.tr,handles.position,displacementThresh);

Lpix = handles.runhaltOpt.Lpix*str2double(get(handles.pixelsizeedit,'String'));
haltTscale = handles.runhaltOpt.haltTscale*str2double(get(handles.framerateedit,'String'));
for i = handles.indices
    [nRuns,~,nHalts,~]=runhalt2(handles.position{i},Lpix,haltTscale);
    eval(['track_line.n' num2str(i) '= plot(handles.movieAxes,  handles.tr{i}(:,1),handles.tr{i}(:,2));']),...
        set(eval(['track_line.n' num2str(i)]),'color','r','DisplayName',['Particle:' num2str(handles.tr{i}(1,4))],...
        'Tag',['nRuns:' num2str(nRuns)...
        ' nHalts:' num2str(nHalts)...
        ]);
end
dcm_obj = datacursormode(trackInfo);
set(dcm_obj,'DisplayStyle','window',...%%'datatip',...
    'SnapToDataVertex','on','Enable','on')
set(dcm_obj,'UpdateFcn',{@get_data, handles})
guidata(hObject, handles);

% --- Executes on slider movement.
function trackSlider_Callback(hObject, eventdata, handles)
% hObject    handle to trackSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
pixel_size = str2double(get(handles.pixelsizeedit,'String'));
frame_no = round(get(hObject,'Value'));
set(handles.sliderValue,'String',['Frame:' num2str(frame_no)])

movie = get(handles.moviePathEdit,'String');
if strcmp(movie(end-3:end),'.avi') == 1
    handles.image = read(mmreader(movie), frame_no);
    set(gcf,'CurrentAxes',handles.trackAxes)
    imshow(handles.image,'Parent',handles.trackAxes);
else
    handles.image = imread([movie num2str(frame_no) '.tif']);
    imshow(single(handles.image-min(min(handles.image)))/single(max(max(handles.image))-min(min(handles.image)))*256,gray(256),'Parent',handles.trackAxes);
end
set(handles.MovieName,'String',[movie num2str(frame_no)])
set(handles.trackAxes,'NextPlot','Replace')

set(handles.trackAxes,'NextPlot','Add')
% nam = get(handles.particlenum,'String');
nam = handles.nam;
num = str2double(nam(10:end));
%Need to find which particle is numbered num
for i=1:length(handles.tr)
    if handles.tr{i}(1,4) == num
        num = i;
        break
    end
end
%Set the axis of the track subplot and insert length scale marker
axis(handles.trackAxes,[min(handles.tr{num}(:,1))-10 max(handles.tr{num}(:,1))+10 ...
    min(handles.tr{num}(:,2))-10 max(handles.tr{num}(:,2))+10])
rectangle('Position',[min(handles.tr{num}(:,1))-9,min(handles.tr{num}(:,2))-9,pixel_size,1]...
    ,'facecolor','w','Parent',handles.trackAxes)
text(double(min(handles.tr{num}(:,1))-6),double(min(handles.tr{num}(:,2))-6)...
    ,'1{\mu}m','fontsize',14,'color','w','Parent',handles.trackAxes);
% set(handles.trackAxes,'Visible','On')
%Mark the start and end positions of each run
plot(handles.trackAxes,handles.tr{num}(1,1),handles.tr{num}(1,2),'yo','markersize',10,'linewidth',2);
text(double(handles.tr{num}(1,1)+3), double(handles.tr{num}(1,2)+3), ...
    'Start','fontsize',14,'color','w','Parent',handles.trackAxes);
plot(handles.trackAxes,handles.tr{num}(end,1),handles.tr{num}(end,2),'y+','markersize',10,'linewidth',2);
text(double(handles.tr{num}(end,1)+3),double(handles.tr{num}(end,2)+3),...
    'End','fontsize',14,'color','w','Parent',handles.trackAxes);

%Plot the position for that frame
plot(handles.trackAxes,handles.tr{num}((handles.tr{num}(:,3) == frame_no),1),handles.tr{num}((handles.tr{num}(:,3) == frame_no),2),'go','markersize',20);

if get(handles.holdtracks,'Value') == 1
    %Plot the runs and rests of the track
    Lpix = handles.runhaltOpt.Lpix*str2double(get(handles.pixelsizeedit,'String'));
    haltTscale = handles.runhaltOpt.haltTscale*str2double(get(handles.framerateedit,'String'));
    [nRuns,runToggle,nHalts,haltToggle]=runhalt2(handles.position{num},Lpix,haltTscale);
    for N=1:nHalts
        indices2=find(haltToggle==N);
        plot(handles.trackAxes,handles.tr{num}(indices2,1),handles.tr{num}(indices2,2),'g','linewidth',2);
    end
    for N=1:nRuns
        indices2=find(runToggle==N);
        plot(handles.trackAxes,handles.tr{num}(indices2,1),handles.tr{num}(indices2,2),'m','linewidth',2);
    end
end

%If the box is checked, plot the distance travelled
if get(handles.distspeedcheck,'Value') == 1
    dot = (frame_no-get(handles.trackSlider,'Min'))+1;
    plotDistSpeed_runhalt(handles,nam,dot,round(get(handles.windowSlider,'Value')));
    set(handles.speedSave,'Visible','On')
    set(handles.speedSaveTxt,'Visible','On')
else cla(handles.distancePlot)
    set(handles.distancePlot,'Visible','Off')
    cla(handles.speedPlot)
    set(handles.speedPlot,'Visible','Off')
    set(handles.speedSave,'Visible','Off')
    set(handles.speedSaveTxt,'Visible','Off')
    set(handles.windowSlider,'Visible','Off')
end



% --- Executes during object creation, after setting all properties.
function trackSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trackSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on button press in holdtracks.
function holdtracks_Callback(hObject, eventdata, handles)
% hObject    handle to holdtracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of holdtracks


%%%
%----------------------------------------------------------------------
%GET_DATA FUNCTION
%This is what runs when you click on a selected track
function txt = get_data(~,event_obj, handles) %txt = get_data(~,event_obj, hObject)
% handles = guidata(trackInfo);
pixel_size = str2double(get(handles.pixelsizeedit,'String'));
frame_rate = str2double(get(handles.framerateedit,'String'));
tar = get(event_obj,'Target');
handles.nam = get(tar,'DisplayName');
nam = handles.nam;
tag = get(tar,'Tag');
num = str2double(nam(10:end));
%Need to find which particle is numbered num
for i=1:length(handles.tr)
    if handles.tr{i}(1,4) == num
        num = i;
        break
    end
end
txt = {nam,...
    tag,...
    ['Length: ' num2str(round((handles.totalDisplacement(num)/pixel_size)*100)/100) 'microns'],...
    ['Lifetime: ' num2str(round((handles.lifetime(num)/frame_rate)*100)/100) 's'],...
    ['avSpeed: ' num2str(round(handles.meanSpeed(num)*(frame_rate/pixel_size)*100)/100) 'microns/s']};

% Make the track appear bold and magenta, reset all other lines
set(findobj('Type','line'),'LineWidth',0.5,'Color','r')
set(gco,'LineWidth',2,'Color','m')

%Now for the individual track...
set(handles.trackAxes,'NextPlot','Replace')
set(handles.trackAxes,'Visible','On')
movie = get(handles.moviePathEdit,'String');

if strcmp(movie(end-3:end),'.avi') == 1
    handles.image = read(mmreader(movie), 1);
    set(handles.figure1,'CurrentAxes',handles.trackAxes)
    imshow(handles.image,'Parent',handles.trackAxes);
else
    handles.image = imread([movie '1.tif']);
    imshow(single(handles.image-min(min(handles.image)))/single(max(max(handles.image))-min(min(handles.image)))*256,gray(256),'Parent',handles.trackAxes);
end

set(handles.trackAxes,'NextPlot','Add')
%Set the axis of the track subplot and insert length scale marker
axis(handles.trackAxes,[min(handles.tr{num}(:,1))-10 max(handles.tr{num}(:,1))+10 ...
    min(handles.tr{num}(:,2))-10 max(handles.tr{num}(:,2))+10])
rectangle('Position',[min(handles.tr{num}(:,1))-9,min(handles.tr{num}(:,2))-9,pixel_size,1]...
    ,'facecolor','w','Parent',handles.trackAxes)
text(double(min(handles.tr{num}(:,1))-6),double(min(handles.tr{num}(:,2))-6)...
    ,'1{\mu}m','fontsize',14,'color','w','Parent',handles.trackAxes);


%Mark the start and end positions of each run
plot(handles.trackAxes,handles.tr{num}(1,1),handles.tr{num}(1,2),'yo','markersize',10,'linewidth',2);
text(double(handles.tr{num}(1,1)+3), double(handles.tr{num}(1,2)+3), ...
    'Start','fontsize',14,'color','w','Parent',handles.trackAxes);
plot(handles.trackAxes,handles.tr{num}(end,1),handles.tr{num}(end,2),'y+','markersize',10,'linewidth',2);
text(double(handles.tr{num}(end,1)+3),double(handles.tr{num}(end,2)+3),...
    'End','fontsize',14,'color','w','Parent',handles.trackAxes);

Lpix = handles.runhaltOpt.Lpix*str2double(get(handles.pixelsizeedit,'String'));
haltTscale = handles.runhaltOpt.haltTscale*str2double(get(handles.framerateedit,'String'));
[nRuns,runToggle,nHalts,haltToggle]=runhalt2(handles.position{num},Lpix,haltTscale);
for N=1:nHalts
    indices2=find(haltToggle==N);
    plot(handles.trackAxes,handles.tr{num}(indices2,1),handles.tr{num}(indices2,2),'g','linewidth',2);
end
for N=1:nRuns
    indices2=find(runToggle==N);
    plot(handles.trackAxes,handles.tr{num}(indices2,1),handles.tr{num}(indices2,2),'m','linewidth',2);
end

%If the box is checked, plot the distance travelled
if get(handles.distspeedcheck,'Value') == 1
    plotDistSpeed_runhalt(handles,nam,1,1);
    set(handles.speedSave,'Visible','On')
    set(handles.speedSaveTxt,'Visible','On')
else cla(handles.distancePlot)
    set(handles.distancePlot,'Visible','Off')
    cla(handles.speedPlot)
    set(handles.speedPlot,'Visible','Off')
    set(handles.windowSlider,'Visible','Off')
    set(handles.windowSize,'Visible','Off')
    set(handles.speedSave,'Visible','Off')
    set(handles.speedSaveTxt,'Visible','Off')
end

%Reset the slider values before changing
set(handles.trackSlider,'Min',1);
set(handles.trackSlider,'Value',1)
set(handles.trackSlider,'Max',99999999999);
%Set the slider values
set(handles.trackSlider,'Value',handles.tr{num}(1,3))
set(handles.trackSlider,'Min',handles.tr{num}(1,3));
set(handles.trackSlider,'Max',handles.tr{num}(end,3));
set(handles.trackSlider,'Sliderstep',[1/(handles.tr{num}(end,3) - handles.tr{num}(1,3)+1) 10/(handles.tr{num}(end,3) - handles.tr{num}(1,3)+1)]);
set(handles.trackSlider,'Visible','On')
%Set the values to be displayed
set(handles.sliderMin,'String',['Min:' num2str(handles.tr{num}(1,3))])
set(handles.sliderMax,'String',['Max:' num2str(handles.tr{num}(end,3))])
set(handles.sliderValue,'String',['Frame:' num2str(get(handles.trackSlider,'Value'))])
set(handles.sliderMin,'Visible','On')
set(handles.sliderMax,'Visible','On')
set(handles.sliderValue,'Visible','On')

%Reset the windowing slider
set(handles.windowSlider,'Min',1);
set(handles.windowSlider,'Value',1)
set(handles.windowSlider,'Max',handles.tr{num}(end,3)-handles.tr{num}(1,3)+1);
set(handles.windowSlider,'Sliderstep',[1/(get(handles.windowSlider,'Max')) 10/(get(handles.windowSlider,'Max'))])


 guidata(handles.figure1, handles);
% set(handles.particlenum,'String',nam)
% set(handles.holdtracks,'Visible','On')
% set(handles.plotPanel,'Visible','On')
% set(handles.distspeedcheck,'Visible','On')




function framerateedit_Callback(hObject, eventdata, handles)
% hObject    handle to framerateedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framerateedit as text
%        str2double(get(hObject,'String')) returns contents of framerateedit as a double


% --- Executes during object creation, after setting all properties.
function framerateedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framerateedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pixelsizeedit_Callback(hObject, eventdata, handles)
% hObject    handle to pixelsizeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pixelsizeedit as text
%        str2double(get(hObject,'String')) returns contents of pixelsizeedit as a double


% --- Executes during object creation, after setting all properties.
function pixelsizeedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pixelsizeedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in distspeedcheck.
function distspeedcheck_Callback(hObject, eventdata, handles)
% hObject    handle to distspeedcheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of distspeedcheck


% --- Executes on slider movement.
function windowSlider_Callback(hObject, eventdata, handles)
% hObject    handle to windowSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
nam = handles.nam;
frame_no = round(get(handles.trackSlider,'Value'));
dot = (frame_no-get(handles.trackSlider,'Min'))+1;
plotDistSpeed_runhalt(handles,nam,dot,round(get(hObject,'Value')));
set(handles.windowSize,'String',['Window size: ' num2str(round(get(hObject,'Value'))) ' frame(s)'])

% --- Executes during object creation, after setting all properties.
function windowSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in helpButton.
function helpButton_Callback(hObject, eventdata, handles)
% hObject    handle to helpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
doc trackInfo


% --- Executes on button press in angleButton.
function angleButton_Callback(hObject, eventdata, handles)
% hObject    handle to angleButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.angleSave,'Visible','Off')
set(handles.angleSaveTxt,'Visible','Off')
set(handles.figure1,'CurrentAxes',handles.anglePlot)
movie = get(handles.moviePathEdit,'String');

if strcmp(movie(end-3:end),'.avi') == 1
    handles.image = read(mmreader(movie), 1);
    imshow(handles.image,'Parent',handles.anglePlot);
else
    handles.image = imread([movie '1.tif']);
    imshow(single(handles.image-min(min(handles.image)))/single(max(max(handles.image))-min(min(handles.image)))*256,gray(256),'Parent',handles.anglePlot);
end
% set(handles.anglePlot,'Visible','On')

tr2 = cell(length(handles.indices),1);
for i =1:length(handles.indices)
    tr2{i} = handles.tr{handles.indices(i)};
end
[direction] = particle_angle2centre(tr2,handles.anglePlot);
handles.direction = direction;
set(handles.angleSave,'Visible','On')
set(handles.angleSaveTxt,'Visible','On')
set(handles.directionPanel,'Visible','On')
set(handles.inwardNo,'Visible','On')
set(handles.inwardNo,'String',['Inward: ' num2str(direction.no_forward)])
set(handles.outwardNo,'Visible','On')
set(handles.outwardNo,'String',['Outward: ' num2str(direction.no_backward)])
set(handles.transverseNo,'Visible','On')
set(handles.transverseNo,'String',['Transverse: ' num2str(direction.no_transverse)])
set(handles.figure1,'CurrentAxes',handles.movieAxes)
dcm_obj = datacursormode(trackInfo);
set(dcm_obj,'DisplayStyle','window',...%%'datatip',...
    'SnapToDataVertex','on','Enable','on')
set(dcm_obj,'UpdateFcn',{@get_data, handles})
guidata(hObject, handles);


% --- Executes on button press in angleSave.
function angleSave_Callback(hObject, eventdata, handles)
% hObject    handle to angleSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
assignin('base','direction',handles.direction)


% --- Executes on button press in speedSave.
function speedSave_Callback(hObject, eventdata, handles)
% hObject    handle to speedSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pixel_size = str2double(get(handles.pixelsizeedit,'String'));
frame_rate = str2double(get(handles.framerateedit,'String'));
window = round(get(handles.windowSlider,'Value'));
speed = cell(length(handles.indices),1);

%Error if there is a particle that is shorter than the window length
for i = 1:length(handles.indices)
    if length(handles.tracksmooth{handles.indices(i)}) < window
        error('trackInfo:windowChk', 'A track exists shorter than the window length. In order to export speeds with constant window length, lower the window value.')
    end
end
%

for i = 1:length(handles.indices)
    % Speed and distance plots
    %Calculate distance and speed
    initial_x = handles.tracksmooth{handles.indices(i)}(1,1);
    initial_y = handles.tracksmooth{handles.indices(i)}(1,2);
    % % % Reset positions
    tr1 = [];
    tr1(:,1) = handles.tracksmooth{handles.indices(i)}(:,1) - initial_x;
    tr1(:,2) = handles.tracksmooth{handles.indices(i)}(:,2) - initial_y;
    %Interpolate stitches
    [tr1] = interpolate_stitches(tr1);
    %Rescale to microns
    tr1(:,1:2) = tr1(:,1:2)/pixel_size;
    % calc separation
    displacement=sqrt(diff(tr1(:,2)).^2 + diff(tr1(:,1)).^2);
    % distance
    tr1(:,3) = [0; cumsum(displacement)];
    % speed - rescaled to seconds
    diffr = diff(tr1(:,3))*frame_rate;
    tr1(:,4) = [0; diffr];
    
    if window <= 2
        speed{i} = tr1(:,4);
    else speed{i} = movingslope(tr1(:,3),window);
    end
end
assignin('base','speed',speed)


function plotDistSpeed_runhalt(handles,nam,dot,window)
%Dot - the location of a dot on the plot showing current time of frame
%Window for 'wavelet' - number of points over which to calculate speed

pixel_size = str2double(get(handles.pixelsizeedit,'String'));
frame_rate = str2double(get(handles.framerateedit,'String'));

num = str2double(nam(10:end));
%Need to find which particle is numbered num
for i=1:length(handles.tr)
    if handles.tr{i}(1,4) == num
        num = i;
        break
    end
end

% % Speed and distance plots
% %Calculate distance and speed
% initial_x = handles.tracksmooth{num}(1,1);
% initial_y = handles.tracksmooth{num}(1,2);
% % % % Reset positions
% tr1 = [];
% tr1(:,1) = handles.tracksmooth{num}(:,1) - initial_x;
% tr1(:,2) = handles.tracksmooth{num}(:,2) - initial_y;
% %Interpolate stitches
% [tr1] = interpolate_stitches(tr1);
% %Rescale to microns
% tr1(:,1:2) = tr1(:,1:2)/pixel_size;
% % calc separation
% displacement=sqrt(diff(tr1(:,2)).^2 + diff(tr1(:,1)).^2);
% % distance
% tr1(:,3) = [0; cumsum(displacement)];
% % speed - rescaled to seconds
% diffr = diff(tr1(:,3))*frame_rate;
% tr1(:,4) = [0; diffr];



Lpix = handles.runhaltOpt.Lpix*str2double(get(handles.pixelsizeedit,'String'));
haltTscale = handles.runhaltOpt.haltTscale*str2double(get(handles.framerateedit,'String'));
[nRuns,runToggle,nHalts,haltToggle]=runhalt2(handles.position{num},Lpix,haltTscale);

set(handles.distancePlot,'NextPlot','Replace')
% time = 1/frame_rate:1/frame_rate:length(tr1(:,3))/frame_rate;
time = 1/frame_rate:1/frame_rate:length(handles.position{num})/frame_rate;

% for i = 1:nRuns
%     plot(handles.distancePlot,time(find(runToggle==i)),tr1(find(runToggle==i),3),'m','LineWidth',1)
%     set(handles.distancePlot,'NextPlot','Add')
% end
for i = 1:nRuns
    plot(handles.distancePlot,time(find(runToggle==i)),handles.position{num}(find(runToggle==i))/pixel_size,'m','LineWidth',1)
    set(handles.distancePlot,'NextPlot','Add')
end

% for i = 1:nHalts
%     plot(handles.distancePlot,time(find(haltToggle==i)),tr1(find(haltToggle==i),3),'g','LineWidth',1)
% end
for i = 1:nHalts
    plot(handles.distancePlot,time(find(haltToggle==i)),handles.position{num}(find(haltToggle==i))/pixel_size,'g','LineWidth',1)
end

% plot(handles.distancePlot,time(dot),tr1(dot,3),'rx','MarkerSize',10,'LineWidth',2)
plot(handles.distancePlot,time(dot),handles.position{num}(dot)/pixel_size,'rx','MarkerSize',10,'LineWidth',2)
ylabel(handles.distancePlot,{'Distance along'; 'track [{\mu}m]'})
xlabel(handles.distancePlot,'Time [s]')
axis(handles.distancePlot,'tight')

set(handles.speedPlot,'NextPlot','Replace')

if window <= 2
    speed = diff(handles.position{num});
    speed = [speed(1); speed];
    for i = 1:nRuns
        plot(handles.speedPlot,time(find(runToggle==i)),speed(find(runToggle==i))*(frame_rate/pixel_size),'m','LineWidth',1)
        set(handles.speedPlot,'NextPlot','Add')
    end
    for i = 1:nHalts
        plot(handles.speedPlot,time(find(haltToggle==i)),speed(find(haltToggle==i))*(frame_rate/pixel_size),'g','LineWidth',1)
    end
    plot(handles.speedPlot,time(dot),speed(dot)*(frame_rate/pixel_size),'rx','MarkerSize',10,'LineWidth',2)
    
else speed = movingslope(handles.position{num},window);
    for i = 1:nRuns
        plot(handles.speedPlot,time(find(runToggle==i)),speed(find(runToggle==i))*(frame_rate/pixel_size),'m','LineWidth',1)
        set(handles.speedPlot,'NextPlot','Add')
    end
    for i = 1:nHalts
        plot(handles.speedPlot,time(find(haltToggle==i)),speed(find(haltToggle==i))*(frame_rate/pixel_size),'g','LineWidth',1)
    end
    plot(handles.speedPlot,time(dot),speed(dot)*(frame_rate/pixel_size),'rx','MarkerSize',10,'LineWidth',2)
end
ylabel(handles.speedPlot,{'Speed'; '[{\mu}ms^{-1}]'})
xlabel(handles.speedPlot,'Time [s]')
axis(handles.speedPlot,'tight')
set(handles.speedPlot,'NextPlot','Add')


set(handles.windowSize,'Visible','On')
set(handles.windowSlider,'Visible','On')

function [direction] = particle_angle2centre(tracks,axis_handle)

% clicktext = text(100,100,'Click on centre','FontSize',14,'Color','w');
temp_txt=xlabel('Click on centre','FontSize',10);
%User then clicks the location of the 'centre' point of reference
[centre_x centre_y] = ginput(1);
delete(temp_txt)

angle=zeros(length(tracks),1);
for i = 1:length(tracks)
    
    Dx = tracks{i}(end,1) - tracks{i}(1,1); %The change between first and last x coordinate (for general direction)
    Dy = tracks{i}(end,2) - tracks{i}(1,2); %The change between first and last y coordinate (for general direction)
    
    Dxc = centre_x - tracks{i}(1,1);    %The distance to the centre point
    Dyc = centre_y - tracks{i}(1,2);
    
    %The angle to the centre from the first tracking point is:
    if Dxc >= 0 && Dyc >= 0
        angle(i) = atan(Dyc/Dxc);
    elseif Dxc >= 0 && Dyc < 0
        angle(i) = (2*pi)-abs(atan(Dyc/Dxc));
    elseif Dxc < 0 && Dyc >= 0
        angle(i) = pi-abs(atan(Dyc/Dxc));
    elseif Dxc < 0 && Dyc < 0
        angle(i) = pi+atan(Dyc/Dxc);
    end
    
    %The subtract this angle from the angle that the particle travels
    if Dx >= 0 && Dy >= 0
        angle(i) = mod(atan(Dy/Dx)-angle(i),2*pi);
    elseif Dx >= 0 && Dy < 0
        angle(i) =  mod(- (2*pi)-abs(atan(Dy/Dx))-angle(i),2*pi);
    elseif Dx < 0 && Dy >= 0
        angle(i) =  mod(- pi-abs(atan(Dy/Dx))-angle(i),2*pi);
    elseif Dx < 0 && Dy < 0
        angle(i) =  mod(- pi+atan(Dy/Dx)-angle(i),2*pi);
    end
    
end

if exist('axis_handle','var') == 1
    axes(axis_handle)
end
h=rose(angle);
x = get(h,'Xdata');
y = get(h,'Ydata');

patch(x,y,'r');
ylabel({'Angle of displacement'; '[degrees]'})

% set(findobj('Type','line'),'LineWidth',1.5)
% set(findobj('Type','line'),'Color','k')
% set(findobj('Type','font'),'FontWeight','normal','FontSize',12)


%For the output:
direction.angle = angle;
direction.no_forward = 0;
direction.no_backward = 0;
direction.no_transverse = 0;


for i = 1:length(angle)
    if angle(i) < pi/4 || angle(i) >= 7*pi/4
        direction.no_forward = direction.no_forward+1;
    elseif angle(i) >= 3*pi/4 && angle(i) < 5*pi/4
        direction.no_backward = direction.no_backward+1;
        %transverse movement
    elseif (angle(i) >= pi/4 && angle(i) < 3*pi/4) || (angle(i) >= 3*pi/2 && angle(i) < 7*pi/4) || (angle(i) >= 5*pi/4 && angle(i) < 3*pi/2)
        direction.no_transverse = direction.no_transverse+1;
    end
end


% --- Executes on button press in pushbuttonSettings.
function pushbuttonSettings_Callback(hObject, eventdata, handles)
%Call the trackInfo_settings GUI
trackInfo_settings();
