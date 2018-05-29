function varargout = singleFPP(varargin)
% SINGLEFPP MATLAB code for singleFPP.fig
%      SINGLEFPP, by itself, creates a new SINGLEFPP or raises the existing
%      singleton*.
%
%      H = SINGLEFPP returns the handle to a new SINGLEFPP or the handle to
%      the existing singleton*.
%
%      SINGLEFPP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SINGLEFPP.M with the given input arguments.
%
%      SINGLEFPP('Property','Value',...) creates a new SINGLEFPP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before singleFPP_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to singleFPP_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help singleFPP

% Last Modified by GUIDE v2.5 16-May-2012 13:40:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @singleFPP_OpeningFcn, ...
    'gui_OutputFcn',  @singleFPP_OutputFcn, ...
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


% --- Executes just before singleFPP is made visible.
function singleFPP_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to singleFPP (see VARARGIN)

% Choose default command line output for singleFPP
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes singleFPP wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = singleFPP_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editFrameRate_Callback(hObject, eventdata, handles)
% hObject    handle to editFrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFrameRate as text
%        str2double(get(hObject,'String')) returns contents of editFrameRate as a double


% --- Executes during object creation, after setting all properties.
function editFrameRate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFrameRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editPixScale_Callback(hObject, eventdata, handles)
% hObject    handle to editPixScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editPixScale as text
%        str2double(get(hObject,'String')) returns contents of editPixScale as a double


% --- Executes during object creation, after setting all properties.
function editPixScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPixScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMinTrackLength_Callback(hObject, eventdata, handles)
% hObject    handle to editMinTrackLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMinTrackLength as text
%        str2double(get(hObject,'String')) returns contents of editMinTrackLength as a double
if isfield(handles,'tracks')
    plotTracks(hObject,handles)
end

% --- Executes during object creation, after setting all properties.
function editMinTrackLength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMinTrackLength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editMovieLoc_Callback(hObject, eventdata, handles)
% hObject    handle to editMovieLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMovieLoc as text
%        str2double(get(hObject,'String')) returns contents of editMovieLoc as a double
handles.MovieLoc = get(handles.editMovieLoc,'String');
guidata(hObject, handles);
loadMovie(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function editMovieLoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMovieLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonBrowseMovie.
function pushbuttonBrowseMovie_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonBrowseMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[MovieFilename, MoviePathname] = uigetfile( ...
    {'*.tif*',  'Tiff files (*.tif*)'; '*.avi*', 'Avi files (*.avi*)'}, ...
    'Select all frames (tif) or avi file','MultiSelect','On');
handles.MoviePath = MoviePathname;
handles.MovieFile = MovieFilename;
guidata(hObject, handles);
loadMovie(hObject, eventdata, handles)


function loadMovie(hObject, eventdata, handles)
%If multiple tif files were selected, store the frames and plot frame 1
if iscell(handles.MovieFile)
    handles.frames = handles.MovieFile;
    handles.image = imread([handles.MoviePath handles.MovieFile{1}]);
    imshow(single(handles.image-min(min(handles.image)))/single(max(max(handles.image))-min(min(handles.image)))*256,gray(256),'Parent',handles.axesMovie);
    set(handles.slider1,'Max',length(handles.frames),'Sliderstep',[1/length(handles.frames) 10/length(handles.frames)]);
    %If an avi movie was 
elseif strcmp(handles.MovieFile(end-3:end),'.avi') == 1
    handles.image = read(mmreader([handles.MoviePath handles.MovieFile]), 1);
    imshow(handles.image,'Parent',handles.axesMovie);
else
    handles.image = imread([handles.MoviePath handles.MovieFile]);
    imshow(single(handles.image-min(min(handles.image)))/single(max(max(handles.image))-min(min(handles.image)))*256,gray(256),'Parent',handles.axesMovie);
end
set(handles.axesMovie,'NextPlot','Add')
guidata(hObject, handles);



function editDataLoc_Callback(hObject, eventdata, handles)
% hObject    handle to editDataLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDataLoc as text
%        str2double(get(hObject,'String')) returns contents of editDataLoc as a double



% --- Executes during object creation, after setting all properties.
function editDataLoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDataLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonBrowseData.
function pushbuttonBrowseData_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonBrowseData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
temp = uiimport;
temp2 = fieldnames(temp);
handles.tracks = eval(['temp.' temp2{1}]);
plotTracks(hObject,handles)
guidata(hObject, handles)


function plotTracks(hObject, handles)
set(handles.axesMovie,'NextPlot','Replace')
imshow(single(handles.image-min(min(handles.image)))/single(max(max(handles.image))-min(min(handles.image)))*256,gray(256),'Parent',handles.axesMovie);
set(handles.axesMovie,'NextPlot','Add')
order = length(num2str(length(handles.tracks)));
for i = 1:length(handles.tracks)
    if rem(i,10^(order-2)) == 0
        h=waitbar(i/length(handles.tracks));
    end
    if ((handles.tracks{i}(end,1)-handles.tracks{i}(1,1))^2 + (handles.tracks{i}(end,2)-handles.tracks{i}(1,2))^2) >= (str2double(get(handles.editMinTrackLength,'String'))*str2double(get(handles.editPixScale,'String')))^2
        eval(['track_line.n' num2str(i) '= plot(handles.axesMovie,  handles.tracks{i}(:,1),handles.tracks{i}(:,2));']),...
            set(eval(['track_line.n' num2str(i)]),'color','r','DisplayName',['Particle:' num2str(handles.tracks{i}(1,4))]);
    end
end
close(h)

dcm_obj = datacursormode(singleFPP);
set(dcm_obj,'DisplayStyle','datatip',...%%'window',...
    'SnapToDataVertex','on','Enable','on')
set(dcm_obj,'UpdateFcn',{@get_data, handles})
guidata(hObject, handles);






%GET_DATA FUNCTION
%This is what runs when you click on a selected track
function txt = get_data(~,event_obj, handles) %txt = get_data(~,event_obj, hObject)
% handles = guidata(trackInfo);
pixel_size = 1/str2double(get(handles.editPixScale,'String'));
frame_rate = str2double(get(handles.editFrameRate,'String'));
tar = get(event_obj,'Target');
handles.nam = get(tar,'DisplayName');
nam = handles.nam;
tag = get(tar,'Tag');
num = str2double(nam(10:end));

%Need to find which particle is numbered num
for i=1:length(handles.tracks)
    if handles.tracks{i}(1,4) == num
        num = i;
        break
    end
end

% check for NaN values
if any(isnan(handles.tracks{num}(:,1)))
    %Interpolate between the values
    x = handles.tracks{num}(:,1);
    y = handles.tracks{num}(:,2);
    x(isnan(x)) = interp1(find(~isnan(x)), x(~isnan(x)), find(isnan(x)),'linear');   %Could change this to cubic...
    y(isnan(y)) = interp1(find(~isnan(y)), y(~isnan(y)), find(isnan(y)),'linear');
    handles.tracks{num}(:,1) = x;
    handles.tracks{num}(:,2) = y;
end

tracklength = sqrt((handles.tracks{num}(end,1) - handles.tracks{num}(1,1)).^2 + (handles.tracks{num}(end,2) - handles.tracks{num}(1,2)).^2)*pixel_size;
lifetime = (handles.tracks{num}(end,3)-handles.tracks{num}(1,3))/frame_rate;
avSpeed = tracklength/lifetime;

txt = {nam,...
    ['Length: ' num2str(tracklength,3) 'microns'],...
    ['Lifetime: ' num2str(lifetime,3) 's'],...
    ['avSpeed: ' num2str(avSpeed,3) 'microns/s']};

% Make the track appear bold and magenta, reset all other lines
set(findobj('Type','line'),'LineWidth',0.5,'Color','r')
set(gco,'LineWidth',2,'Color','g')

% Need to set the slider to only allow frames that the particle exists in
set(handles.slider1,'Min',1,...
    'Max',length(handles.tracks{num}(:,3)),...
    'Value',1,...
    'sliderstep',[1/(length(handles.tracks{num}(:,3))-1) 10/(length(handles.tracks{num}(:,3))-1)])

updatePlots(handles)
% guidata(hObject,handles)



function [S,FPP]=calcFPP(tr_raw,framerate,pixscale,Lscale)

tracktoggle=true(size(tr_raw));

% calc SD for each particle at timescale tau

NEvents=0;
end2end=zeros(size(tr_raw));

for n=find(tracktoggle)
    NEvents=NEvents+length(tr_raw{n}(:,1));
    end2end(n)=sqrt(sum(diff(tr_raw{n}([1 end],1:2)).^2));
end
S=cell(size(Lscale));
FPP=cell(size(Lscale));
for L=1:length(Lscale)
    timeWithin=cell(length(tr_raw),1);
    %     Lscale1=(Lscale(L)/pixscale);
    Lscale2=(Lscale(L)/pixscale)^2;
    ttoggleLhit=false(length(tr_raw),1);
    for n=find(tracktoggle)
        timeWithin{n}=NaN(size(tr_raw{n}(:,1)));
        for t=1:size(tr_raw{n},1)
            dr2end=(tr_raw{n}(end,1)-tr_raw{n}(t,1)).^2+(tr_raw{n}(end,2)-tr_raw{n}(t,2)).^2;
            if dr2end>Lscale2 %if particle has made passage at least once
                %                 %find time at which particle leaves Lscale diamond
                %                 absdrRoughDia=abs(tr_raw{n}(t:end,1)-tr_raw{n}(t,1))+abs(tr_raw{n}(t:end,2)-tr_raw{n}(t,2));
                %                 indices=cumsum(~(absdrRoughDia<Lscale1));
                %                 timeWithinDia=sum(~indices)-1;
                %
                %                 %find time at which particle leaves Lscale square
                %                 absdrRoughSq=max(abs(tr_raw{n}(t:end,1)-tr_raw{n}(t,1)),abs(tr_raw{n}(t:end,2)-tr_raw{n}(t,2)));
                %                 indices=cumsum(~(absdrRoughSq<Lscale1));
                %                 timeWithinSquare=sum(~indices)-1;
                %
                %                 timesCrossing=t+(timeWithinDia:timeWithinSquare);
                %                 dr2=(tr_raw{n}(timesCrossing,1)-tr_raw{n}(t,1)).^2+(tr_raw{n}(timesCrossing,2)-tr_raw{n}(t,2)).^2;
                dr2=(tr_raw{n}(t:end,1)-tr_raw{n}(t,1)).^2+(tr_raw{n}(t:end,2)-tr_raw{n}(t,2)).^2;
                indices=cumsum(~(dr2<Lscale2)); %DAK-this appears to be the number of passages of length L
                if indices(end) %if particle has made passage at least once
                    %timeWithin{n}(t)=sum(~indices)-1;
                    timeWithin{n}(t)=sum(~indices);
                    ttoggleLhit(n)=true;
                end
            end
        end
    end
    timeWithin=cell2mat(timeWithin);
    timeWithin=timeWithin(~isnan(timeWithin)&(timeWithin>0));
    %timeWithin=timeWithin/framerate; % first passage time
    timeWithin=timeWithin/Lscale(L)/framerate; % first passage inv vel
    %     timeWithin=1./(timeWithin/Lscale(L)/framerate); % first passage vel
    
    if ~isempty(timeWithin)
        switch 'a'
            case 'a'
                %points=log(0.01)+(0:0.05:1)*(log(max(timeWithin))-log(0.1));
                nbars=28;
                points=log(min(timeWithin)/2)+(0:1/nbars:1)*(log(max(timeWithin)*2)-log(min(timeWithin)/2));
                edges=points(1:2:end);
                S{L}=exp(points(2:2:end));
                N=histc(log(timeWithin),edges);
                N=N(:);
                %                 countingError=sqrt(N');
                %%%%%%%%%%%%%%%%%%%%
                FPP{L}=N(1:end-1)'/NEvents; %normalise with total number of time points
                %FPP{L}=N/length(timeWithin); %normalise with events that reach L only
                disp(['L=' num2str(Lscale(L)) ', number of passages =' num2str(length(timeWithin))])
                %                 h = waitbar(L/length(Lscale));
                %%%%%%%%%%%%%%%%%%%%
                FPP{L}=FPP{L}./diff(exp(edges)); %rescale in terms of prob density
                
            case 'b'
                [FPP{L},S{L}]=hist(timeWithin,20); disp('test')
                
        end
    else
        S{L}=[];
        FPP{L}=[];
    end
end
% close(h)



function editLpix_Callback(hObject, eventdata, handles)
% hObject    handle to editLpix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editLpix as text
%        str2double(get(hObject,'String')) returns contents of editLpix as a double
if isfield(handles,'tracks')
    updatePlots(handles)
end

% --- Executes during object creation, after setting all properties.
function editLpix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editLpix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkboxDist.
function checkboxDist_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxDist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxDist
updatePlots(handles)


% --- Executes on button press in radioFPP.
function radioFPP_Callback(hObject, eventdata, handles)
% hObject    handle to radioFPP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioFPP
if get(handles.radioFPPall,'Value') == 1
    set(handles.radioFPPall,'Value',0)
end
updatePlots(handles)
guidata(hObject,handles)


function updatePlots(handles)
% handles=hObject;
pixel_size = 1/str2double(get(handles.editPixScale,'String'));
frame_rate = str2double(get(handles.editFrameRate,'String'));

nam = handles.nam;
num = str2double(nam(10:end));

%Need to find which particle is numbered num
for i=1:length(handles.tracks)
    if handles.tracks{i}(1,4) == num
        num = i;
        break
    end
end

if get(handles.checkboxDist,'Value') == 1
    % Calculate and plot the distance along a smoothed contour
    Lpix = str2double(get(handles.editLpix,'String'));
    [position,contour,~]=positionAveraging(handles.tracks{num},Lpix);
    time  = 1/frame_rate:1/frame_rate:length(position)/frame_rate;
    plot(handles.axesContour,time,position*pixel_size)
    set(handles.axesContour,'NextPlot','Add')
    plot(handles.axesContour,time(get(handles.slider1,'Value')),position(get(handles.slider1,'Value'))*pixel_size,'rx','MarkerSize',10,'LineWidth',2)
    xlabel(handles.axesContour,'Time (s)')
    ylabel(handles.axesContour,{'Distance along';'contour (\mum)'})
    set(handles.axesContour,'NextPlot','Replace')
end

%Now for the individual track...
set(handles.axesTrack,'NextPlot','Replace')
set(handles.axesTrack,'Visible','On')
% movie = handles.MovieLoc;

if strcmp(handles.MovieFile(end-3:end),'.avi') == 1
    handles.image = read(mmreader(handles.MovieFile), get(handles.slider1,'Value'));
    %     set(handles.figure1,'CurrentAxes',handles.trackAxes)
    imshow(handles.image,'Parent',handles.axesTrack);
else
    handles.image = imread([handles.MoviePath handles.frames{get(handles.slider1,'Value')}]);
    imshow(single(handles.image-min(min(handles.image)))/single(max(max(handles.image))-min(min(handles.image)))*256,gray(256),'Parent',handles.axesTrack);
end

set(handles.axesTrack,'NextPlot','Add')
%Set the axis of the track subplot and insert length scale marker
axis(handles.axesTrack,[min(handles.tracks{num}(:,1))-10 max(handles.tracks{num}(:,1))+10 ...
    min(handles.tracks{num}(:,2))-10 max(handles.tracks{num}(:,2))+10])
rectangle('Position',[min(handles.tracks{num}(:,1))-9,min(handles.tracks{num}(:,2))-9,1/pixel_size,1]...
    ,'facecolor','w','Parent',handles.axesTrack)
text(double(min(handles.tracks{num}(:,1))-6),double(min(handles.tracks{num}(:,2))-6)...
    ,'1{\mu}m','fontsize',14,'color','w','Parent',handles.axesTrack);


%If the distance vs time box has been selected, plot the smooth contour
%first
if get(handles.checkboxDist,'Value') == 1
    plot(handles.axesTrack,contour(:,1),contour(:,2),'m','linewidth',3,'DisplayName','Contour')
end

%Plot the track
plot(handles.axesTrack,handles.tracks{num}(:,1),handles.tracks{num}(:,2),'g','DisplayName','Track')
plot(handles.axesTrack,handles.tracks{num}(get(handles.slider1,'Value'),1),handles.tracks{num}(get(handles.slider1,'Value'),2),'ro','Markersize',20,'linewidth',2)

%Mark the start and end positions of each run
plot(handles.axesTrack,handles.tracks{num}(1,1),handles.tracks{num}(1,2),'yo','markersize',10,'linewidth',2);
text(double(handles.tracks{num}(1,1)+3), double(handles.tracks{num}(1,2)+3), ...
    'Start','fontsize',14,'color','w','Parent',handles.axesTrack);
plot(handles.axesTrack,handles.tracks{num}(end,1),handles.tracks{num}(end,2),'y+','markersize',10,'linewidth',2);
text(double(handles.tracks{num}(end,1)+3),double(handles.tracks{num}(end,2)+3),...
    'End','fontsize',14,'color','w','Parent',handles.axesTrack);



% If the FPP box has been selected, calculate and plot the FPP
if get(handles.checkboxFPP,'value') == 1
    % Calculate and plot the FPP
    if get(handles.radioFPP,'value') == 1
        tempcell{1} = handles.tracks{num};
        titletxt = ['Single track'];
    elseif get(handles.radioFPPall,'value') == 1
        titletxt = 'All tracks';
        n=1;
        for i = 1:length(handles.tracks)
            if ((handles.tracks{i}(end,1)-handles.tracks{i}(1,1))^2 + (handles.tracks{i}(end,2)-handles.tracks{i}(1,2))^2) >= (str2double(get(handles.editMinTrackLength,'String'))*str2double(get(handles.editPixScale,'String')))^2
                tempcell{n} = handles.tracks{i};
                n=n+1;
            end
        end
    end
    % Lscale = logspace(-1,2);
    % Lscale = Lscale(Lscale<tracklength);
    Lscale = [0.01 0.05 0.1 0.5 1 5 10];
    disp(titletxt)
    [S,FPP]=calcFPP(tempcell,frame_rate,pixel_size,Lscale);
    cla(handles.axesFPP)
    colours=jet(length(Lscale));
    for i = 1:length(S)
        plot(handles.axesFPP,S{i},FPP{i},'displayname',['L = ' num2str(Lscale(i)) ' \mum'],'color',colours(i,:),'LineWidth',2)
        set(handles.axesFPP,'NextPlot','Add')
    end
    xlabel(handles.axesFPP,'Inverse speed S  (s\mum^{-1})')
    ylabel(handles.axesFPP,'Probability density  (s{\mu}m^{-1})')
    title(handles.axesFPP,titletxt)
    set(handles.axesFPP,'xscale','log')
    set(handles.axesFPP,'yscale','log')
    legend(handles.axesFPP,'location','best')
end
guidata(handles.figure1, handles);


% --- Executes on button press in radioFPPall.
function radioFPPall_Callback(hObject, eventdata, handles)
% hObject    handle to radioFPPall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radioFPPall
% if get(handles.radioFPP,'Value') == 1
%     set(handles.radioFPP,'Value',0)
% end
updatePlots(handles)
guidata(hObject,handles)


% --- Executes on button press in checkboxFPP.
function checkboxFPP_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxFPP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxFPP
updatePlots(handles)


% --- Executes when selected object is changed in uipanelFPP.
function uipanelFPP_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanelFPP 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
updatePlots(handles)


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
updatePlots(handles)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
