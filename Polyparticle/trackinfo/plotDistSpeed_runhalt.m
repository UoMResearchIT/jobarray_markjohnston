function plotDistSpeed_runhalt(handles,nam,dot,window)
%Dot - the location of a dot on the plot showing current time of frame
%Window for 'wavelet' - number of points over which to calculate speed

pixel_size = str2double(get(handles.pixelsizeedit,'String'));
frame_rate = str2double(get(handles.framerateedit,'String'));

% Speed and distance plots
%Calculate distance and speed
initial_x = handles.tracksmooth{str2double(nam(10:end))}(1,1);
initial_y = handles.tracksmooth{str2double(nam(10:end))}(1,2);
% % % Reset positions
tr1 = [];
tr1(:,1) = handles.tracksmooth{str2double(nam(10:end))}(:,1) - initial_x;
tr1(:,2) = handles.tracksmooth{str2double(nam(10:end))}(:,2) - initial_y;
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

[nRuns,runToggle,nHalts,haltToggle]=runhalt(handles.position{str2double(nam(10:end))});

set(handles.distancePlot,'NextPlot','Replace')
time = 1/frame_rate:1/frame_rate:length(tr1(:,3))/frame_rate;

for i = 1:nRuns
    plot(handles.distancePlot,time(find(runToggle==i)),tr1(find(runToggle==i),3),'m','LineWidth',1)
    set(handles.distancePlot,'NextPlot','Add')
end
for i = 1:nHalts
    plot(handles.distancePlot,time(find(haltToggle==i)),tr1(find(haltToggle==i),3),'g','LineWidth',1)
end

plot(handles.distancePlot,time(dot),tr1(dot,3),'rx','MarkerSize',10,'LineWidth',2)
ylabel(handles.distancePlot,{'Distance along'; 'track [{\mu}m]'})
xlabel(handles.distancePlot,'Time [s]')
axis(handles.distancePlot,'tight')

set(handles.speedPlot,'NextPlot','Replace')
if window <= 2
    
    for i = 1:nRuns
        plot(handles.speedPlot,time(find(runToggle==i)),tr1(find(runToggle==i),4),'m','LineWidth',1)
        set(handles.speedPlot,'NextPlot','Add')
    end
    for i = 1:nHalts
        plot(handles.speedPlot,time(find(haltToggle==i)),tr1(find(haltToggle==i),4),'g','LineWidth',1)
    end
    plot(handles.speedPlot,time(dot),tr1(dot,4),'rx','MarkerSize',10,'LineWidth',2)
else speed = movingslope(tr1(:,3),window);
    
    for i = 1:nRuns
    plot(handles.speedPlot,time(find(runToggle==i)),speed(find(runToggle==i))*frame_rate,'m','LineWidth',1)
    set(handles.speedPlot,'NextPlot','Add')
    end
    for i = 1:nHalts
        plot(handles.speedPlot,time(find(haltToggle==i)),speed(find(haltToggle==i))*frame_rate,'g','LineWidth',1)
    end
    plot(handles.speedPlot,time(dot),speed(dot)*frame_rate,'rx','MarkerSize',10,'LineWidth',2)
end
ylabel(handles.speedPlot,{'Speed'; '[{\mu}ms^{-1}]'})
xlabel(handles.speedPlot,'Time [s]')
axis(handles.speedPlot,'tight')
set(handles.speedPlot,'NextPlot','Add')


set(handles.windowSize,'Visible','On')
set(handles.windowSlider,'Visible','On')