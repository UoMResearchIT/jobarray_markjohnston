%function [tracksmooth neighbourhoodmap]=coarsecontour2(track,Lpix)
%---
%calculate a smooth version of a particle track, by taking positions separated by a coarsening length Lpix
%
function [tracksmooth neighbourhoodmap]=coarsecontour2(track,Lpix)

%Lpix=4;
%n=str2num(get(gco,'displayname'));
t=1;
tmax=length(track(:,1));
x=track(t,1);
y=track(t,2);
T=NaN(1,100);
N=1;
T(N)=t;
while t<tmax
    N=N+1;
    dr=0;
    while isnan(dr)||((dr<Lpix)&&(t<tmax))
        t=t+1;
        dr=sqrt((x-track(t,1))^2+(y-track(t,2))^2);
    end
    T(N)=t;
    x=track(t,1);
    y=track(t,2);
end
if N<=3
    x=track(end,1)-track(1,1);
    y=track(end,2)-track(1,2);
    norm=sqrt(x^2+y^2);
    x=x/norm*Lpix/2;
    y=y/norm*Lpix/2;
    tracksmooth(:,1)=[track(1,1)-x; track(end,1)+x];
    tracksmooth(:,2)=[track(1,2)-y; track(end,2)+y];
    neighbourhoodmap=true(length(track(:,1)),length(tracksmooth(:,1)));
end
if N>3
    coarsetrackT=T(1:N);
    coarsetrackT=[coarsetrackT(1:end-2) coarsetrackT(end)];
    tracksmooth=zeros(length(coarsetrackT),2);
    neighbourhoodmap=false(length(track(:,1)),length(tracksmooth(:,1)));
    % produce smooth track by mean in window at each point
    for N=2:length(coarsetrackT)-1;
        indices=coarsetrackT(N)+(ceil((coarsetrackT(N-1)-coarsetrackT(N))/2):floor((coarsetrackT(N+1)-coarsetrackT(N))/2));
        %window=ones(size(indices))/length(indices);
        tracksmooth(N,1)=meanNotNaN(track(indices,1));
        tracksmooth(N,2)=meanNotNaN(track(indices,2));
        %neighbourhoodmap(indices,N-1:N+1)=true;
        neighbourhoodmap(indices,N:N+1)=true;
        %neighbourhoodmap(indices,N)=true;
    end
    % extrapolate endpoints
    for N=1;
        indices=coarsetrackT(N)+(0:floor((coarsetrackT(N+1)-coarsetrackT(N))/2));
        tracksmooth(N,1)=3*meanNotNaN(track(indices,1))-2*tracksmooth(N+1,1);
        tracksmooth(N,2)=3*meanNotNaN(track(indices,2))-2*tracksmooth(N+1,2);
        neighbourhoodmap(indices,N:N+1)=true;
    end
    for N=length(coarsetrackT);
        indices=coarsetrackT(N)+(ceil((coarsetrackT(N-1)-coarsetrackT(N))/2):0);
        tracksmooth(N,1)=3*meanNotNaN(track(indices,1))-2*tracksmooth(N-1,1);
        tracksmooth(N,2)=3*meanNotNaN(track(indices,2))-2*tracksmooth(N-1,2);
        neighbourhoodmap(indices,N-1:N)=true;
    end
end
%plot(xsmooth,ysmooth,'y')


function y = meanNotNaN(x)
toggle=~isnan(x);
y = sum(x(toggle))/sum(toggle);
