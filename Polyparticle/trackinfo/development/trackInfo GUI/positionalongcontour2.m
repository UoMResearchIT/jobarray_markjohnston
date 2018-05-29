%function position=positionalongcontour2(track,tracksmooth,neighbourhoodmap)
%---------
%for particle track, calculate position along smoothened contour tracksmooth
%neighbourhoodmap is logical array of length(track) by length(tracksmooth), showing which points on track belong to which on tracksmooth
%
function position=positionalongcontour2(track,tracksmooth,neighbourhoodmap)

timeindices=find(~isnan(track(:,1)));
contourLength=[0; cumsum(sqrt(diff(tracksmooth(:,1)).^2+diff(tracksmooth(:,2)).^2))];
position=NaN(size(track(:,1)));
for t=timeindices'
    %find nearest coarse points on smoothed contour
    %indices=find((abs(track(t,1)-tracksmooth(:,1))<neighbourhood)&(abs(track(t,2)-tracksmooth(:,2))<neighbourhood))';
    %find nearest corresponding points on smoothed contour
    indices=find(neighbourhoodmap(t,:));
    dr=zeros(size(indices));
    for N=1:length(indices)
        dr(N)=sqrt((track(t,1)-tracksmooth(indices(N),1))^2+(track(t,2)-tracksmooth(indices(N),2))^2);
    end
    drIndices=sortrows([dr' indices']);
    x1=contourLength(drIndices(1,2));
    d1=drIndices(1,1);
    x2=contourLength(drIndices(2,2));
    d2=drIndices(2,1);
    position(t)=(x2^2-x1^2-d2^2+d1^2)/2/(x2-x1);
end
position=position-position(1);