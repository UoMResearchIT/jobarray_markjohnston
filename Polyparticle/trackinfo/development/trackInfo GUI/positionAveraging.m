%function [position,contour,tracksmooth]=positionAveraging(track,Lpix)
%---------
%for a particle track, calculate position at each point as mean position within Lpix
%
function [position,contour,tracksmooth]=positionAveraging(track,Lpix)

nanSections=isnan(track(:,3));
timeindices=find(~nanSections);
nanStart=find(diff(nanSections)==1)+1;
nanEnd=find(diff(nanSections)==-1);
gap=nanEnd-nanStart+1;
filledTrack=track(:,1:2);
for n=1:length(nanStart)
    filledTrack(nanStart(n):nanEnd(n),:)=ones(gap(n),1)*filledTrack(nanStart(n)-1,:)+(1:gap(n))'*(filledTrack(nanEnd(n)+1,:)-filledTrack(nanStart(n)-1,:))/(gap(n)+2);
end


%produce smooth track by averaging positions within Lpix
Lpix2=Lpix^2;
tracksmooth=zeros(length(filledTrack),2);
%times=zeros(length(filledTrack),1);
for t=1:length(nanSections)
    dr2=(filledTrack(:,1)-filledTrack(t,1)).^2+(filledTrack(:,2)-filledTrack(t,2)).^2;
    indices=cumsum(dr2>Lpix2);
    indices=indices==indices(t);
    tracksmooth(t,:)=mean(filledTrack(indices,:),1);
    %times(t)=mean(find(dr2<Lpix2));
end

%produce contour by taking course contour along smoothened track
[contour, neighbourhoodmap]=coarsecontour2(tracksmooth,Lpix);
if isnan(sum(contour(:)))
    contour=[2*track(1,1)-track(end,1) 2*track(1,2)-track(end,2); 2*track(end,1)-track(1,1) 2*track(end,2)-track(1,2)];
end
%position=positionalongcontour2(tracksmooth,contour,neighbourhoodmap);
position=positionalongcontour2(track,contour,neighbourhoodmap);

%position=[0; cumsum(sqrt(sum(diff(tracksmooth).^2,2)))];
