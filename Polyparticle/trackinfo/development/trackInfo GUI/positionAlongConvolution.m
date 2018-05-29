%function [displacement,tracksmooth]=positionAlongConvolution(track,w)
%---------
%for a particle track, calculate postion as cumulative sum of displacements
%in the direction of a smoothened contour made by convolution with a window on the scale w
%
function [displacement,tracksmooth]=positionAlongConvolution(track,w)

nanSections=isnan(track(:,3));
timeindices=find(~nanSections);
nanStart=find(diff(nanSections)==1)+1;
nanEnd=find(diff(nanSections)==-1);
gap=nanEnd-nanStart+1;
filledTrack=track(:,1:2);
for n=1:length(nanStart)
    filledTrack(nanStart(n):nanEnd(n),:)=ones(gap(n),1)*filledTrack(nanStart(n)-1,:)+(1:gap(n))'*(filledTrack(nanEnd(n)+1,:)-filledTrack(nanStart(n)-1,:))/(gap(n)+2);
end

%window=gausswin(w,2.5);
window=ones(1,w);
window=window/sum(window);
tracksmooth=[conv(filledTrack(:,1),window) conv(filledTrack(:,2),window)];
tracksmooth(1:w-1,1)=(-w+1:-1)*diff(tracksmooth(w:w+1,1))+tracksmooth(w,1);
tracksmooth(end-w+2:end,1)=(1:w-1)*diff(tracksmooth(end-w:end-w+1,1))+tracksmooth(end-w+1,1);
tracksmooth(1:w-1,2)=(-w+1:-1)*diff(tracksmooth(w:w+1,2))+tracksmooth(w,2);
tracksmooth(end-w+2:end,2)=(1:w-1)*diff(tracksmooth(end-w:end-w+1,2))+tracksmooth(end-w+1,2);

displacement=NaN(size(nanSections));
displacement(1)=0;
for t=timeindices(2:end)'
    contourVector=diff(tracksmooth(t-1:t,:));
    displacement(t)=contourVector*diff(track(t-1:t,1:2))'/sqrt(sum(contourVector.^2));
end
