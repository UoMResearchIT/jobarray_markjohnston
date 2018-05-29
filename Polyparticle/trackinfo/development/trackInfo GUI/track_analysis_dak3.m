function [meanSpeed,totalDisplacement,lifetime,indices] = track_analysis_dak3(stitched_tracks,position,displacementThresh)

tr = stitched_tracks;

% 1. select particles with contour length greater than totalDisplacementThresh pixels
totalDisplacementThresh=displacementThresh; %23.44=2um, 11.72=1um or Set value,cascade

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