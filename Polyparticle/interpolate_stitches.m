% INTERPOLATE STITCHES
% Tracks that are stitched together (using Salman Rogers'
% stitchTracksTogether.m function) have NaN values at times between the
% stitches. This can cause some issues with certain calculations (e.g. FPP
% distance calculations). This function linearly interpolates between the preceeding and
% proceeding coordinates. It also fills in the missing frame numbers if there is a 3rd column.
%
% Usage: [track] = interpolate_stitches(track);
% OR
% for i = 1:length(tracks)
%   tracks{i} = interpolate_stitches(tracks{i});
% end

function [track] = interpolate_stitches(track)

x = track(:,1);
y = track(:,2);


x(isnan(x)) = interp1(find(~isnan(x)), x(~isnan(x)), find(isnan(x)),'linear');   %Could change this to cubic...
y(isnan(y)) = interp1(find(~isnan(y)), y(~isnan(y)), find(isnan(y)),'linear');

track(:,1) = x;
track(:,2) = y;

if size(track,2)>2
t = track(:,3);
t(isnan(t)) = interp1(find(~isnan(t)), t(~isnan(t)), find(isnan(t)),'linear');
track(:,3) = t;
end

end