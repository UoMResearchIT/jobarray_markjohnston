% INTERPOLATE STITCHES
% Tracks that are stitched together (using Salman Rogers'
% stitchTracksTogether.m function) have NaN values at times between the
% stitches. This can cause some issues with certain calculations (e.g. FPP
% distance calculations). This replaces the NaNs in x and y coordinates
% (first and second columns) with the mean value of the preceeding and
% proceeding coordinates. Therefore, if there are several frames between
% stitches, the particle is stationary half-way between for the duration.
%
% Usage: [particle] = interpolate_stitches(track);

function [particle] = interpolate_stitches_mid(particle)

for i = 2:length(particle)-1
    if isnan(particle(i,1))
        j=0;
        while i+j+1 < length(particle) && isnan(particle(i+j+1,1)) && j<length(particle)
            j = j+1;
        end
        particle(i:i+j,1) = mean([particle(i-1,1) particle(i+j+1,1)]);
        particle(i:i+j,2) = mean([particle(i-1,2) particle(i+j+1,2)]);
    end
end

end