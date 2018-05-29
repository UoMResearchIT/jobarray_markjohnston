%Program that goes through the particle tracks and extracts only those that
%are considered to be long runs.
%Need to determine what is a 'long' run (in terms of pixels)

% tracks = tracks2;

function [runs] = find_runs(tracks,run_min)

% run_min = 10;
runs = [];
lengths = [];

for i = 1:size(tracks,2)    %i.e. the number of tracked particles
    % a^2 + b^2 = c^2
    % sqrt((x2-x1) + (y2-y1)) >= run_min
    length = sqrt((abs((tracks{1,i}(end,1) - tracks{1,i}(1,1))))^2 + (abs((tracks{1,i}(end,2) - tracks{1,i}(1,2))))^2);
    lengths = [lengths length];
    if length >= run_min
        runs = [runs i];
        
    end
end

%Can then plot the trajectories of those particles considered to make a run

% for i = runs
%     figure,
%     plot(tracks{1,i}(:,1),tracks{1,i}(:,2)), hold on
%     plot(tracks{1,i}(1,1),tracks{1,i}(1,2),'g+','LineWidth',3)
%     plot(tracks{1,i}(end,1),tracks{1,i}(end,2),'r+','LineWidth',3)
%     title(['Track no. ' num2str(i)])
% end