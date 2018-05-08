function MVR2 (path, video_num)
% MVR2 Process video files using PolyParticleTracker
% This function has been rewritten to make it suitable for calling from a
% job array on the CSF.
% Inputs:
%   path: directory containing the video files to process
%   video_num: which video number to process (This is the SGE_TASK_ID
%   when using a job array).
% Requires:
%   Poly particle tracker toolbox

    % Video settings
    nu_frame = 300;         % Set the number of frames for each video
    material = 'data';      % Set the name of material
    
    % Job-array variables
    videos = dir([path, '*.avi']);          % Get list of all videos
    video_name = videos(video_num).name;    % Get current video name
    filenaming{video_num}.avifile = [path video_name]; % Full path to video file
    
    %%%%%%%%%%%%%%%%%%%%% PolyParticleTracker settings %%%%%%%%%%%%%%%%%%%%%
    settings.lnoise = 1;                  %Smoothing lengthscale lnoise
    settings.N = 2;                       %Look for new particles every N frames
    settings.mintracklength = 10;         %framerate/2;  %Discard tracks less than so many frames (usually 0.5sec)
    settings.Sthresh = 1;                 %Exclude particles with greater skewness parameter
    settings.eccthresh = 0.85;            %Eccentricity threshold
    settings.darkness = 1;                %Search for light and dark particles, 1 is light on dark, 0 dark on light and 2 is for both.
    settings.radiusmin = 1;               %Set min radius in pixels
    settings.radiusmax = 30;
    settings.roi = [1 512 1 512];
    movieLength = nu_frame;
    data.tr = cell(1,1);
    
    
    %name the file
    disp (filenaming{video_num}.avifile);
    
    [tr,tr_lst,mmov]=polyparticletracker_parallelx(filenaming{video_num},...
        1:movieLength,...               %Specifies which frames are used (in this case all of them)
        settings,...                    %The settings as defined above
        0);                             %Interactive value - set to zero to not display movie
    
    tracks{video_num} = tr;
    num_tracks= length(tracks{video_num});
    
    for track_num = 1:num_tracks
        vidnum = repmat(video_num,[length(tracks{video_num}{track_num}(:,1)),1]);
        data.tr{track_num} = [tracks{video_num}{track_num},vidnum];
    end
    
    video_name = replace(video_name, '.avi', '');
    videodataname = [material '_' video_name '.mat'];
    save(videodataname);
end
