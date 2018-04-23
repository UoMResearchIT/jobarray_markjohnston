function MVR2 (path,number_vids,nu_frame,framerate,material,framerate_name)
    
    %%%%%%%%%%%%%%%%%%%%% PolyParticleTracker settings %%%%%%%%%%%%%%%%%%%%%
    settings.lnoise = 1;                  %Smoothing lengthscale lnoise
    settings.N = 2;                      %Look for new particles every N frames
    settings.mintracklength = 10;         %framerate/2;  %Discard tracks less than so many frames (usually 0.5sec)
    settings.Sthresh = 1;                   %Exclude particles with greater skewness parameter
    settings.eccthresh = 0.85;               %Eccentricity threshold
    settings.darkness = 1;                  %Search for light and dark particles, 1 is light on dark, 0 dark on light and 2 is for both.
    settings.radiusmin = 1;                 % set min radius in pixels
    settings.radiusmax = 30;
    settings.roi = [1 512 1 512];
    movieLength = nu_frame;
    data.tr = cell(1,1);
    for i=1:number_vids
        
        if number_vids == 1
            movie = framerate_name;
            filenaming{i}.avifile = [path movie '.avi'];
        else movie = [framerate_name '_' num2str(i)];
            filenaming{i}.avifile = [path movie '.avi'];
        end
        
    end
    
    
    %matlabpool('open','4');
    
    for i= 1:number_vids
        
        %name the file
        disp (filenaming{i}.avifile);
        movieObj = VideoReader(filenaming{i}.avifile);
        %settings.N = movieLength;                                        %specifies number of frames
        %settings.mintracklength = movieLength;
        [tr,tr_lst,mmov]=polyparticletracker_parallelx(filenaming{i},...
            1:movieLength,...               %Specifies which frames are used (in this case all of them)
            settings,...                    %The settings as defined above
            0);                             %Interactive value - set to zero to not display movie
        
        %nu_tracks = length(data.tr);
        tracks{i} = tr;
        
    end
    %matlabpool('close');
    data.tr=[];
    for i= 1:number_vids
        
        nu_tracks = length(data.tr);
        nu_tracks2= length(tracks{i});
        
        for j = 1:nu_tracks2
            nu_tracks = nu_tracks+1;
            vidnum = repmat([i],[length(tracks{i}{j}(:,1)),1]);
            data.tr{nu_tracks} = [tracks{i}{j},vidnum];
        end
        
    end
    
    
    
    videodataname = [material '_' framerate_name];
    save (videodataname);
    
    clear s j d how_many mmov tr_lst tracks tr path settings movieObj movieLength filenaming number_vids nu_tracks nu_tracks2 movie i;
end