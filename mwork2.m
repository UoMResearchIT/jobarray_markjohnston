path = '/home/mbexegc2/Dropbox (The University of Manchester)/application_support/matlab/code_enquiries/mark_johnston/';      %set the path of videos
videos = dir([path, '*.avi']);
number_vids = length(videos);      %set the number of videos

nu_frame=300;         %set the number of frames for each video
material='data';       %set the name of material

for ndex = 1:number_vids    
    if ndex == 1
        framerate_name='HeLaM Lysobrite 1mL 4hrs after staining5_t1';
    else
        framerate_name=['HeLaM Dextran 0.5mgml overnight' num2str(ndex-1) '_t1'];
    end

    MVR2 (path,ndex,nu_frame,material,framerate_name);
end
