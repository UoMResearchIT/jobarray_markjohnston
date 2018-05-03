path = '/home/mbexegc2/Dropbox (The University of Manchester)/application_support/matlab/code_enquiries/mark_johnston/';      %set the path of videos
videos = dir([path, '*.avi']);
number_vids = length(videos);      %set the number of videos

nu_frame=300;         %set the number of frames for each video
material='data';       %set the name of material

for video_num = 1:number_vids    
    video_name = videos(video_num).name;
    MVR2 (path,video_num,nu_frame,material,video_name);
end
