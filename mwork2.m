path = '/home/mbexegc2/Dropbox (The University of Manchester)/application_support/matlab/code_enquiries/mark_johnston/';      %set the path of videos
videos = dir([path, '*.avi']);
number_vids = length(videos);      %set the number of videos



for video_num = 1:number_vids    
    MVR2 (path, video_num);
end
