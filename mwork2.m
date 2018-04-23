path = ['/home/mbexegc2/Dropbox (The University of Manchester)/application_support/matlab/code_enquiries/mark_johnston/'];      %set the path of videos
number_vids=1;      %set the number of videos

for ndex = 1:number_vids
    nu_frame=300;         %set the number of frames for each video
    framerate=30;         %set the frame rate (FPS) 
    if ndex == 1
        framerate_name=['HeLaM Lysobrite 1mL 4hrs after staining5_t1'];
    else
        framerate_name=['HeLaM Dextran 0.5mgml overnight' num2str(ndex-1) '_t1'];
    end
    material='data';       %set the name of material

	MVR2 (path,number_vids,nu_frame,framerate,material,framerate_name);
end