% Plot tracks on image

%% For a .TIF file

image = imread(['Control-90225-15a_t1.tif']);   %Loads the tif file
figure
imshow(single(image-min(min(image)))/single(max(max(image))-min(min(image)))*256,gray(256));
% hold on

%% For an .avi file

image = read(mmreader('10kHz_2.avi'), 1);  %Loads the first frame of the avi file
figure
imshow(image);
hold on

%% Plot the tracks on the image

% load('data_10kHz_1_Thrashed','tracks') %Where 'movie_tracks' is the MAT file and 'tracks' is the variable
% tracks = data.tr;
for i = 1:length(tracks) % Or choose which particle you want to show
    plot(tracks{i}(:,1),tracks{i}(:,2),'r')
end

%% Plot a scale bar

x = 50;
y = 25;
w = 5.998*5;
h = 3;

rectangle('Position',[x,y,w,h]...
    ,'facecolor','w')
text(x,y+10,...
    '5{\mu}m','fontsize',14,'color','w');