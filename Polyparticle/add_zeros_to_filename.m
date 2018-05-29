% Add zeros to numbered image files
%
% I had a sequence of tiff images in such a format:
%
% image1.tif, image2.tif...image200.tif, image201.tif etc.
%
% What I want is the following:
%
% image00001.tif, image00002.tif...image00200.tif, image00201.tif
%
% such that all files are the same size.

string = 'Control-90225-11a_t';

 for i = 1:1000
image = imread([string num2str(i) '.tif']);
imwrite(image, ['test\' strcat([string '_'], sprintf('%04d',i)) '.tif'],'tif','Compression','None')
end