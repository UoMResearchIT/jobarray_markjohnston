% Follow Cheezum et al. (2001) procedure for simulated particle images
% Generate 1000 frames
% (Script written by SS Rogers (2007))
%
% N.B. A new directory, 'sim_images' is generated in the current directory,
% in which the tiff images are saved

% first make high res "precise" image - 9nm pixel size
% add circular particle of defined radius
% convolve with PSF
highres=zeros(1001);
%one pix particle - "point source"
highres(501,501)=10000;

%convolve with point spread function (PSF)
x=(-200:200)*9e-9;
o=ones(size(x));
radius=sqrt((x'*o).^2+(o'*x).^2);
a=2*pi*1.3/570e-9;
PSF=(2*besselj(1,radius*a)./radius).^2;
PSF(201,201)=a^2;
highrestemp=conv2(highres,PSF);
highres2=sparse([],[],[],3401,3401);
highres2(1001:2401,1001:2401)=highrestemp;
imagemask=1201:2201;

% Create 1000 low resolution frames:
% - shift data by whole number of high res pixels
% - integrate over CCD cells
% - scale signal and background
% - add Poisson noise
% - write image file

%generate a particle track - in whole displacements of the fine pixels
%initial particle position is image centre
shiftx=zeros(1,1000);
shifty=zeros(1,1000);
for t=2:1000
    % - use Gaussian distribution for random numbers, mean=0, stdev=3
    shiftx(t)=shiftx(t-1)+round(random('Normal',0,3,1));
    shifty(t)=shifty(t-1)+round(random('Normal',0,3,1));
end

figure
plot(shiftx,shifty,'r')
drawnow

%scale to given signal (and background=10)
SignaltoBackground=32;
Signal=10*SignaltoBackground;
SignaltoNoise=(Signal-10)/sqrt(Signal);
disp(['SignaltoNoise = ' num2str(SignaltoNoise)])

%generate images
mkdir sim_images
for t=1:1000
    highres3=highres2(imagemask-shiftx(t),imagemask-shifty(t));

    %integrate over blocks (CCD 'cells') of 11x11 pixels - such that unshifted image has peak in the centre of a block
    lowres=zeros(91);
    for n=1:91
        for m=1:91
            lowres(n,m)=mean2(highres3((n-1)*11+(1:11),(m-1)*11+(1:11)));
        end
    end
    
    lowres2=lowres / max(max(lowres)) * (Signal-10) + 10;
    
    %add poisson noise
    %use Image Toolbox function IMNOISE - note: input values are interpreted scaled by 1e12
    lowres3 = imnoise(lowres2/1e12,'poisson')*1e12;

    %write image
    filename=['sim_images/frame' num2str(t,'%0.4d') '.tif'];
    imwrite(uint16(lowres3),filename,'tif');
    disp(t)
end

%rescale x and y trajectory for CCD cell grid
exacttrajectoryx=shifty/11+46;
exacttrajectoryy=shiftx/11+46;
