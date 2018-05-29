function [tr,tr_lst,mmov]=polyparticletracker(filenaming,frameseries,settings,interactive)
% ----------------------------------
% -- Function POLYPARTICLETRACKER --
% ----------------------------------
% by Salman S. Rogers, University of Manchester (2007)
% last edited 5 July 2007
% Documentation included at <a href="polyparticletracker_help.html">polyparticletracker_help.html</a>
%
% Recommended usage:
%   via graphical user interface: see HELP POLYPARTICLETRACKER_GUI
%
% Command line usage: 
%   [tr,tr_lst,mmov]=polyparticletracker(filenaming,frameseries,settings,interactive)
%
% Inputs:
%   FILENAMING: specifies filenames. Input files may either be (1) a single 
%               AVI file (readable by Matlab) or (2) a series of frame images.
%
%           (1) - Input single AVI file:
%               FILENAMING should be a structure with single field AVIFILE
%               specifying input filename, ending in ".avi".
%
%           (2) - Input series of frame images: 
%               FILENAMING should be a structure with fields 
%               FILEPREFIX, FILESUFFIX and NZEROS. 
%  
%               fileprefix - e.g. 'directory' or 'directory/name'
%  
%               filesuffix - e.g. '.tif'
%  
%               nzeros - e.g. '7' for 0000001.tif. Note: nzeros is a string!
%  
%   FRAMESERIES: one-dimensional array specifying frame numbers, e.g. (1:100)
%
%   SETTINGS: structure with allowed members:
%  
%               lnoise - smoothing lengthscale lnoise. Default: 1
%  
%               N - look for new particles every N frames. Default: 200
%  
%               mintracklength - discard any tracks lasting for fewer 
%                     frames. Default: 200
%  
%               brightnessmin - brightness minimum threshhold. Default: 0
%  
%               brightnessmax - brightness maximum threshhold. Default: 64 
%                     (intensities are scaled to a floating point number between 0 and 64)
%  
%               Sthresh - exclude particles with greater skewness parameter. 
%                     Default: 1.0
%  
%               eccthresh - exclude particles with greater eccentricity.
%                     Default: 0.85
%  
%               radiusmax - exclude particles with greater radius. Default: 10 
%  
%               radiusmin - exclude particles with lesser radius. Default: 1 
%  
%               subpixpeak - use peak of fitted surface to find sub-pixel 
%                     centre correction. Set to true or false. Default: true for PFGW (gives better results). 
%                     If false, find the centroid of the image within the particle, with the same Gaussian 
%                     weighting function. 
%  
%               darkness - set to 1 if particles are bright on a dark 
%                     background or 0 for dark on light background. Set to
%                     2 to search for both (default).
%  
%               intlims - two element array containing lowest and highest 
%                     brightness values to scale image. If unset, brightness 
%                     scaled automatically.
%
%   INTERACTIVE: set to 1 to produce movie output and show tracking. Set 
%               to 0 otherwise
%
%
% Outputs:
%   TR: cell array, where each cell is an array containing the track and 
%               measurements for each particle. (Format of array detailed 
%               below.)
%
%   TR_LST: complete array with all identified particles in all frames.
%               (Format of array detailed below.)
%
%   MMOV: movie of tracked particles with background image. The movie is 
%               created only when INTERACTIVE (above) is set to 1.
%
% Output array format: 
%   Each output array is a table with the following columns:
%       Column 1: X position in whole image (pixels)
%       Column 2: Y position in whole image (pixels)
%       Column 3: frame number
%       Column 4: particle number - unique number identifying particle
%       Column 5: radius (pixels) 
%       Column 6: J
%       Column 7: eccentricity - equals 0 for circles, 1 for lines/ridges
%       Column 8: rotation (radians) of particle's major axis
%       Column 9: average brightness of pixels within particle
%       Column 10: skewness
%
%
% See also
%
% <a href="polyparticletracker_help.html">polyparticletracker_help</a> and <a href="polyparticletracker_reference.html">polyparticletracker_reference</a>
tic
%Settings:
limmult=1.6; %multiplies size of box around particle for fitting -- not a sensitive parameter - should be slightly over one... say 1.6
maskmult=1; %multiplies fall-off of weighting exponential in fine fit  -- should be about 1.0
improverthresh=0.5; %repeat centre refinement cycle if x or y correction is greater than improverthresh

% smoothing lengthscale lnoise
if ~ismember('lnoise',fieldnames(settings));  lnoise=1; else lnoise=settings.lnoise; end  
% look for new particles every N frames
if ~ismember('N',fieldnames(settings));  N=200; else N=settings.N; end  
% discard any tracks shorter than mintracklength
if ~ismember('mintracklength',fieldnames(settings));  mintracklength=200; else mintracklength=settings.mintracklength; end
% brightness min threshhold
if ~ismember('brightnessmin',fieldnames(settings));  brightnessmin=0; else brightnessmin=settings.brightnessmin; end  
% brightness max threshhold
if ~ismember('brightnessmax',fieldnames(settings));  brightnessmax=64; else brightnessmax=settings.brightnessmax; end  
% exclude particles with Skewness greater than Sthresh
if ~ismember('Sthresh',fieldnames(settings));  Sthresh=1; else Sthresh=settings.Sthresh; end  
% exclude particles with eccentricity greater than eccthresh
if ~ismember('eccthresh',fieldnames(settings));  eccthresh=0.85; else eccthresh=settings.eccthresh; end  
% exclude particles with radius greater than radiusmax
if ~ismember('radiusmax',fieldnames(settings));  radiusmax=10; else radiusmax=settings.radiusmax; end  
% exclude particles with radius less than radiusmin
if ~ismember('radiusmin',fieldnames(settings));  radiusmin=1; else radiusmin=settings.radiusmin; end  
% darkness = 0 or 1 if particles are dark or bright against background ---- usu darkness = 2 for both bright and dark particles
if ~ismember('darkness',fieldnames(settings));  darkness=2; else darkness=settings.darkness; end  
% intlims --- two element array containing lowest and highest brightness values to scale image
if ~ismember('intlims',fieldnames(settings)); intlims=[]; else intlims=double(settings.intlims); end  
% subpixmethod --- default: subpixel centre is peak of fitted image
if ~ismember('subpixpeak',fieldnames(settings)); subpixpeak=true; else subpixpeak=settings.subpixpeak; end  
% starthandle --- use startbutton userdata to send messages to running program
if ~ismember('starthandle',fieldnames(settings)); starthandle=false; else starthandle=settings.starthandle; end  

% EDGE - number of pixels to disregard from edge
EDGE=max([1 6*round(lnoise)]);

% nearthresh - distance to search for same particle moved out of its own radius
nearthresh=radiusmax;


% A) initialise
% do rough particle find
% do fine particle find
% select good particles
% allow plotting and manual selecting if interactive mode is on

%turn off all warnings
warning off all


%select first frame
% - either fileprefix/suffix/nzeros specified (maybe just the fileprefix)
% - or avifile specified
n=1;
if ismember('fileprefix',fieldnames(filenaming))
    avifilereader=false;
    filename=[filenaming.fileprefix num2str(frameseries(n),['%0' filenaming.nzeros '.0f']) filenaming.filesuffix];
    info=imfinfo(filename);
elseif ismember('avifile',fieldnames(filenaming))
    avifilereader=true;
    avibatch=100;
    info=aviinfo(filenaming.avifile);
    nseq=1;
    %Nseq=ceil(length(frameseries)/avibatch);
    %
    %check for different avifilenumbers
    if ismember('avifilenum',fieldnames(filenaming))
        multiavi=true;
        avisequences=[];
        cumavifileframes=[0 cumsum(filenaming.avifileframes)];
        firstavi=find(diff(frameseries(1)<=cumavifileframes));
        lastavi=find(diff(frameseries(end)<=cumavifileframes));
        for currentavi=firstavi:lastavi
            framesincurrent=frameseries((frameseries>cumavifileframes(currentavi))&(frameseries<=cumavifileframes(currentavi+1)))-cumavifileframes(currentavi);
            padframes=mod(length(framesincurrent)-1,avibatch)+1;
            avisequences=[avisequences framesincurrent NaN(1,avibatch-padframes)];
        end
        avisequences=reshape(avisequences,avibatch,[]);
        %
        currentavi=firstavi;
        if isempty(currentavi)
            warning('Frame number out of bounds')
        end
        filenaming.avifile(end-5:end-4)=num2str(filenaming.avifilenum(currentavi),'%02d');
    else
        multiavi=false;
        padframes=mod(length(frameseries),avibatch);
        avisequences=reshape([frameseries NaN(1,avibatch-padframes)],avibatch,[]);
    end
    % define n_nextseq - values of n where the next set of frames are read
    temp=avisequences;
    temp(~isnan(temp))=0;
    temp(1,:)=1;
    n_nextseq=find(temp(~isnan(temp)));
    %
    avisequences_now=avisequences(~isnan(avisequences(:,nseq)),nseq);
    avibatch_now=length(avisequences_now);
    disp(['Reading ' filenaming.avifile ' ...'])
    disp(['Extracting AVI file frames ' num2str(avisequences_now(1)) ' to ' num2str(avisequences_now(end)) ' ...'])
    avidata=aviread(filenaming.avifile,avisequences_now);
    %temp=mmread(filenaming.avifile,avisequences_now);avidata=temp.frames;
else
    error('Filenaming does not specify correct filenames')
end

% select region of interest
if ismember('roi',fieldnames(settings))
    roi=settings.roi;
else
    roi=[1 info.Width 1 info.Height];
end

frame=readandfilter;

% set sizes and coordinates
sizeroi=[roi(4)-roi(3)+1 roi(2)-roi(1)+1];
x=(1:sizeroi(2));
y=(1:sizeroi(1));
theta=(0:0.05:1)*2*pi;

% pick noise level - exclude peaks smaller than mean
if ~ismember('noiselevel',fieldnames(settings))
    %noiselevel=mean(frame(:));
    littleframe=frame(EDGE:sizeroi(1)-EDGE,EDGE:sizeroi(2)-EDGE);
    %fftframe=fft(littleframe(:));
    %noiselevel=sum(abs(fftframe(1:max(sizeroi))))/sizeroi(1)/sizeroi(2);
    noiselevel=mean(littleframe(:))+mean(abs(diff(littleframe(:))))*10;
else
    noiselevel=settings.noiselevel;
end  

if interactive;
    %m(n)={uint8(frame)};
    if ~ismember('outputaxes',fieldnames(settings));
        figure
        axeshandle=gca;
    else
        axeshandle=settings.outputaxes;
    end
    imhandle=image(frame,'Parent',axeshandle);
    %colormap('jet')
    colormap('gray')
    set(axeshandle,'XLim',[1 sizeroi(2)])
    set(axeshandle,'YLim',[1 sizeroi(1)])
    hold on
    drawnow
    axis manual
end

% rough particle find
[npeaks,centreroughx,centreroughy,radiusrough,lightnotdark,Jrough,brightnesspeakrough]=roughparticlefind_wrap;
particleson=1:npeaks;
good=initgoodness;
particleson=find(good);
% if interactive
%     for i=particleson
%         [xcirc,ycirc]=pol2cart(theta,radiusrough(i));
%         plot(xcirc+centreroughx(i),ycirc+centreroughy(i),'g')
%     end
%     drawnow
% end

% fine particle find
% and select good particles
[centrex,centrey,J,eccentricity,brightness,rotation,skewness,radius,lightnotdark,brightnesspeak]=fineparticlefind(centreroughx,centreroughy,radiusrough,lightnotdark,Jrough,brightnesspeakrough);
good=goodness(good);
particleson=find(good);

% initialise particlenumbers for particles that will be output
particlenumbers=particleson;
switchedoff=false(1,length(particlenumbers));
if isempty(switchedoff)
    switchedoff=[];
end

%initialise tracks array: tr_lst
max_tr_lst_size=1000000;
tr_lst=single(zeros(max_tr_lst_size,10)); %1000000*10 double array takes 80 Mb
if ~isempty(particleson)
    tr_lst(1:sum(good),:)=[centrex(particleson)' centrey(particleson)' frameseries(n)*ones(length(particleson),1) particleson' radius(particleson)' J(particleson)' eccentricity(particleson)' rotation(particleson)' brightness(particleson)' skewness(particleson)'];
end
lastrow_tr_lst=sum(good);

% use initplotter to plot particles for the first time
particlehandles=initplotter;

%***clear cache file***
delete('polyparticletracker_cache_*.dat')
cachenum=1;

% B) loop
% in each iteration:
% - check particles still good
% - check particle hasn't changed much
for n=2:length(frameseries)
    
    % use avi reader or individual frames
    if avifilereader
        %if mod(n,avibatch)==1
        if ismember(n,n_nextseq)
            %nseq=ceil(n/avibatch);
            nseq=nseq+1;
            avisequences_now=avisequences(~isnan(avisequences(:,nseq)),nseq);
            %avibatch_now=length(avisequences_now);
            %temp=mmread(filenaming.avifile,avisequences_now);avidata=temp.frames;            
            %
            if multiavi
                %check whether to load next file
                if avisequences_now(1)==1
                    currentavi=currentavi+1;
                    filenaming.avifile(end-5:end-4)=num2str(filenaming.avifilenum(currentavi),'%02d');
                end
            end
            disp(['Reading ' filenaming.avifile ' ...'])
            disp(['Extracting AVI file frames ' num2str(avisequences_now(1)) ' to ' num2str(avisequences_now(end)) ' ...'])
            avidata=aviread(filenaming.avifile,avisequences_now);
        end
    else
        filename=[filenaming.fileprefix num2str(frameseries(n),['%0' filenaming.nzeros '.0f']) filenaming.filesuffix];
    end
    frame=readandfilter;
    
    disp(n)
    % in each N interations, look for new particles outside radii of existing particles
    if rem(n,N)==0
        %**** store existing particle values ****
        exgood=good;
        exparticleson=particleson;
        excentrex=centrex;
        excentrey=centrey;
        exradius=radius;
        exlightnotdark=lightnotdark;
        exJ=J;
        exbrightnesspeak=brightnesspeak;
        spacer=zeros(1,length(good)-length(centrex));
        clear good particleson
        
        % now look for new particles
        [npeaks,centreroughx,centreroughy,radiusrough,lightnotdark,Jrough,brightnesspeakrough]=roughparticlefind_wrap;
        % discard peaks which overlap with existing particles
        %clear iswithin
        iswithin=zeros(1,length(npeaks));
        %iswithin=[];
        for i=1:npeaks
            %iswithin(i)=1;
            iswithin(i)=sum(sqrt((centrex(exparticleson)-centreroughx(i)).^2+(centrey(exparticleson)-centreroughy(i)).^2)<1.2*max([radius(exparticleson) radiusrough(i)]));
        end
        particleson=find(~iswithin);
        good=initgoodness;
        particleson=find(good);
        % fine particle find
        % and select good particles
        if ~isempty(particleson)
            [centrex,centrey,J,eccentricity,brightness,rotation,skewness,radius,lightnotdark,brightnesspeak]=fineparticlefind(centreroughx,centreroughy,radiusrough,lightnotdark,Jrough,brightnesspeakrough);
            good=goodness(good);
            particleson=find(good);
        end

        % use initplotter to plot particles for the first time
        if interactive
            particlehandlesnew=initplotter;
            particlehandles.centre=[particlehandles.centre particlehandlesnew.centre];
            particlehandles.ellipse=[particlehandles.ellipse particlehandlesnew.ellipse];
        end        
        
        % update particlenumbers, particleson, good for NEW particles
        centrex=[excentrex spacer centrex(particleson)];
        centrey=[excentrey spacer centrey(particleson)];
        radius=[exradius spacer radius(particleson)];
        lightnotdark=[exlightnotdark spacer lightnotdark(particleson)];
        J=[exJ spacer J(particleson)];
        brightnesspeak=[exbrightnesspeak spacer brightnesspeak];
        %
        newparticlenumbers=length(exgood)+find(particleson);
        newswitchedoff=false(1,length(newparticlenumbers));
        particlenumbers=[particlenumbers newparticlenumbers];
        switchedoff=[switchedoff newswitchedoff];
        particleson=[exparticleson newparticlenumbers];
        good=[exgood ones(1,length(newparticlenumbers))];
        
        
    end

    % now track all good particles - numbered by particleson
    if ~isempty(particleson)
        [centrex,centrey,J,eccentricity,brightness,rotation,skewness,radius,lightnotdark,brightnesspeak]=fineparticlefind(centrex,centrey,radius,lightnotdark,J,brightnesspeak);
        good=goodness(good);
        particleson=particleson(ismember(particleson,find(good)));
    end
    if ~isempty(particleson)
        % select good particles
        tr_lst(lastrow_tr_lst+(1:sum(good)),:)=[centrex(particleson)' centrey(particleson)' frameseries(n)*ones(length(particleson),1) particleson' radius(particleson)' J(particleson)' eccentricity(particleson)' rotation(particleson)' brightness(particleson)' skewness(particleson)'];
    end
    numevents=sum(good);
    lastrow_tr_lst=lastrow_tr_lst+numevents;
    %if lastrow_tr_lst greater than 1e6, then save data to hard disk
    if lastrow_tr_lst > (max_tr_lst_size-2*numevents)
        save_tr_lst=tr_lst(1:lastrow_tr_lst,:);
        cachename=['polyparticletracker_cache_' num2str(cachenum) '.dat'];
        disp(['Saving cache file ''' cachename '''...'])
        save(cachename,'save_tr_lst','-mat')
        cachelength(cachenum)=lastrow_tr_lst;
        cachenum=cachenum+1;
        clear save_tr_lst
        tr_lst=zeros(max_tr_lst_size,10);
        lastrow_tr_lst=0;
    end
    
    plotter;
    %pause momentarily to allow interrupt if necessary
    pause(0.00000001)
    if starthandle
        if get(starthandle,'userdata')
            disp('Execution stopped!')
            set(starthandle,'Userdata',0) %switch toggle off again
            tr={};
            tr_lst=tr_lst(1:lastrow_tr_lst,:);
            return
        end
    end
end

% switch warnings back on
warning on all

save_tr_lst=tr_lst(1:lastrow_tr_lst,:);
clear tr_lst;

tr={};
%particles={};
retainer=[];
discarder=[];
% create output matrices
if cachenum>1
    % cache final data
    cachename=['polyparticletracker_cache_' num2str(cachenum) '.dat'];
    disp(['Saving final cache file ''' cachename '''...'])
    save(cachename,'save_tr_lst','cachenum','cachelength','lastrow_tr_lst','particlenumbers','roi','mintracklength','-mat')
    clear save_tr_lst

    % loop through cache
    %tr_lst=single(zeros(sum(cachelength)+lastrow_tr_lst,10));
    cachelengthcum=[0 cumsum(cachelength) sum(cachelength)+lastrow_tr_lst];
    tr_lst=cell(1,cachenum);
    for m=1:cachenum
        cachename=['polyparticletracker_cache_' num2str(m) '.dat'];
        disp(['Loading cache file ''' cachename '''...'])
        load(cachename,'save_tr_lst','-mat');
        %tr_lst((1+cachelengthcum(m)):cachelengthcum(m+1),:)=save_tr_lst;
        
        %discard short tracks as each file is loaded
        particlenumberstotest=setdiff(unique(save_tr_lst(:,4)),union(retainer,discarder));
        for n=1:length(particlenumberstotest)
            
            events_in_current=sum(save_tr_lst(:,4)==particlenumberstotest(n));
            %1. does particle have enough events in range to retain? yes retain no (2) (3)
            %2. is particle present at begninning => does it have enough events in both ranges to keep? yes retain no discard
            %3. is particle at end? yes leave no discard
            if events_in_current>=mintracklength
                retainer=union(retainer,particlenumberstotest(n));
            else
                exist_at_beginning=ismember(particlenumberstotest(n), save_tr_lst(save_tr_lst(:,3)==save_tr_lst(1,3),4));  
                if exist_at_beginning 
                    if m>1
                        events_in_current=events_in_current+sum(tr_lst{m-1}(:,4)==particlenumberstotest(n));
                    end
                    if events_in_current>=mintracklength
                        retainer=union(retainer,particlenumberstotest(n));
                    else
                        discarder=union(discarder,particlenumberstotest(n));
                    end
                elseif ~ismember(particlenumberstotest(n), save_tr_lst(save_tr_lst(:,3)==save_tr_lst(end,3),4));  
                    %i.e. if doesn't exist at end then discard
                    discarder=union(discarder,particlenumberstotest(n));
                end
            end
                    
        end
        %remove tracks with less than mintracklength points
        rows=ismember(save_tr_lst(:,4),discarder);
        tr_lst(m)={save_tr_lst(~rows,:)};
    end
    tr_lst=vertcat(tr_lst{1:end});
else
    tr_lst=save_tr_lst;
    clear save_tr_lst
end

disp('Creating output arrays...')
if ~isempty(tr_lst)
    %put x and y into proper coordinate frame
    tr_lst(:,1)=tr_lst(:,1)+roi(1);
    tr_lst(:,2)=tr_lst(:,2)+roi(3);

    %extract individual tracks
    discarder=[];
    for n=1:length(particlenumbers)
        rows=tr_lst(:,4)==particlenumbers(n);
        if sum(rows)>=mintracklength
            tr=[tr {tr_lst(rows,:)}];
        else
            discarder=[discarder particlenumbers(n)];
        end
    end
    %remove tracks with less than mintracklength points
    rows=ismember(tr_lst(:,4),discarder);
    tr_lst=tr_lst(~rows,:);

    %extract individual frames
    %for n=0:length(frameseries)-1
    %    particles(n+1)={tr_lst(tr_lst(:,3)==n,:)};
    %end
    %particles='not set';
else
    disp('Warning: No particles found!')
end

disp('Done!')
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to read and filter each frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function frameb=readandfilter
        if avifilereader
            frameb=frame2im(avidata(n-n_nextseq(nseq)+1));
        else
            frameb=imread(filename);
        end
        if length(roi)==1
            frameb=mean(single(frameb),3);
        else
            frameb=mean(single(frameb(roi(3):roi(4),roi(1):roi(2),:)),3);
        end
        %frameb=(frameb-min(min(frameb)))/(max(max(frameb))-min(min(frameb)))*64;
        if isempty(intlims)
            intlims=[min(min(frameb)) max(max(frameb))].*[0.8 1.2];
        end
        frameb=(frameb-intlims(1))/(intlims(2)-intlims(1))*64;
        %if darkness==1
        %    %frameb=frameb;     %dark background
        %else
        %    frameb=64-frameb;   %light background... try imcomplement
        %end
        %if lnoise>0
            frameb=lpass(double(frameb),lnoise); % lnoise = length scale of noise -- normally = 1
        %end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function lpass 
% - Implements a real-space low pass filter which suppresses pixel noise (from Crocker/Grier/Weeks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function arr_g = lpass(arr,lnoise)
  if lnoise>0
      w = 3*round(lnoise);

      % Gaussian Convolution kernel
      sm = 0:2*w;
      r = (sm - w)/(2 * lnoise);
      gx = exp( -r.^2) / (2 * lnoise * sqrt(pi));
      gy = gx';

      res = arr;
      g = conv2(res,gx,'valid');
      tmpg = g;
      g = conv2(tmpg,gy,'valid');
      arr_g = zeros(size(arr));
      arr_g((w+1):end-w,(w+1):end-w) = g; 
  else
      arr_g=arr;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wrapper functions - choose light/dark particles or both
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [npeaks,centreroughx,centreroughy,radiusrough,lightnotdark,Jrough,brightnesspeakrough]=roughparticlefind_wrap
    switch darkness
        case 2 %track both light and dark
            [npeaks1,centreroughx1,centreroughy1,radiusrough1,Jrough1,brightnesspeakrough1]=roughparticlefind(frame);
            [npeaks2,centreroughx2,centreroughy2,radiusrough2,Jrough2,brightnesspeakrough2]=roughparticlefind(64-frame);
            npeaks=npeaks1+npeaks2;
            centreroughx=[centreroughx1 centreroughx2];
            centreroughy=[centreroughy1 centreroughy2];
            radiusrough=[radiusrough1 radiusrough2];
            Jrough=[Jrough1 Jrough2];
            brightnesspeakrough=[brightnesspeakrough1 brightnesspeakrough2];
            lightnotdark=[true(1,npeaks1) false(1,npeaks2)];
        case 0 %track dark particles (i.e. normal particles in focus in bright field)
            [npeaks,centreroughx,centreroughy,radiusrough,Jrough,brightnesspeakrough]=roughparticlefind(64-frame);
            lightnotdark=false(1,npeaks);
        case 1 %track light particles (i.e. normal in fluorescence)
            [npeaks,centreroughx,centreroughy,radiusrough,Jrough,brightnesspeakrough]=roughparticlefind(frame);
            lightnotdark=true(1,npeaks);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to find particles (rough)
% - go through all peaks, fit with quartic surface in square with limits given by 4-way turning points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [npeaks,centreroughx,centreroughy,radiusrough,Jrough,brightnesspeakrough]=roughparticlefind(frame)

    %pick local maxima
    peaksind=find(imregionalmax(frame));
    [peaksy,peaksx]=ind2sub(sizeroi,peaksind);
    peakbrightness=frame(peaksind);
    %[peaksy,peaksx,peakbrightness]=find(imregionalmax(frame));
    %take only those further than EDGE from edge
    %take only peaks which have a brightness greater than noise level and
    %greater than brightnessmin, less than brightnessmax
    peakson=and(peaksx>EDGE,peaksx<sizeroi(2)-EDGE) ...
            & and(peaksy>EDGE,peaksy<sizeroi(1)-EDGE) ...
            ...& (peakbrightness>noiselevel) ...
            & (peakbrightness>=brightnessmin) ...
            & (peakbrightness<=brightnessmax);
    peaksx=peaksx(peakson);
    peaksy=peaksy(peakson);
    npeaks=sum(peakson);
    disp(['Examining ' num2str(npeaks) ' peaks...'])
    
    centreroughx=zeros(1,npeaks);
    centreroughy=zeros(1,npeaks);
    radiusrough=zeros(1,npeaks);
    Jrough=zeros(1,npeaks);
    brightnesspeakrough=zeros(1,npeaks);
    radiusmax3=3*radiusmax;
    
    for i=1:npeaks
        slicex=frame(peaksy(i),:);
        slicey=frame(:,peaksx(i));
        tpointsxl=find(imregionalmax(diff(slicex)));
        tpointsxr=find(imregionalmin(diff(slicex)))+1;
        tpointsyl=find(imregionalmax(diff(slicey)));
        tpointsyr=find(imregionalmin(diff(slicey)))+1;

        %find adjacent turning points
        limx=[max(tpointsxl(tpointsxl<peaksx(i))) min(tpointsxr(tpointsxr>peaksx(i)))];
        limy=[max(tpointsyl(tpointsyl<peaksy(i))) min(tpointsyr(tpointsyr>peaksy(i)))];
        
        if and(diff(limx)<radiusmax3, diff(limy)<radiusmax3)
            
            %define coordinates around peaks (xp,yp) and (xpp,ypp)
            xp=limx(1):limx(2);
            xpp=xp-peaksx(i);
            yp=limy(1):limy(2);
            ypp=yp-peaksy(i);

            %fit quartic surface within limits, weighted by difference between intensity and minimum within limits, squared
            minwithinlim=min([frame(limy(1),limx(1)) frame(limy(1),limx(2)) frame(limy(2),limx(1)) frame(limy(2),limx(2))]);
            p=polyfitweighted2(xpp,ypp,frame(yp,xp),4,(frame(yp,xp)-minwithinlim).^2+1);

            %calculate centre based on quadric part p(1:6)
            %note p = [p00 p10 p01 p20 p11 p02 p30 p21 p12 p03...]
            %but following Wolfram Mathworld's ellipse equation: a=p20 b=p11/2 c=p02 d=p10/2 f=p01/2 g=p00
            a=p(4); b=p(5)/2; c=p(6); d=p(2)/2; f=p(3)/2; g=p(1);
            Jrough(i)=a*c-b^2;
            xc=(b*f-c*d)/Jrough(i);
            yc=(b*d-a*f)/Jrough(i);
            centreroughx(i)=peaksx(i)+xc;
            centreroughy(i)=peaksy(i)+yc;
            %centreroughx(i)=peaksx(i)+sign(xc)*min([abs(xc) improverthresh]);
            %centreroughy(i)=peaksy(i)+sign(yc)*min([abs(yc) improverthresh]);

            %decide limits of particle by when intensity is sufficiently different from quadric, and quadric>0
            quadric=polyval2(p(1:6),xpp,ypp);
            %inparticle=and(frame(yp,xp)<(polyval2(p(1:6),xcorrection,ycorrection)+quadric)/1.9, quadric>0);
            inparticle=and(frame(yp,xp)<polyval2(p(1:6),xc,yc)/3+2*quadric/2, quadric>0);
            lxpp=length(xpp);
            lypp=length(ypp);
            ypp2=ypp'*ones(1,lxpp);
            xpp2=ones(lypp,1)*xpp;
            pixradiisquare=(xpp2(inparticle)-xc).^2+(ypp2(inparticle)-yc).^2;
            radiusrough(i)=sqrt(sum(sum(pixradiisquare.^2))/sum(sum(pixradiisquare))); %weighted by the distance to the pixel - prob better for finding edge                
            brightnesspeakrough(i)=g;
        else
            centreroughx(i)=peaksx(i);
            centreroughy(i)=peaksy(i);
            radiusrough(i)=radiusmax3;
        end 
    end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to find particles (fine)
% - now we have radiusrough, calculate again using gaussian weighting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [centrex,centrey,J,eccentricity,brightness,rotation,skewness,radius,lightnotdark,brightnesspeak]=fineparticlefind(centreroughx,centreroughy,radiusrough,lightnotdark,Jrough,brightnesspeakrough)
    centrex=zeros(1,max(particleson));
    centrey=zeros(1,max(particleson));
    J=zeros(1,max(particleson));
    eccentricity=zeros(1,max(particleson));
    brightness=zeros(1,max(particleson));
    brightnesspeak=zeros(1,max(particleson));
    rotation=zeros(1,max(particleson));
    skewness=zeros(1,max(particleson));
    radius=zeros(1,max(particleson));
    lightnotdark=lightnotdark(1:max(particleson));
    %for i=1:find(radiusrough<radiusmax)
    for i=particleson
        % repeat subpixel correction until satisfactory
        doextracycles=3;
        stillgood=true;
        refinementcount=0;
        maxrefinements=limmult*radiusrough(i)/improverthresh; %if it takes more than maxrefinements to find centre correction, then move on
        while stillgood && doextracycles
            refinementcount=refinementcount+1;
            [p,quadric,frameypxp,weightexp,xpp2,ypp2]=polyfitgaussweight(centreroughx(i),centreroughy(i),radiusrough(i));

            %calculate centre based on quadric part p(1:6)
            %note p = [p00 p10 p01 p20 p11 p02 p30 p21 p12 p03...]
            %but following Wolfram Mathworld's ellipse equation: a=p20 b=p11/2 c=p02 d=p10/2 f=p01/2 g=p00
            a=p(4); b=p(5)/2; c=p(6); d=p(2)/2; f=p(3)/2; g=p(1);
            J(i)=a*c-b^2;

            %if J(i) negative, try to rescue track by looking for a neighbouring peak of similar J, nearer than next nearest particle
            if J(i)<0
                indices=setdiff(particleson,i);
                nearestparticle=sqrt(min((centreroughx(indices)-centreroughx(i)).^2+(centreroughy(indices)-centreroughy(i)).^2)/2);
                nearx=min(abs([nearthresh nearestparticle-radiusrough(i)]));
                neary=min(abs([nearthresh nearestparticle-radiusrough(i)]));
                nearlimx=round([max([centreroughx(i)-nearx 1]) min([centreroughx(i)+nearx x(end)])]);
                nearlimy=round([max([centreroughy(i)-neary 1]) min([centreroughy(i)+neary y(end)])]);
                nearxp=nearlimx(1):nearlimx(2);
                nearyp=nearlimy(1):nearlimy(2);
                sizenearroi=[length(nearyp) length(nearxp)];
                                
                %indices = sub2ind(sizeroi,nearyp,nearxp);
                if lightnotdark(i)
                    nearpeakind=find(imregionalmax(frame(nearyp,nearxp)));
                else
                    nearpeakind=find(imregionalmin(frame(nearyp,nearxp)));
                end
                [nextnearesty,nextnearestx]=ind2sub(sizenearroi,nearpeakind);                
                nextnearestx=nextnearestx+nearlimx(1)-1;
                nextnearesty=nextnearesty+nearlimy(1)-1;
                %[peakbrightness,nearbest]=max(frame(sub2ind(sizeroi,nextnearesty,nextnearestx)));
                peakbrightness=frame(sub2ind(sizeroi,nextnearesty,nextnearestx));
                [temp,nearbest]=min((abs(peakbrightness-brightnesspeakrough(i))+1).*((nextnearestx-centreroughx(i)).^2+(nextnearesty-centreroughy(i)).^2));
                centreroughx(i)=nextnearestx(nearbest);
                centreroughy(i)=nextnearesty(nearbest);
                [p,quadric,frameypxp,weightexp,xpp2,ypp2]=polyfitgaussweight(centreroughx(i),centreroughy(i),radiusrough(i));
                a=p(4); b=p(5)/2; c=p(6); d=p(2)/2; f=p(3)/2; g=p(1);
                J(i)=a*c-b^2;
                %reject if J has appreciably changed
                %if abs((J(i)-Jrough(i))/Jrough(i)) > 5
                %    J(i)=-1;
                %end
            end
                
            %choose method for subpixel calculation
            %use peak of fitted surface
            xc=(b*f-c*d)/J(i);
            yc=(b*d-a*f)/J(i);
            %refine with 3rd order terms using 1st order perturbation
            %d1=d+(3*p(7)*xc^2 + 2*p(8)*xc*yc + p(9)*yc^2)/2;
            %f1=f+(p(8)*xc^2 + 2*p(9)*xc*yc + 3*p(10)*yc^2)/2;
            %refine with 3rd, 4th order terms using 1st order perturbation
            %d1=d+(3*p(7)*xc^2 + 2*p(8)*xc*yc + p(9)*yc^2 + 4*p(11)*xc^3 + 3*p(12)*xc^2*yc + 2*p(13)*xc*yc^2 + p(14)*yc^3)/2;
            %f1=f+(p(8)*xc^2 + 2*p(9)*xc*yc + 3*p(10)*yc^2 + p(12)*xc^3 + 2*p(13)*xc^2*yc + 3*p(14)*xc*yc^2 + 4*p(15)*yc^3)/2;
            %xc=(c*d1-b*f1)/(b^2-a*c);
            %yc=(a*f1-b*d1)/(b^2-a*c);
        
            centreroughx(i)=centreroughx(i)+sign(xc)*min([abs(xc) improverthresh]);
            centreroughy(i)=centreroughy(i)+sign(yc)*min([abs(yc) improverthresh]);

            rotation(i)=0.5*acot((c-a)/2/b)+(a-c<0)*pi/2;
            ct=cos(rotation(i));
            st=sin(rotation(i));
            P1=p(11)*ct^4-p(12)*ct^3*st+p(13)*ct^2*st^2-p(14)*ct*st^3+p(15)*st^4;
            P2=p(11)*st^4+p(12)*st^3*ct+p(13)*st^2*ct^2+p(14)*st*ct^3+p(15)*ct^4;
            Q1=p(4)*ct^2-p(5)*ct*st+p(6)*st^2;
            Q2=p(4)*st^2+p(5)*st*ct+p(6)*ct^2;
            %radiusrough(i)=sqrt(mean([-Q1/6/P1 -Q2/6/P2]));
            radiusrough(i)=abs(sqrt(sqrt(Q1*Q2/P1/P2/36))); %geometric mean
            
            % if |xc| or |yc| greater than improverthresh, move centrerough by up to a pixel and calculate again
            % stop recalculating if J becomes negative or if finding is taking too long
            stillgood=(refinementcount<=maxrefinements)&&(J(i)>0); %if not still good, stop at once... 
            improverswitch=(abs(xc)>improverthresh)||(abs(yc)>improverthresh); % check if xc,yc above thresh or extra cycles required
            doextracycles = doextracycles - ~improverswitch; %if still good and improverswitch turns off, do extra cycles
        
        end
        
        radius(i)=radiusrough(i);
        
        %decide limits of particle by when intensity is sufficiently different from quadric, and quadric>0
        %inparticle=and(frameypxp<(polyval2(p(1:6),xcorrection,ycorrection)+quadric)/1.9, quadric>0);
        %inparticle=and(frameypxp<(polyval2(p(1:6),xcorrection,ycorrection)+3*quadric)/4, quadric>0);
        %inparticle=and((frameypxp-quadric)<2, quadric>0);

        %pixradiisquare=(xpp2(inparticle)-xcorrection).^2+(ypp2(inparticle)-ycorrection).^2;
        %radius(i)=sqrt(sum(sum(pixradiisquare.*frame(inparticle)))/brightness(i)); %weighted by the brightness of the pixel - prob better for finding peak
        %radius(i)=sqrt(sum(sum(pixradiisquare.^2))/sum(sum(pixradiisquare))); %weighted by the distance to the pixel - prob better for finding edge    
        %radius(i)=sqrt(abs(a+c)/sum(abs(p(7:15))));
        %radius(i)=sqrt(J(i)/mean(abs(p(11:15)./[1 1 1 1 1])));
        %radius(i)=sqrt(mean([-a/(6*p(11)+3*p(12)+p(13)) -c/(6*p(15)+3*p(14)+p(13))]));
        %radius(i)=sqrt(mean([-a/6/p(11) -c/6/p(15)]));
        semiaxes=sort(sqrt(2*(a*f^2+c*d^2+g*b^2-2*b*d*f-a*c*g)/-J(i)./([c-a a-c]*sqrt(1+4*b^2/(a-c)^2)-c-a)));
        eccentricity(i)=sqrt(1-semiaxes(1)^2/semiaxes(2)^2);
        
        xedge=pol2cart(0,2*radius(i)./sqrt(cos(rotation(i)).^2+sin(rotation(i)).^2/(1-eccentricity(i).^2))/(1+sqrt(1-eccentricity(i)^2)));
        inparticle=quadric>a*xedge^2+2*d*xedge+g;

        if ~subpixpeak
            %calculate sub-pixel corrections based on centroid of image within particle
            weight=frameypxp.*inparticle.*weightexp;
            sumsumweight=sum(sum(weight));
            xc=sum(sum(xpp2.*weight))/sumsumweight;
            yc=sum(sum(ypp2.*weight))/sumsumweight;
        end

        centrex(i)=centreroughx(i)+xc;
        centrey(i)=centreroughy(i)+yc;
        %centrex(i)=centreroughx(i)+min([xc xc2]);
        %centrey(i)=centreroughy(i)+min([yc yc2]);

        %brightness = mean of fitted surface in particle, weighted by brightness at each pixel
        brightness(i)=sum(sum(quadric(inparticle).^2))/sum(sum(quadric(inparticle)));
        brightnesspeak(i)=g;
        %brightness = RMS of fitted surface in particle
        %brightness(i)=sqrt(mean2(quadric.^2.*inparticle)); 
        if ~lightnotdark(i)
            brightness(i)=64-brightness(i);
            brightnesspeak(i)=64-g;
        end
        %the third order terms are all skewed => use as rough measure of skewness, defined as ratio of cubic to quadric terms at x,y=radius
        %skewness(i)=sum(abs(p(7:10)))*radius(i); 
        skewness(i)=sum(abs(p(7:10)))*radius(i)/J(i); 
    
    end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate the Gaussian-weighted polynomial fit of the particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [p,quadric,frameypxp,weightexp,xpp2,ypp2]=polyfitgaussweight(roughx,roughy,roughradius)
        %define coordinates around centrerough - (xp,yp) and (xpp,ypp)
        limx=round([max([roughx-limmult*roughradius 1]) min([roughx+limmult*roughradius x(end)])]);
        limy=round([max([roughy-limmult*roughradius 1]) min([roughy+limmult*roughradius y(end)])]);
        xp=limx(1):limx(2);
        xpp=xp-roughx;
        yp=limy(1):limy(2);
        ypp=yp-roughy;
        %fit quartic surface within limits, weighted by exponential around centrerough
        lxpp=length(xpp);
        lypp=length(ypp);
        xpp2=ones(lypp,1)*xpp;
        ypp2=ypp'*ones(1,lxpp);
        weightexp=exp(-(xpp2.^2+ypp2.^2)/roughradius.^2*maskmult);
        
        if lightnotdark(i)
            frameypxp=frame(yp,xp);
        else
            frameypxp=64-frame(yp,xp);
        end
        
        %p=polyfitweighted2(xpp,ypp,frameypxp,4,weightexp.*frameypxp);
        p=polyfitweighted2(xpp,ypp,frameypxp,4,weightexp);

        quadric=polyval2(p(1:6),xpp,ypp);
        %quartic=polyval2(p(1:15),xpp,ypp);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to rate initial goodness
% perform tests to evaluate goodness of peak
% - peak must be further than EDGE from edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function good=initgoodness
    %clear good
        if isempty(particleson)
            good=[];
        else
            good(particleson)=(centreroughx(particleson)>EDGE) & (centreroughx(particleson)<x(end)-EDGE) ...
                & (centreroughy(particleson)>EDGE).*(centreroughy(particleson)<y(end)-EDGE) ... 
                & isreal(radiusrough(particleson)) ...
                & (radiusrough(particleson)<radiusmax);
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to rate goodness
% perform tests to evaluate goodness of peak
% - peak must have positive J=ac-b^2 ... in fact set this to some small value to exclude the tiniest peaks
% - peak must be further than EDGE from edge
% - skewness or egginess must be small compared to J
% - discard if within 2 times the radius of a brighter particle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function good=goodness(good)
    if ~isempty(particleson)
        good=false(1,length(good));
        %for i=particleson
        %    isradiusreal(i)=isreal(radius(i));
        %end
        isradiusreal=~imag(radius);
        good(particleson)= ...
             (centrex(particleson)>EDGE) & (centrex(particleson)<x(end)-EDGE) ...
             & (centrey(particleson)>EDGE) & (centrey(particleson)<y(end)-EDGE) ...
             & (J(particleson)>0) ...
             & (skewness(particleson)<Sthresh) ...
             & (eccentricity(particleson)<eccthresh) ...
             & (radius(particleson)>radiusmin) ...
             & (radius(particleson)<radiusmax) ...
             & (brightness(particleson)>=brightnessmin) ...
             & (brightness(particleson)<=brightnessmax) ...
             & isradiusreal(particleson) ...
           ...  & (skewness(particleson)<sqrt(J(particleson))) ...
             ;
        for i=particleson
            ibr = good(particleson) & ((lightnotdark(i) & (brightness(particleson)>brightness(i))) | (~lightnotdark(i) & (brightness(particleson)<brightness(i))));
            %iswithin(i)=sum(sqrt((centrex(ibr)-centrex(i)).^2+(centrey(ibr)-centrey(i)).^2)<radius(i)+radius(ibr));
            iswithin(i)=sum(sqrt((centrex(particleson(ibr))-centrex(i)).^2+(centrey(particleson(ibr))-centrey(i)).^2)<1.2*radius(i));
        end
        good(particleson)=good(particleson) & ~iswithin(particleson);
        %good=true(1,length(good));
        %radius(particleson)
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial function to plot particles
% - if interactive is on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function particlehandles=initplotter
    if interactive
        set(imhandle,'cdata',frame);
        particlehandles.centre=[];
        particlehandles.ellipse=[];
        for i=1:length(particleson)
            %contour(xp,yp,frame(yp,xp))
            particlehandles.centre(i)=plot(centrex(particleson(i)),centrey(particleson(i)),'.');
            [xellipse,yellipse]=pol2cart(theta-rotation(particleson(i)), ...
                2*radius(particleson(i)) ...
                    ./sqrt(cos(theta).^2+sin(theta).^2/(1-eccentricity(particleson(i)).^2))...
                    /(1+sqrt(1-eccentricity(particleson(i))^2)));
            particlehandles.ellipse(i)=plot(xellipse+centrex(particleson(i)),yellipse+centrey(particleson(i)),'k');
            %contour(xp,yp,inparticle,[0.5 0.5])
        end
        drawnow
        mmov(1)=getframe(axeshandle);
    else
        particlehandles=[];
    end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subsequent function to plot particles
% - if interactive is on
% - for all particleson, set new plot data
% - delete particles that go off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function plotter
    if interactive
        if ~isempty(switchedoff)
            %disp(particleson)
            set(imhandle,'cdata',frame);
            for i=1:length(particleson)
                ii=(particlenumbers==particleson(i));
                %contour(xp,yp,frame(yp,xp))
                set(particlehandles.centre(ii),'xdata',centrex(particleson(i)));
                set(particlehandles.centre(ii),'ydata',centrey(particleson(i)));
                [xellipse,yellipse]=pol2cart(theta-rotation(particleson(i)),2*radius(particleson(i))./sqrt(cos(theta).^2+sin(theta).^2/(1-eccentricity(particleson(i)).^2))/(1+sqrt(1-eccentricity(particleson(i))^2)));
                set(particlehandles.ellipse(ii),'xdata',xellipse+centrex(particleson(i)));
                set(particlehandles.ellipse(ii),'ydata',yellipse+centrey(particleson(i)));
                %contour(xp,yp,inparticle,[0.5 0.5])
            end
            switchoff=and(~ismember(particlenumbers,particleson),~switchedoff); %switch off those not already switched off
            switchedoff=or(switchedoff,switchoff);
            delete(particlehandles.centre(switchoff))
            delete(particlehandles.ellipse(switchoff))
            %particlehandles.centre(switchoff)=0;
            %particlehandles.ellipse(switchoff)=0;
            drawnow
            mmov(n)=getframe(axeshandle);
        end
    else
        mmov(n)=0;
    end
    end

end