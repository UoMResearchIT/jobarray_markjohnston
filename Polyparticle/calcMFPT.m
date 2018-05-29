%calcMFPT - calculates the mean first passage time DAK
% usage [MFPT, FPT]=calcMFPT(tracks,framerate,pixscale,Lscale)
% input: tracks - cell array of tracks
%        framerate - frames per second
%        pixscale - microns per pixel
%        Lscale - array of length scales at which to calculate MFPT
%
function [MFPT, FPT]=calcMFPT(tr_raw,framerate,pixscale,Lscale)

tracktoggle=true(size(tr_raw));

MFPT=zeros(size(Lscale));
FPT =cell(size(Lscale));

parfor L=1:length(Lscale)
    timeWithin=cell(length(tr_raw),1);
    Lscale2=(Lscale(L)/pixscale)^2;
    ttoggleLhit=false(length(tr_raw),1);
    for n=find(tracktoggle)
        timeWithin{n}=NaN(size(tr_raw{n}(:,1)));
        for t=1:size(tr_raw{n},1)
            dr2end=(tr_raw{n}(end,1)-tr_raw{n}(t,1)).^2+(tr_raw{n}(end,2)-tr_raw{n}(t,2)).^2;
            if dr2end>Lscale2 %if particle has made passage at least once
                dr2=(tr_raw{n}(t:end,1)-tr_raw{n}(t,1)).^2+(tr_raw{n}(t:end,2)-tr_raw{n}(t,2)).^2;
                indices=cumsum(~(dr2<Lscale2));
                if indices(end) %if particle has made passage at least once
                    timeWithin{n}(t)=sum(~indices);
                    ttoggleLhit(n)=true;
                end
            end
        end
    end
    timeWithin=cell2mat(timeWithin);
    timeWithin=timeWithin(~isnan(timeWithin)&(timeWithin>0));
    timeWithin=timeWithin/framerate; % first passage time
    
    if ~isempty(timeWithin)
        MFPT(L) = mean(timeWithin);
        FPT{L} = timeWithin;
    end
    
    disp(L)
end
