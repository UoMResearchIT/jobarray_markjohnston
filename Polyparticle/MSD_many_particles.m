% mean square displacement for many particles
% DAK
% Inputs:
% tracks - cell array of tracked particles a la PolyParticleTracker
% pixscale - microns per pixel
% framerate - frames per second
% logtime - set to 1 for approximately logarithmically-spaced time (else
% zero). Logtime is faster as there will be fewer data points
% N.B. May need to manually adjust the time parameter to be more
% appropriate to your time frame
function [MSD, time] = MSD_many_particles(tracks,pixscale,framerate,logtime)

% find the maximum time
endtime = zeros(length(tracks),1);
for i = 1:length(tracks)
    endtime(i) = tracks{i}(end,3);
end

if logtime==1
time = space125(1,max(endtime));
else
    time = 1:max(endtime);
end

MSD = nan(length(tracks),length(time));

%         For a range of times, what is the average displacement squared?
% Range of times:
for j = 1:length(tracks)
    disp(num2str(j))
    n=1; %Counter
    x = tracks{j}(:,1)*pixscale;
    y = tracks{j}(:,2)*pixscale;
    
    for  t = time(time < length(tracks{j}));
        %     disp(num2str(t))
        displacement2 = [];
        for i = 1:length(tracks{j})-t
            
            displacement2(i) = ((x(i+t) - x(i)).^2 + (y(i+t) - y(i)).^2);
        end
        MSD(j,n) = mean(displacement2);
        n=n+1;
    end
    
end

figure
plot(time(1:size(MSD,2))/framerate,nanmean(MSD),'.-')
xlabel('Time [s]')
ylabel('MSD [\mum^2]')