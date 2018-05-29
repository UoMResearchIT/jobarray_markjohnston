function MSD_single_particle(track,framerate,pixscale)

% mean square displacement for a single particle
% framerate = 1;
% pixscale = 0.167;

        x = track(:,1)*pixscale;
        y = track(:,2)*pixscale;
        
                n=1; %Counter
                MSD = [];
                time = 1:200;
%                 time = [1 2 3 4 5 6 7 8 9 10 12 14 16 18 20 25 30 40 50 60 70 80 90 100 120 140 160 180 200 250 300 400 500 600 700 800 900 1000 1200 1400 1600 1800 2000 2500 3000 4000 5000 6000 7000 8000 9000 10000];
% time = space125(1,10000);
%         For a range of times, what is the average displacement squared?
        % Range of times:

for  t = time;
%     disp(num2str(t))
    displacement2 = [];
    for i = 1:length(track)-t

        displacement2(i) = ((x(i+t) - x(i)).^2 + (y(i+t) - y(i)).^2);
    end
    MSD(n,1) = mean(displacement2);
    n=n+1;
end

figure
loglog(time/framerate,MSD,'.-')
xlabel('Time')
ylabel('MSD [\mum^2]')