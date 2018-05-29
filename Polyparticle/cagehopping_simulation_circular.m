% Cage-hopping simulation
% DAK 2012
% A particle is confined within a box. The particle moves randomly unless
% it hits the box boundary, at which there is a chance it will hop to another box.
% Otherwise the particle bounces.
%
% settings.radius - cage radius
% settings.timeLength - the number of times the particle moves (frames)
% interactive = 0;% settings.initial_hop_prob - the initial probability of a hop each time the
%                               boundary is reached (increases with number
%                               of failed hops (bounces)
% settings.xstep - the particle step length in the x direction
% settings.ystep - the particle step length in the y direction
% interactive  - set to 1 to display ther particle bouncing 'live'

function [x,y] = cagehopping_simulation_circular(settings,interactive)

initial_hop_prob = settings.initial_hop_prob; %Probability that a hop will occur at boundary

% cage edges (initial starting points)
radius = settings.radius;
box_width = radius;
box_height = radius;
xmin_box = -box_width/2;
xmax_box = box_width/2;
ymin_box = -box_height/2;
ymax_box = box_height/2;

% Determine step size of the particle
xstep = settings.xstep;
ystep = settings.ystep;
% % Determine box step size (for cage hopping there is not step)
% xbox_step = 0;
% ybox_step = 0;

%Create a time variable (can have intervals less than or greater than 1)
time = 1:settings.timeLength;

x = nan(length(time),1);
y = nan(length(time),1);
x(1) = 0;  %The start position
y(1) = 0;

% interactive = 0;

if interactive == 1
    figure
    
    % Plot the box
    boxhandle = rectangle('Position',[xmin_box, ymin_box, (xmax_box-xmin_box), (ymax_box-ymin_box)],'Curvature',[1,1],'edgecolor','r','linewidth',2); %Box outline, x y width height
    hold on
    pointhandle = plot(x,y);
    axis([-100 100 -100 100]);
%     gridxy([-5*box_width/2:box_width:5*box_width/2],[-5*box_height/2:box_height:5*box_height/2],'Linestyle','--')
    drawnow
end

colours = lines(100);
c=1;

hop_prob = initial_hop_prob;
hopped = 0; %marker for whether the particle has hopped
n = 1; % Number of times it hasn't hopped
bounced = 0;
centre = [0,0];
for i = 2:length(time)
    
    
    location = [x(i-1)-centre(1) y(i-1)-centre(2)];
    [theta,rho] = cart2pol(location(1),location(2));
    
    %Check if the particle has passed the  edge
    if rho > radius
        %Check if the particle hops
        if rand < hop_prob
            centre = centre+location;
            hopped = 1;
            c=c+1;
            x(i) = x(i-1) + randn*xstep;      %Otherwise random +/- displacement
            y(i) = y(i-1) + randn*ystep;      %Otherwise random +/- displacement
        %Otherwise bounce
        else
        [xbounce,ybounce] = pol2cart(theta-pi,abs(randn(1))/4*xstep);
        x(i) = x(i-1)+xbounce;
        y(i) = y(i-1)+ybounce;
        bounced = 1;
        end
    %Otherwise randomly diffuse
    else
        x(i) = x(i-1) + randn*xstep;      %Otherwise random +/- displacement
        y(i) = y(i-1) + randn*ystep;      %Otherwise random +/- displacement
    end

    
    if hopped == 0
        if bounced == 1
        hop_prob = hop_prob+(hop_prob/n);
        n = n+1;
        end
        bounced = 0;
    elseif hopped == 1
        hop_prob = initial_hop_prob; %Reset probability
        n = 1; %Reset count
        hopped = 0; %Reset hop
    end
    
    %     %Move the boundary
    %     xmin_box = xmin_box + xbox_step;
    %     xmax_box = xmax_box + xbox_step;
    %     ymin_box = ymin_box + ybox_step;
    %     ymax_box = ymax_box + ybox_step;
    
    if interactive == 1
%         delete(boxhandle) %Clear the previous box
        delete(pointhandle) %Clear the previous trace
        % Plot the box
        boxhandle = rectangle('Position',[centre(1)-radius, centre(2)-radius, radius*2, radius*2],'Curvature',[1,1],'edgecolor',colours(c,:),'linewidth',2); %Box outline, x y width height
        hold on
        pointhandle = plot(x,y);
        drawnow
%             pause(0.01)
        
    end
    
end