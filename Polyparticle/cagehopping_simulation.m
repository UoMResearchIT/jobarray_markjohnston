% Cage-hopping simulation
% DAK 2012
% A particle is confined within a box. The particle moves randomly unless
% it hits the box boundary, at which there is a chance it will hop to another box.
% Otherwise the particle bounces.
%
% Inputs:
% settings.cageWidth - width of the cage
% settings.cageHeight - height of the cage
% settings.timeLength - the number of times the particle moves (frames)
% settings.initial_hop_prob - the initial probability of a hop each time the
%                               boundary is reached (increases with number
%                               of failed hops (bounces)
% settings.xstep - the particle step length in the x direction
% settings.ystep - the particle step length in the y direction
% interactive  - set to 1 to display ther particle bouncing 'live'


function [x,y] = cagehopping_simulation(settings,interactive)

initial_hop_prob = settings.initial_hop_prob; %Probability that a hop will occur at boundary

% cage edges (initial starting points)
box_width = settings.cageWidth;
box_height = settings.cageHeight;
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
    boxhandle = rectangle('Position',[xmin_box,ymin_box, (xmax_box-xmin_box), (ymax_box-ymin_box)],'edgecolor','r','linewidth',2); %Box outline, x y width height
    hold on
    pointhandle = plot(x,y);
    axis([-30 30 -30 30]);
    gridxy([-5*box_width/2:box_width:5*box_width/2],[-5*box_height/2:box_height:5*box_height/2],'Linestyle','--')
    drawnow
end

hop_prob = initial_hop_prob;
hopped = 0; %marker for whether the particle has hopped
n = 1; % Number of times it hasn't hopped
bounced = 0;
for i = 2:length(time)
    
    if x(i-1) < xmin_box        %If the minimum boundary has been reached, chance that it will hop to next cage
        if rand < hop_prob, x(i) = xmin_box - abs(randn)*xstep;    %Hop to next cage
            xmax_box = xmin_box;
            xmin_box = xmin_box - box_width;
            hopped = 1;
        else
            x(i) = xmin_box + abs(randn)*xstep;    %Bounce
            bounced = 1;
        end
    elseif x(i-1) >= xmax_box   %If the maximum boundary has been reached, chance that it will hop to next cage
        if rand < hop_prob, x(i) = xmax_box + abs(randn)*xstep;    %Hop to next cage
            xmin_box = xmax_box;    %Reset cage position
            xmax_box = xmax_box + box_width;
            hopped = 1;
        else
            x(i) = xmax_box - abs(randn)*xstep;    %Bounce
            bounced = 1;
        end
    else  x(i) = x(i-1) + randn*xstep;      %Otherwise random +/- displacement
    end
    
    if y(i-1) < ymin_box        %If the minimum boundary has been reached, chance that it will hop to next cage
        if rand < hop_prob, y(i) = ymin_box - abs(randn)*ystep;    %Hop to next cage
            ymax_box = ymin_box;
            ymin_box = ymin_box - box_height; %Reset the box position
            hopped = 1;
        else
            y(i) = ymin_box + abs(randn)*ystep; %Bounce
            bounced = 1;
        end
    elseif y(i-1) >= ymax_box    %If the maximum boundary has been reached, chance that it will hop to next cage
        if rand < hop_prob, y(i) = ymax_box + abs(randn)*ystep;    %Hop to next cage
            ymin_box = ymax_box; %Reset box position
            ymax_box = ymax_box + box_height;
            hopped = 1;
        else
            y(i) = ymax_box - abs(randn)*ystep;    %Bounce
            bounced = 1;
        end
    else y(i) = y(i-1) + randn*ystep;        %Otherwise random +/- displacement
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
        delete(boxhandle) %Clear the previous box
        delete(pointhandle) %Clear the previous trace
        % Plot the box
        boxhandle = rectangle('Position',[xmin_box,ymin_box, (xmax_box-xmin_box), (ymax_box-ymin_box)],'edgecolor','r','linewidth',2); %Box outline, x y width height
        hold on
        pointhandle = plot(x,y);
        drawnow
        %     pause(0.05)
        
    end
    
end