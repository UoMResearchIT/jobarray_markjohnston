% Brownian motion simulation


%% Generate x and y values for i particles
N = 1000;

tracks = cell(1,100);

for i = 1:length(tracks);
tracks{i}(:,1) = cumsum( randn(N, 1) );
tracks{i}(:,2) = cumsum( randn(N, 1) );
% figure
% plot(tracks{i}(:,1), tracks{i}(:,2));
% ylabel('Y Position');
% xlabel('X Position');
% title(['Track ' num2str(i) ' position versus time in 2D']);
end


%% Compute the Displacement Squared

% The displacement squared is equal to the x coordinate squared plus 
% the y coordinate squared. Since the simulated particle always start 
% at (0,0), it is unnecessary to subtract off the initial position 

for i = 1:length(tracks)
dsquared = tracks{i}(:,1).^ 2 + tracks{i}(:,2).^ 2;
% figure; plot(dsquared);
end

%% Theoretical value of D

% The theoretical value of the diffusion coefficient, D, 
% is given by  where T = temperature (Kelvin), kB = Boltzmann's 
% constant, eta = viscosity, and d = radius.

% Note that the units of D are length squared divided by time. Let's compute D for a 1 micron particle in water at 293 degrees Kelvin.

d    = 1.0e-6;              % radius in meters
eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
kB   = 1.38e-23;            % Boltzmann constant
T    = 293;                 % Temperature in degrees Kelvin

D    = 1;%kB * T / (3 * pi * eta * d);

%% A more realistic particle

dimensions = 2;         % two dimensional simulation
tau = .1;               % time interval in seconds
time = tau * 1:N;       % create a time vector for plotting

k = sqrt(D * dimensions * tau);
dx = k * randn(N,1);
dy = k * randn(N,1);

x = cumsum(dx);
y = cumsum(dy);

dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
 squaredDisplacement = ( x .^ 2) + ( y .^ 2);

plot(x,y);
title('Particle Track of a Single Simulated Particle');

%% Displacement squared plot

clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3);      % plot theoretical line

plot(time, squaredDisplacement);
hold off;
xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time for 1 Particle in 2 Dimensions');

%% Estimated D for simulated data

simulatedD = mean( dSquaredDisplacement ) / ( 2 * dimensions * tau );

%% Uncertainty in the Estimate

%The likely error of this measurement decreases as the square root of the number of samples. This will be discussed in more detail later.

standardError = std( dSquaredDisplacement ) / ( 2 * dimensions * tau * sqrt(N) );
actualError = D - simulatedD;

%% Adding bulk flow to the model

dx = dx + 0.2 * k;
dy = dy + 0.05 * k;

x = cumsum(dx);
y = cumsum(dy);

dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
 squaredDisplacement = ( x .^ 2) + ( y .^ 2);

simulatedD    = mean( dSquaredDisplacement ) / ( 2 * dimensions * tau );
standardError = std(  dSquaredDisplacement ) / ( 2 * dimensions * tau * sqrt(N) );
actualError = D - simulatedD;

plot(x,y);
title('Particle Track of a Single Simulated Particle with Bulk Flow');

%% Displacement Squared in the Presence of Bulk Flow

% Notice how the plot of displacement squared diverges from the theoretical 
% value. It has a distinct quadratic component. The magnitude of this error 
% increases dramatically with time. This suggests that the error caused by 
% bulk flow can be minimized by using the shortest possible sampling period. 
% But there's a catch. As you increase the sampling rate, the amount of 
% noise from the motion tracking algorithm goes up. 

clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3);      % plot theoretical line
plot(time, squaredDisplacement);
hold off;

xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time with Bulk Flow');

%% Multiple particles

tracksCount = 10;
N = 50;
tau = .1;
time = 0:tau:(N-1) * tau;
tracks = { };             % create an empty cell array to hold the results

for i = 1:tracksCount
    tracks{i} = struct();
    tracks{i}.dx = randn(1,N) + 50;
    tracks{i}.x = cumsum(tracks{i}.dx);
    tracks{i}.dy = randn(1,N);
    tracks{i}.y = cumsum(tracks{i}.dy);
    tracks{i}.drsquared = tracks{i}.dx .^2 + tracks{i}.dy .^ 2;
    tracks{i}.rsquared = tracks{i}.x .^ 2 + tracks{i}.y .^ 2;
    tracks{i}.D = mean( tracks{i}.drsquared ) / ( 2 * dimensions * tau );
    tracks{i}.standardError = std( tracks{i}.drsquared ) / ( 2 * dimensions * tau * sqrt(N) );
end

%% Plot the tracks

clf;
hold on;
for i = 1:tracksCount
    plot(tracks{i}.x, tracks{i}.y, 'color', rand(1,3));
end

xlabel('X position (m)');
ylabel('Y position (m)');
title('Combined Particle Tracks');
hold off;

%% Multiple particles with bulk flow

tracksCount = 10;
N = 50;
tau = .1;
time = 0:tau:(N-1) * tau;
tracks = { };             % create an empty cell array to hold the results

for i = 1:tracksCount
    tracks{i} = struct();
    tracks{i}.dx = k * randn(1,N);
    tracks{i}.x = cumsum(tracks{i}.dx);
    tracks{i}.dy = k * randn(1,N);
    tracks{i}.y = cumsum(tracks{i}.dy);
    tracks{i}.drsquared = tracks{i}.dx .^2 + tracks{i}.dy .^ 2;
    tracks{i}.rsquared = tracks{i}.x .^ 2 + tracks{i}.y .^ 2;
    tracks{i}.D = mean( tracks{i}.drsquared ) / ( 2 * dimensions * tau );
    tracks{i}.standardError = std( tracks{i}.drsquared ) / ( 2 * dimensions * tau * sqrt(N) );
end

%% Plot the tracks with bulk flow

clf;
hold on;
for i = 1:tracksCount
    plot(tracks{i}.x, tracks{i}.y, 'color', rand(1,3));
end

xlabel('X position (m)');
ylabel('Y position (m)');
title('Combined Particle Tracks');
hold off;


%% Cross correlation of two

clf;
c = xcorr(tracks{1}.dx, tracks{2}.dx, 'coeff');
xaxis = (1-length(c))/2:1:(length(c)-1)/2;
plot(xaxis, c);
xlabel('Lag');
ylabel('Correlation');
title('Particle 1, 2 x-axis Displacement Cross Correlation');

%% Cross correlation of all particles

% create an array whose columns contain the dx values for each particle

for i = 1:tracksCount
    allDx(:,i) = tracks{i}.dx';
end

% compute all possible auto and cross correlations

c = xcorr(allDx, 'coeff');

% plot the results
clf;
hold on;
for i=1:size(c,1)
    plot(xaxis, c(:,i),'color',rand(1,3));
end
hold off;

xlabel('Lag');
ylabel('Correlation Coefficient');
title('All Possible Auto and Cross Correlations in the x Dimension');
