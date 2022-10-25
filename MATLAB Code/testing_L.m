% ----------------------------------------------------------------------- %
% testing_L
% ----------------------------------------------------------------------- %
% The purpose of this script is to determine the maximum flow through
% velocity for a given range of L
% ----------------------------------------------------------------------- %
% Establishing the main geometry and channel parameters of the structure
geom_param(1) = 3; % option, geometry chosen
geom_param(2) = 0.1; % D, characteristic length of the geometry
geom_param(3) = 0.04973; % r, outlet radius
geom_param(4) = 1; % N, the number of suns
chan_param(1) = 0.000263; % A, channel width
chan_param(2) = 10*chan_param(1); % B, channel length
chan_param(5) = 0.000263; % S, channel spacing
chan_param(6) = 50*10^-9; % ALD thickness, t
t = chan_param(6); %t
% Establishing the number of channels, X
chan_param(4) = ceil((chan_param(2)-chan_param(5))/(chan_param(1)+chan_param(5))); % X, number of channels in cell

% Calculating the fill factor using the channel parameters
A = chan_param(1);
B = chan_param(2);
X = chan_param(4);
S = chan_param(5);
phi = X*B*A/(X*B*A + S*B*X); % phi, fill factor

% Establishing the altitude vector
altitude = 0:5:80;

% Optimization range
vec = (10*50*10^-9):(1*10^-8):(10^-2);

% Matrix to store results
M = zeros(5,length(vec));

% Performing a local optimization to find the maximum vft for each L
for k = 1:length(vec)

        chan_param(3) = vec(k);
        
        % Calling the calculating force function
        [net_lift,fit,vft,deltaP,deltaT,~,vout,~] = calc_F(altitude,geom_param,chan_param);
        M(1,k) = fit(10);
        M(2,k) = vft(10);
        M(3,k) = deltaP(10);
        M(4,k) = deltaT(10);
        M(5,k) = vout(10);
  
end
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Plotting
% ----------------------------------------------------------------------- %
tiledlayout(2,2)
% Force Sub-Plot
nexttile
plot(vec,M(1,:))
title(['Altitude = ',num2str(altitude),' km, t = 50 nm'])
set(gca, 'XScale', 'log')
xlabel('L (m)')
ylabel('Force (N)')
set(gca,'FontSize',15)

% Vft Sub-Plot
nexttile
plot(vec,M(2,:)) 
title(['Altitude = ',num2str(altitude),' km, t = 50 nm'])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('L (m)')
ylabel('Vft (m/s)')
set(gca,'FontSize',15)

% Pressure Sub-Plot
nexttile
plot(vec,M(3,:))
title(['Altitude = ',num2str(altitude),' km, t = 50 nm'])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('L (m)')
ylabel('\Delta P')
set(gca,'FontSize',15)

% Temperature Sub-Plot
nexttile
plot(vec,M(4,:))
title(['Altitude = ',num2str(altitude),' km, t = 50 nm'])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlabel('L (m)')
ylabel('\Delta T')
set(gca,'FontSize',15)

figure
loglog(vec,M(2,:),'LineWidth', 2)
hold on
loglog(vec,M(5,:),'LineWidth', 2)
xlabel('L (m)')
ylabel('Velocities (m/s)')
legend('Vft','Vout')
set(gca, 'FontSize', 18, 'box', 'off')
% ----------------------------------------------------------------------- %