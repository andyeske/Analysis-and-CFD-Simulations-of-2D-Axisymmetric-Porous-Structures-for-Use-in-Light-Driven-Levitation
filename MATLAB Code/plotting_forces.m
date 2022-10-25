% ----------------------------------------------------------------------- %
% Plotting Forces
% ----------------------------------------------------------------------- %
% The purpose of this code is to plot the results from the calc_F function
% ----------------------------------------------------------------------- %
% Establishing the main geometry and channel parameters of the structure
geom_param(1) = 2; % option, geometry chosen
geom_param(2) = 0.1; % D, characteristic length of the geometry
geom_param(3) = 0.0202; % r, outlet radius
geom_param(4) = 1; % N, the number of suns
chan_param(1) = 0.000263; % A, channel width
chan_param(2) = 10*chan_param(1); % B, channel length
chan_param(3) = 2.121*10^-6; % L, channel thickness
chan_param(5) = chan_param(1); % S, channel spacing
chan_param(4) = ceil((chan_param(2)-chan_param(5))/(chan_param(1)+chan_param(5))); % X, number of channels in cell
chan_param(6) = 50*10^-9; % t, the ALD thickness

% Calculating the fill factor using the channel parameters
A = chan_param(1);
B = chan_param(2);
X = chan_param(4);
S = chan_param(5);
t = chan_param(6);
phi = X*B*A/(X*B*A + S*B*X); % phi, fill factor

% Establishing the altitude vector
altitude = 0:5:80;

% Calling the calculating force function
 [net_lift,fit,vft,deltaP,deltaT,~,~,aerial] = calc_F(altitude,geom_param,chan_param);
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Plotting
% ----------------------------------------------------------------------- %
% Force against Altitude Plot
% The following plot describes how force changes as a function of altitude
tiledlayout(1,3)
nexttile
semilogy(altitude, fit, '-r', 'LineWidth', 2)
ylabel('Force (N)');
xlabel('Altitude (km)');
if geom_param(1) == 1
    title({'Reaction forces for various altitudes','Arbitrary Geometry'})
elseif geom_param(1) == 2
    title({'Reaction forces for various altitudes','Cone Geometry'})
elseif geom_param(1) == 3
    title({'Reaction forces for various altitudes','Sphere Geometry'})     
else
    title({'Reaction forces for various altitudes','Rocket Geometry'})
end
l = legend(['Maximum force = ' num2str(max(fit)) 'N, at ' num2str(altitude(find(fit == max(fit)))) 'km'],'EdgeColor','none','Location','Southeast');
set(gca,'XTick',0:20:80)
set(gca, 'FontSize', 18, 'box', 'off')
l.FontSize = 12;
% ----------------------------------------------------------------------- %
% Net Lift Plot
nexttile
% The following plot describes how lift changes as a function of altitude
loglog(altitude,abs(net_lift),'LineWidth',2)
% Finding the first positive net_lift and its corresponding altitude
pos_lift = net_lift(net_lift > 0);
pos_lift_alt = altitude(net_lift > 0);
first_lift = pos_lift(1);
first_lift_alt = pos_lift_alt(1);
% Finding the maximum payload and its corresponding altitude
max_lift = net_lift(net_lift == max(net_lift));
max_lift_alt = altitude(net_lift == max(net_lift));
hold on
loglog(max_lift_alt,max_lift,'o','LineWidth',2)
ylabel('Net Lift (N)')
xlabel('Altitude (km)')
set(gca,'XTick',0:20:80)
if geom_param(1) == 1
    title({'Payload for various altitudes','Arbitrary Geometry',['\phi (fill factor) = ' num2str(phi)]})
elseif geom_param(1) == 2
    title({'Payload for various altitudes','Cone Geometry',['\phi (fill factor) = ' num2str(phi)]})
elseif geom_param(1) == 3
    title({'Payload for various altitudes','Sphere Geometry',['\phi (fill factor) = ' num2str(phi)]})     
else
    title({'Payload for various altitudes','Rocket Geometry',['\phi (fill factor) = ' num2str(phi)]})
end
hold on
logic = find(net_lift > 0);
pa = patch([altitude(logic(1)-1);altitude(logic(1)-1);pos_lift_alt'], [pos_lift(end);abs(net_lift(logic(1)-1));pos_lift], 'y');
pa.FaceAlpha = 0.05;
l = legend('Structure Net Lift',['Maximum Net Lift (',...
    num2str(max_lift), ' N at ', num2str(max_lift_alt),...
    ' km)'],'Positive Payload Zone','EdgeColor','none','Location','Southwest');
set(gca, 'FontSize', 18, 'box', 'off')
l.FontSize = 12;
% ----------------------------------------------------------------------- %
% Vft Plot
nexttile
semilogy(altitude, vft, '-b', 'LineWidth', 2)
ylabel('V_{ft} (m/s)');
xlabel('Altitude (km)');
if geom_param(1) == 1
    title({'V_{ft} for various altitudes','Arbitrary Geometry'})
elseif geom_param(1) == 2
    title({'V_{ft} for various altitudes','Cone Geometry'})
elseif geom_param(1) == 3
    title({'V_{ft} for various altitudes','Sphere Geometry'})     
else
    title({'V_{ft} for various altitudes','Rocket Geometry'})
end
set(gca,'XTick',0:20:80)
l = legend(['Maximum V_{ft} = ' num2str(max(vft)) 'm/s, at ' num2str(altitude(find(vft == max(vft)))) 'km'],'EdgeColor','none','Location','SouthEast');
set(gca, 'FontSize', 18, 'box', 'off')
l.FontSize = 12;
% ----------------------------------------------------------------------- %