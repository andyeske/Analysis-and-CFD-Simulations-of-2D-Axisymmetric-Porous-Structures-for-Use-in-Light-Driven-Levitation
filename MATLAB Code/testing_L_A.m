% ----------------------------------------------------------------------- %
% testing_L_A
% ----------------------------------------------------------------------- %
% The purpose of this script is to determine the maximum flow through
% velocity for a given range of A and L
% ----------------------------------------------------------------------- %
% Establishing the main geometry and channel parameters of the structure
geom_param(1) = 3; % option, geometry chosen
geom_param(2) = 0.01; % Ra, characteristic radius
geom_param(3) = 0.005; % Outlet radius
geom_param(4) = 5; % N, the number of suns
chan_param(6) = 50*10^-9; % ALD thickness, t
t = chan_param(6); %t

% Establishing the altitude vector
altitude = 0;

% Optimization range
L_vec = logspace(log10(1*10^-6),log10(10^-2),100);
A_vec = logspace(log10(1*10^-8),log10(500*10^-5),100);

% Matrices to store results
P_mat = zeros(length(A_vec),length(L_vec));
T_mat = zeros(length(A_vec),length(L_vec));
v_mat = zeros(length(A_vec),length(L_vec));
fit_mat = zeros(length(A_vec),length(L_vec));

% Performing a local optimization to find the maximum vft for each A and L
for q = 1:length(A_vec)
    for k = 1:length(L_vec)

        chan_param(3) = L_vec(k);
        chan_param(1) = A_vec(q); 
        
        chan_param(2) = 10*chan_param(1); % B, channel length
        chan_param(5) = chan_param(1); % S, channel spacing
        
        chan_param(4) = ceil((chan_param(2)-chan_param(5))/(chan_param(1)+chan_param(5))); % X, number of channels in cell

        % Calculating the fill factor using the channel parameters
        A = chan_param(1);
        B = chan_param(2);
        X = chan_param(4);
        S = chan_param(5);
        phi = X*B*A/(X*B*A + S*B*X); % phi, fill factor

        % Calling the calculating force function
        [net_lift,fit,vft,deltaP,deltaT,~,~,~] = calc_F(altitude,geom_param,chan_param);
        fit_mat(q,k) = fit;
        v_mat(q,k) = vft;
        P_mat(q,k) = deltaP;
        T_mat(q,k) = deltaT;

    end
end
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Plotting
% ----------------------------------------------------------------------- %
tiledlayout(1,3)
% Vft Sub-Plot
nexttile
s1 = surf(L_vec,A_vec,v_mat,'FaceAlpha',0.5);
c = colorbar;
s1.EdgeColor = 'none';
title({'Vft',['Altitude = ',num2str(altitude),' km, t = ',num2str(t*10^9),' nm'],'Sphere (R = 1cm, Rout = 0.5cm)'})
set(gca, 'XScale', 'log')
ylabel('A (m)')
xlabel('L (m)')
zlabel(' V_{ft} (m/s)')
xticks([10^-8 10^-6 10^-4 10^-2])
yticks([10^-5 10^-4 10^-3 10^-2])
set(gca,'FontSize',15)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
c.Title.String = 'm/s';

% Pressure Sub-Plot
nexttile
s = surf(L_vec,A_vec,P_mat,'FaceAlpha',0.5);
c = colorbar;
s.EdgeColor = 'none';
title({'Pressure',['Altitude = ',num2str(altitude),' km, t = ',num2str(t*10^9),' nm'],'Sphere (R = 1cm, Rout = 0.5cm)'})
set(gca, 'XScale', 'log')
ylabel('A (m)')
xlabel('L (m)')
xticks([10^-8 10^-6 10^-4 10^-2])
yticks([10^-5 10^-4 10^-3 10^-2])
zlabel('Pressure (Pa)')
set(gca,'FontSize',15)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
c.Title.String = 'Pa';

% Temperature Sub-Plot
nexttile
s = surf(L_vec,A_vec,T_mat,'FaceAlpha',0.5);
c = colorbar;
s.EdgeColor = 'none';
title({'Temperature',['Altitude = ',num2str(altitude),' km, t = ',num2str(t*10^9),' nm'],'Sphere (R = 1cm, Rout = 0.5cm)'})
set(gca, 'XScale', 'log')
ylabel('A (m)')
xlabel('L (m)')
xticks([10^-8 10^-6 10^-4 10^-2])
yticks([10^-5 10^-4 10^-3 10^-2])
zlabel('\DeltaT')
set(gca,'FontSize',15)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
c.Title.String = 'K';
% ----------------------------------------------------------------------- %