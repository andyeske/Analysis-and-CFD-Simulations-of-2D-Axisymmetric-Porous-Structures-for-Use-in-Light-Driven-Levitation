% ----------------------------------------------------------------------- %
% Computing the forces 
% ----------------------------------------------------------------------- %
% The purpose of this script is to determine the maximum flow through
% velocity for a given range of A and t
% ----------------------------------------------------------------------- %
% Establishing the main geometry and channel parameters of the structure
geom_param(1) = 3; % option, geometry chosen
geom_param(2) = 0.01; % Ra, characteristic radius
geom_param(3) = 0.005; % q, outlet sphere radius or cone/rocket length
geom_param(4) = 5; % N, the number of suns
chan_param(3) = 50*10^-6; % L

% Establishing the altitude vector
altitude = 0;

% Optimization range
vec_t = logspace(log10(1*10^-9),log10(100*10^-9),100);
vec_A = logspace(log10(1*10^-7),log10(50*10^-5),100);

% Matrices to store results
P_mat = zeros(length(vec_A),length(vec_t));
T_mat = zeros(length(vec_A),length(vec_t));
v_mat = zeros(length(vec_A),length(vec_t));
fit_mat = zeros(length(vec_A),length(vec_t));

% Performing a local optimization to find the maximum vft for each t and A
for q = 1:length(vec_A)
    for k = 1:length(vec_t)

        % Calculating the fill factor using the channel parameters
        chan_param(6) = vec_t(k);
        chan_param(1) = vec_A(q); % A
        chan_param(2) = 10*chan_param(1); % B, channel length
        chan_param(5) = chan_param(1); % S
        chan_param(4) = ceil((chan_param(2)-chan_param(5))/(chan_param(1)+chan_param(5))); % X, number of channels in cell
        
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
% Pressure Sub-Plot
nexttile
s = surf(vec_t,vec_A,P_mat,'FaceAlpha',0.5);
c = colorbar;
s.EdgeColor = 'none';
title({'Pressure',['Altitude = ',num2str(altitude),' km'], 'Sphere (R = 1cm, Rout = 0.5cm)'})
set(gca, 'XScale', 'log')
ylabel('A (m)')
xlabel('t (m)')
zlabel('\DeltaP')
xticks([10^-8 10^-6 10^-4 10^-2])
yticks([10^-5 10^-4 10^-3 10^-2])
set(gca,'FontSize',15)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
c.Title.String = 'Pa';

% Vft Sub-Plot
nexttile
s1 = surf(vec_t,vec_A,v_mat,'FaceAlpha',0.5);
c = colorbar;
s1.EdgeColor = 'none';
title({'Vft',['Altitude = ',num2str(altitude),' km'], 'Sphere (R = 1cm, Rout = 0.5cm)'})
set(gca, 'XScale', 'log')
ylabel('A (m)')
xlabel('t (m)')
xticks([10^-9 10^-8 10^-7])
yticks([10^-5 10^-4 10^-3 10^-2])
zlabel('V_{ft}')
set(gca,'FontSize',15)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
c.Title.String = 'm/s';

% Temperature Sub-Plot
nexttile
s = surf(vec_t,vec_A,T_mat,'FaceAlpha',0.5);
c = colorbar;
s.EdgeColor = 'none';
title({'Temperature',['Altitude = ',num2str(altitude),' km'], 'Sphere (R = 1cm, Rout = 0.5cm)'})
set(gca, 'XScale', 'log')
ylabel('A (m)')
xlabel('t (m)')
xticks([10^-8 10^-6 10^-4 10^-2])
yticks([10^-5 10^-4 10^-3 10^-2])
zlabel('\DeltaT')
set(gca,'FontSize',15)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
c.Title.String = 'K';
% ----------------------------------------------------------------------- %