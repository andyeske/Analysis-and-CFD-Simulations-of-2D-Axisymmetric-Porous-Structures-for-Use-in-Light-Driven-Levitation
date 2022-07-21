% ----------------------------------------------------------------------- %
% testing_Isun
% ----------------------------------------------------------------------- %
% The purpose of this script is to calculate the dependency of the maximum
% flow-through velocity as a function of the sun intensity
% ----------------------------------------------------------------------- %
% Establishing the main geometry and channel parameters of the structure
geom_param(1) = 3; % option, geometry chosen
geom_param(2) = 0.01; % Ra, characteristic radius
geom_param(3) = 0.005; % Outlet radius l
chan_param(6) = 0*10^-9; %ALD thickness t
t = chan_param(6); %t

% Establishing the altitude vector
altitude = 20;

% Optimization range
L_vec = logspace(log10(1*10^-6),log10(10^-2),100); % for channel thickness L
A_vec = logspace(log10(1*10^-8),log10(500*10^-5),100); % for channel width A
I_vec = logspace(log10(1),log10(100),50); % for number of suns N

% Matrices to store results
v_mat = zeros(length(A_vec),length(L_vec));
v_mat2 = zeros(length(A_vec),length(L_vec));
v_mat3 = zeros(length(A_vec),length(L_vec));
v_mat_max = zeros(5,length(I_vec));
% ----------------------------------------------------------------------- %

% Iterating through each sun intensity
for h = 1:length(I_vec)
    
    Isun = I_vec(h);
    geom_param(4) = Isun;
    
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
            [net_lift,fit,vft,deltaP,deltaT,vft2,vft3,~] = calc_F(altitude,geom_param,chan_param);
            % Storing values in the matrices
            v_mat(q,k) = vft;
            v_mat2(q,k) = vft2;
            v_mat3(q,k) = vft3;
            
        end
    end

% Finding for each sun intensity the maximum vft
v_mat_max(1,h) = max(max(v_mat));
v_mat_max(2,h) = max(max(v_mat2));
v_mat_max(3,h) = max(max(v_mat3));
[v_mat_max(4,h),v_mat_max(5,h)] = find(v_mat(:,:) == max(max(v_mat)));

    
end
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Plotting
% ----------------------------------------------------------------------- %
semilogy(I_vec(:),v_mat_max(1,:),'LineWidth',1)
hold on
semilogy(I_vec(:),v_mat_max(2,:),'LineWidth',1)
hold on
semilogy(I_vec(:),v_mat_max(3,:),'LineWidth',1)
legend('Full Analytical','Simplification','Simplification w/usual \DeltaT','Location','NorthWest')
xlabel('# of Isun')
ylabel('Maximum Vft')
title({'Dependence of Vft on # of Isun',['ALD thickness = ',num2str(t*10^9),' nm']})
set(gca,'FontSize',15)
% ----------------------------------------------------------------------- %