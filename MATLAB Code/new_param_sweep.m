% ----------------------------------------------------------------------- %
% Parametric sweep
% ----------------------------------------------------------------------- %
% The purpose of this code is to consider a range for the various variables
% and find the most optimal combination of them. This would be the
% combination that yields the highest lift at the lowest altitude
% ----------------------------------------------------------------------- %
% Recall the inputs of the calc_F function. These are as follows:
% geom_param(1): option, geometry chosen
% geom_param(2): Ra, characteristic radius
% geom_param(3): l, length of cone (option 2), outlet radius of the sphere 
% (option 3), or length of the rocket (option 4)
% geom_param(4): N, number of suns
% chan_param(1): A, channel width
% chan_param(2): B, channel length
% chan_param(3): L, channel thickness
% chan_param(4): X, number of channels in cell
% chan_param(5): S, channel spacing
% chan_param(6): t, ALD thickness
% ----------------------------------------------------------------------- %
% Establishing the altitude vector
altitude = 0:5:80;

geom_param(1) = 3; % Geometry option
geom_param(2) = 0.01; % Geometry characteristic radius
geom_param(4) = 5; % N, the nummber of suns
chan_param(6) = 50*10^-9; %t, ALD thickness
        

% Variables to vary:
% In essence, we are establishing vectors that contain the different
% combinations of variables we want to vary, including the A lengths,
% the channel thickness, and the channel spacing

% This is the vector that contains the channel widths A
% We are only modifying A, and setting B to always be 10 times A
min_A = 1*10^-8;
max_A = 500*10^-5;
num_A = 100; % number of points in the vector
A_vector = logspace(log10(min_A),log10(max_A),num_A);

% Channel Thickness L
min_L = 1*10^-6;
max_L = 10^-2;
num_L = 100; % number of points in the vector
L = logspace(log10(min_L),log10(max_L),num_L);

% Outlet Radius
min_ra = 10^-4;
max_ra = 10^-2;
num_ra = 100; % number of points in the vector
ra = logspace(log10(min_ra),log10(max_ra),num_ra);

% Creating the matrices that will store all of the solutions
% For all the the particular combinations, these first two matrices
% store the minimum altitude at which the net payload became positive,
% as well as the corresponding payload
minimum_altitude = zeros(num_A,num_L,num_ra);
minimum_altitude_payload = zeros(num_A,num_L,num_ra);
% These two other matrices store the maximum payload recorded in each
% simulation, and the altitude corresponding to that payload. The last matrix
% stores the corresponding computed aerial densities for each simulation
maximum_payload = zeros(num_A,num_L,num_ra);
maximum_payload_altitude = zeros(num_A,num_L,num_ra);
maximum_payload_aerial = zeros(num_A,num_L,num_ra);

% Next, we will create a triple for-loop, where we will iterate for each
% one of these parameter combinations, and calculate the force and payload
% capabilities for each
for qq = 1:num_L
    for q = 1:num_A
        for qqq = 1:num_ra

            % Fixing the magnitude of the channel width A
            chan_param(1) = A_vector(q);
            % We obtain the channel length B by multiplying A by 10
            chan_param(2) = 10*A_vector(q);
            % We set S = A
            chan_param(5) = A_vector(q);
            % We write the channel thickness L
            chan_param(3) = L(qq);
            % Number of channels in a unit cell
            chan_param(4) = ceil((10*chan_param(1)-chan_param(5))/(chan_param(1)+chan_param(5))); 
            % Outlet radius
            geom_param(3) = ra(qqq);
           
            % If statement that ensures we have at least 1 channel in the
            % cell
            if chan_param(4) > 0

                % Then, we call calc_F, which will take this variable
                % combination and output two vectors: the net_lift 
                % generated and the force corresponding to that
                [net_lift,fit,~,~,~,~,~,aerial] = calc_F(altitude,geom_param,chan_param);

                % Next, we can use this information to find the first point
                % at which the payload becomes positive, together with its 
                % corresponding altitude. We can also find for the same
                % simulation the maximum net-payload, together with its
                % altitude
                if sum(net_lift > 0) > 0
                    % First, we isolate the cases where we actually have a
                    % positive payload
                    positive_payload = net_lift(net_lift > 0);
                    pos_payload_alt = altitude(net_lift > 0);
                    % In the matrix, we only store the first of these
                    % positive payloads
                    minimum_altitude_payload(q,qq,qqq) = positive_payload(1);
                    % Together with the corresponding payload
                    minimum_altitude(q,qq,qqq) = pos_payload_alt(1);
                    % Now, we can find the maximum payload and its 
                    % corresponding altitude
                    maximum_payload(q,qq,qqq) = net_lift(net_lift == max(net_lift));
                    maximum_payload_altitude(q,qq,qqq) = altitude(net_lift == max(net_lift));
                    maximum_payload_aerial(q,qq,qqq) = aerial(net_lift == max(net_lift));
                else
                    % Otherwise, if the returned vector with the
                    % net-payloads has no single positive value to begin
                    % with, then we just set set values of the matrices to
                    % be trivial
                    minimum_altitude_payload(q,qq,qqq) = 0;
                    minimum_altitude(q,qq,qqq) = 100;
                    maximum_payload(q,qq,qqq) = 0;
                    maximum_payload_altitude(q,qq,qqq) = 100;
                    maximum_payload_aerial(q,qq,qqq) = 0;
                end
            else
                minimum_altitude_payload(q,qq,qqq) = 0;
                minimum_altitude(q,qq,qqq) = 100;
                maximum_payload(q,qq,qqq) = 0;
                maximum_payload_altitude(q,qq,qqq) = 100;
                maximum_payload_aerial(q,qq,qqq) = 0;
            end

        end
    end
end

% Now that the triple for loop is done, we can begin looking within
% these matrices for the cases that we are interested in. In
% particular, we want to obtain from "minimum_altitude" the minimum
% altitude (at which there was a positive payload), and the max payload 
% corresponding to it. From "maximum_payload", we want to obtain the 
% maximum payload generated, and the minimum altitude corresponding to 
% it.

% First, we start by finding the ultimate minimum altitude and the 
% ultimate maximum payload 
ultimate_minimum_altitude = min(min(min(minimum_altitude(:,:,:))));
ultimate_maximum_payload = max(max(max(maximum_payload(:,:,:))));

% Next, we continue by looking at the coordinates of the payload 
% corresponding to the ultimate minimum altitude 
vec = find(minimum_altitude(:,:,:) ==  ultimate_minimum_altitude);
% Setting the corresponding minimum payload
minimum_payload = max(minimum_altitude_payload(vec));
max1 = find(minimum_altitude_payload(vec) == max(minimum_altitude_payload(vec)));
% Finding the optimization parameters from above (BA ratio, L and S)
% corresponding to the minimum payload just found
[c1,c2,c3] = ind2sub(size(minimum_altitude_payload),vec(max1));

% Now, we can repeat this same process, but by looking at the other
% matrix. We want to find the coordinate of the altitude corresponding 
% to the ultimate maximum payload
vec2 = find(maximum_payload(:,:,:) == ultimate_maximum_payload);
% Setting the maximum altitude
maximum_altitude = min(maximum_payload_altitude(vec2));
max2 = find(maximum_payload_altitude(vec2) == min(maximum_payload_altitude(vec2)));
% Finding the specific area parameters
[cc1,cc2,cc3] = ind2sub(size(maximum_payload),vec2(max2));

% Finally, after having performed all of this procedure, we are able to
% display a message indicating the results
disp('-----------------------------------------------------------------------------------------------------')
disp(['For the given constraints, the minimum altitude at which the net payload becomes positive is ',num2str(ultimate_minimum_altitude),' km'])
disp(['and has the potential of carrying ',num2str((10^6)*minimum_payload),' mg as its maximum payload'])
disp(['The corresponding parameters that yielded these results are A = ',num2str(min(A_vector(c1))),' m, L = ',num2str(max(L(c2))),' m, and Ra_{out} = ',num2str(min(ra(c3))),' m'])
disp('-----------------------------------------------------------------------------------------------------')
disp(['Additionally, the maximum payload found was ',num2str((10^6)*ultimate_maximum_payload),' mg'])
disp(['and corresponded to an altitude of ',num2str(maximum_altitude),' km'])
disp(['The corresponding parameters that yielded these results are A = ',num2str(min(A_vector(cc1))),' m, L = ',num2str(max(L(cc2))),' m, and Ra_{out} = ',num2str(min(ra(cc3))),' m'])
disp('-----------------------------------------------------------------------------------------------------')
disp(['Finally, the aerial density corresponding to this maximum payload capability was ',num2str(1000*maximum_payload_aerial(cc1,cc2,cc3)),' g/m2'])
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Plotting
% ----------------------------------------------------------------------- %
% Making the 3D Cloud Plot for the Minimum Altitudes 
surf0 = isosurface(L, A_vector, ra, minimum_altitude, ultimate_minimum_altitude);
p0 = patch(surf0);
isonormals(L, A_vector, ra, minimum_altitude, p0);
set(p0,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.2);

surf1 = isosurface(L, A_vector, ra, minimum_altitude, ultimate_minimum_altitude+5);
p1 = patch(surf1);
isonormals(L, A_vector, ra, minimum_altitude, p1);
set(p1,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2);

surf2 = isosurface(L, A_vector, ra, minimum_altitude, ultimate_minimum_altitude+10);
p2 = patch(surf2);
isonormals(L, A_vector, ra, minimum_altitude, p2);
set(p2,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.2);

surf3 = isosurface(L, A_vector, ra, minimum_altitude, ultimate_minimum_altitude+15);
p3 = patch(surf3);
isonormals(L, A_vector, ra, minimum_altitude, p3);
set(p3,'FaceColor','red','EdgeColor','none','FaceAlpha',0.2);

view(3); 
camlight; lighting gouraud
xlabel('{\it L} [m]')
ylabel('{\it A} [m]')
zlabel('{\it l} [m]')
Leg = legend([num2str(ultimate_minimum_altitude),' km'],[num2str(ultimate_minimum_altitude+5),' km'],[num2str(ultimate_minimum_altitude+10),' km'],[num2str(ultimate_minimum_altitude+15),' km']);
Leg.Location = 'eastoutside';
if geom_param(1) == 1
    title({'3D cloud plot for various altitudes','Arbitrary Geometry'})
elseif geom_param(1) == 2
    %title({'3D cloud plot for various altitudes','Cone Geometry'})
    title(Leg,'Minimum Altitudes: Cone Geometry');
elseif geom_param(1) == 3
    %title({'3D cloud plot for various altitudes','Sphere Geometry'}) 
    title(Leg,'Minimum Altitudes: Sphere Geometry');
else
    %title({'3D cloud plot for various altitudes','Rocket Geometry'})
    title(Leg,'Minimum Altitudes: Rocket Geometry');
end
set(gca,'FontSize',18)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
set(gca,'FontName','Times New Roman')


% Making the 3D Cloud Plot for the Maximum Payloads
figure
surf0 = isosurface(L, A_vector, ra, (10^6)*maximum_payload, (10^6)*0.99*ultimate_maximum_payload);
p0 = patch(surf0);
isonormals(L, A_vector, ra, (10^6)*maximum_payload, p0);
set(p0,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.2);

surf1 = isosurface(L, A_vector, ra, (10^6)*maximum_payload, (10^6)*0.9*ultimate_maximum_payload);
p1 = patch(surf1);
isonormals(L, A_vector, ra, (10^6)*maximum_payload, p1);
set(p1,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2);

surf2 = isosurface(L, A_vector, ra, (10^6)*maximum_payload, (10^6)*0.5*ultimate_maximum_payload);
p2 = patch(surf2);
isonormals(L, A_vector, ra, (10^6)*maximum_payload, p2);
set(p2,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.2);

surf3 = isosurface(L, A_vector, ra, (10^6)*maximum_payload, (10^6)*0.3*ultimate_maximum_payload);
p3 = patch(surf3);
isonormals(L, A_vector, ra, (10^6)*maximum_payload, p3);
set(p3,'FaceColor','red','EdgeColor','none','FaceAlpha',0.2);

view(3); 
camlight; lighting gouraud
xlabel('{\it L} [m]')
ylabel('{\it A} [m]')
zlabel('{\it l} [m]')
Leg = legend([num2str(0.99*(10^6)*ultimate_maximum_payload),' mg',' | 99% max. payload'],[num2str(0.9*(10^6)*ultimate_maximum_payload),' mg',' | 90% max. payload'],[num2str((10^6)*ultimate_maximum_payload/2),' mg',' | 50% max. payload'],[num2str((10^6)*0.3*ultimate_maximum_payload),' mg',' | 30% max. payload']);
Leg.Location = 'eastoutside';
if geom_param(1) == 1
    title({'3D cloud plot for various payloads','Arbitrary Geometry'})
elseif geom_param(1) == 2
    %title({'3D cloud plot for various payloads','Cone Geometry'})
    title(Leg,'Maximum Payloads: Cone Geometry');
elseif geom_param(1) == 3
    %title({'3D cloud plot for various payloads','Sphere Geometry'})
    title(Leg,'Maximum Payloads: Sphere Geometry');
else
   % title({'3D cloud plot for various payloads','Rocket Geometry'})
    title(Leg,'Maximum Payloads: Rocket Geometry');
end
set(gca,'FontSize',18)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
set(gca,'FontName','Times New Roman')


% Making the 3D Cloud Plot for the Aerial Densities
figure
VEC = reshape(maximum_payload_aerial,1,{});
surf0 = isosurface(L, A_vector, ra, maximum_payload_aerial, prctile(VEC(VEC~=0),99,'all'));
p0 = patch(surf0);
isonormals(L, A_vector, ra, minimum_altitude, p0);
set(p0,'FaceColor','cyan','EdgeColor','none','FaceAlpha',0.2);

surf1 = isosurface(L, A_vector, ra, maximum_payload_aerial, prctile(VEC(VEC~=0),90,'all'));
p1 = patch(surf1);
isonormals(L, A_vector, ra, minimum_altitude, p1);
set(p1,'FaceColor','yellow','EdgeColor','none','FaceAlpha',0.2);

surf2 = isosurface(L, A_vector, ra, maximum_payload_aerial, prctile(VEC(VEC~=0),50,'all'));
p2 = patch(surf2);
isonormals(L, A_vector, ra, minimum_altitude, p2);
set(p2,'FaceColor','blue','EdgeColor','none','FaceAlpha',0.2);

surf3 = isosurface(L, A_vector, ra, maximum_payload_aerial, prctile(VEC(VEC~=0),25,'all'));
p3 = patch(surf3);
isonormals(L, A_vector, ra, minimum_altitude, p3);
set(p3,'FaceColor','red','EdgeColor','none','FaceAlpha',0.2);

view(3); 
camlight; lighting gouraud
xlabel('{\it L} [m]')
ylabel('{\it A} [m]')
zlabel('{\it l} [m]')
Leg = legend([num2str(1000*prctile(VEC(VEC~=0),99,'all')),' g/m^{2}',' | 99% percentile density'],[num2str(1000*prctile(VEC(VEC~=0),90,'all')),' g/m^{2}',' | 90% percentile density'],[num2str(1000*prctile(VEC(VEC~=0),50,'all')),' g/m^{2}',' | 50% percentile density'],[num2str(1000*prctile(VEC(VEC~=0),25,'all')),' g/m^{2}',' | 25% percentile density']);
Leg.Location = 'eastoutside';
if geom_param(1) == 1
    title({'3D cloud plot for various altitudes','Arbitrary Geometry'})
elseif geom_param(1) == 2
    %title({'3D cloud plot for various altitudes','Cone Geometry'})
    title(Leg,'Aerial Densities: Cone Geometry');
elseif geom_param(1) == 3
    %title({'3D cloud plot for various altitudes','Sphere Geometry'}) 
    title(Leg,'Aerial Densities: Sphere Geometry');
else
    %title({'3D cloud plot for various altitudes','Rocket Geometry'})
    title(Leg,'Aerial Densities: Rocket Geometry');
end
set(gca,'FontSize',18)
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'zscale','log')
set(gca,'FontName','Times New Roman')
% ------------------------------------------------------------------- %