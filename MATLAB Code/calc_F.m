% ----------------------------------------------------------------------- %
% Calculating vft and vout - function
% ----------------------------------------------------------------------- %
% The purpose of this code is to provide a framework for calculating the
% flow-through velocity (vft) and outlet velocity (vout) of an arbitrary 3D
% geometry. 
% ----------------------------------------------------------------------- %

function [net_lift,fit,vft,deltaP,deltaT,vft2,vout,aerial] = calc_F(altitude,geom_param,chan_param)

% option: different geometries available to choose from. These correspond
% to 1 - manual set-up, 2 - cone, 3 - sphere, and 4 - rocket 

% D: characteristic length of the structure (the diameter of the sphere
% and cone, the length of the rocket)
% (units of m)

% r: outlet radius
% (units of m)

% phi: geometric fill factor. This is simply the percentage that the 
% structure's porous area occupies out of the entire surface area
% (unitless)

% A: channel width
% (units of m)

% B: channel length
% (units of m)

% L: channel depth/thickness. 
% (units of m)

% X: number of channels
% (unitless)

% S: spacing between channels
% (units of m)

option = geom_param(1); 
D = geom_param(2); 
r = geom_param(3); 
A = chan_param(1);
B = chan_param(2);
L = chan_param(3);
X = chan_param(4);
S = chan_param(5);

% Calculating the fill factor using the channel parameters
phi = X*B*A/(X*B*A + S*B*X);

% ----------------------------------------------------------------------- %
% Altitude-dependent variables 
% ----------------------------------------------------------------------- %
% In the section below, we proceed to introduce scripts that will store 
% information about air pressure, air temperature, air density and air 
% viscosity at various altitudes, for the earth's atmosphere
% Most of these equations where derived in the "temperature_fit" script
% ----------------------------------------------------------------------- %
% Altitude vector 
% (units in km)
altitude = altitude';
alt = altitude;
% ----------------------------------------------------------------------- %
% Temperature
% (units of K)
p15 = [-4.59195520629687e-29
    4.02322550536528e-27
    1.49047601370574e-23
    -7.94159585089666e-21
    2.02100883624221e-18
    -3.15169237730408e-16
    3.27061799378819e-14
    -2.33217472878003e-12
    1.14960632799406e-10
    -3.86235307877923e-09
    8.52520967699876e-08
    -1.14955512755429e-06
    8.15422681932410e-06
    -2.28247698757490e-05
    9.91066193190049e-05
    0.00347312785053852];
T = (p15(1).*alt.^15 + p15(2).*alt.^14 + p15(3).*alt.^13 + p15(4).*alt.^12 +...
     p15(5).*alt.^11 + p15(6).*alt.^10 + p15(7).*alt.^9 + p15(8).*alt.^8 + p15(9).*alt.^7 +...
     p15(10).*alt.^6 + p15(11).*alt.^5 + p15(12).*alt.^4 + p15(13).*alt.^3 +...
     p15(14).*alt.^2 + p15(15).*alt + p15(16)).^-1;
% ----------------------------------------------------------------------- %
% Pressure
% (units of Pa)
inte = (1/16)*p15(1).*alt.^16 + (1/15)*p15(2).*alt.^15 + (1/14)*p15(3).*alt.^14 + (1/13)*p15(4).*alt.^13 +...
     (1/12)*p15(5).*alt.^12 + (1/11)*p15(6).*alt.^11 + (1/10)*p15(7).*alt.^10 + (1/9)*p15(8).*alt.^9 + (1/8)*p15(9).*alt.^8 +...
     (1/7)*p15(10).*alt.^7 + (1/6)*p15(11).*alt.^6 + (1/5)*p15(12).*alt.^5 + (1/4)*p15(13).*alt.^4 +...
     (1/3)*p15(14).*alt.^3 + (1/2)*p15(15).*alt.^2 + p15(16).*alt;
% R: Ideal gas constant for air
% (units of J/Kg*K)
% g: earth's gravitational constant
% (units of m/s2)
R = 287;
g  = 9.81;
P = 101300*exp(-1000*g*inte/R);
% ----------------------------------------------------------------------- %
% Density
% (units in kg/m3)
% We will use the ideal gas law
rho = P./(R.*T);
% ----------------------------------------------------------------------- %
% Viscosity
% (units in Pa*s)
% We will use Sutherland's law
T_ref = 293;
mu_ref = 0.000018205;
S0 = 110.4;
mu = mu_ref*((T/T_ref).^1.5).*((T_ref+S0)./(T+S0));
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Geometrical properties
% ----------------------------------------------------------------------- %
% In this section, we define the geometrical properties of the 3D structure
% that we are intending to use for analysis. The main parameters that have
% to be established include: 

% A_in: inlet area through which air flows into the 3D structure. In a
% porous geometry with channels in its surface, A_in would be the sum of 
% all of the channel walls (since that's where we define vft)
% (units of m2)

% A_out: outlet area through which the air flows out of the 3D structure.
% In general, this is the cross-sectional area of the structure's nozzle
% (units of m2)

% A_total: total surface area of the 3D structure. 
% (units of m2)

% A_solid: surface area where the sun's irradiance is absorbed, given by
% A_total - A_in.
% (units of m2)
% ----------------------------------------------------------------------- %
% There are 4 options to choose from: 
% | 1 - manual set-up | 2 - cone | 3 - sphere | 4 - rocket |
% ----------------------------------------------------------------------- %
if option == 1
    % Manual Set-up
    Ra = 0.31; 
    A_in = 0.157;
    A_out = 0.0000005;
    A_total = 0.2;
elseif option == 2
    % Cone geometry
    Ra = D/2;
    h1 = sqrt((Ra - r)^2 + D^2);
    h3 = (D^2)/(D - 2*r);
    h2 = sqrt(Ra^2 + h3^2);
    A_out = 3.14*r^2;% outlet area
    A_total = 3.14*Ra^2 + 3.14*Ra*h2 - 3.14*r*(h2 - h1);
    A_in = phi*A_total;
    A_solid = (1-phi)*A_total;
elseif option == 3
    % Sphere geometry
    Ra = D/2;
    h = Ra - sqrt(Ra^2 - r^2); % sphere cap height
    A_out = 3.14*r^2; % outlet area
    A_total = 4*3.14*Ra^2 - 2*Ra*h*3.14; %inlet area
    A_in = phi*A_total;
    A_solid = (1-phi)*A_total;
else
    % Rocket geometry
    A_out = 3.14*r^2; % outlet area
    A_total = 2*3.14*r^2 + 2*3.14*r*D;
    A_in = phi*A_total;
    A_solid = (1-phi)*A_total;
end
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Important physical variables
% ----------------------------------------------------------------------- %
% In the next section, we begin introducing the governing variables used
% throughout the equations that this code implements. These are:

% kb: Boltzmann constant
% (units of J/kg)

% m: mass of a single air molecule. This was obtained by dividing the molar
% mass of air by Avogadro's constant
% (units of kg)

% k_0: conductivity of air at standard atmospheric conditions
% (units of W/m*K)

% I_sun = Intensity of the sun. This is approximated using a linear
% interporlation formula, knowing that at 0km (sea level), the irradiation
% from the sun reaching earth is 1000 W/m2, while at 100km, the 
% corresponding irradiation is 1380 W/m2
% (units of W/m2)

% M = molecular mass of air
% (units of kg/mol)
% ----------------------------------------------------------------------- %
kb = 1.38*10^-23; 
m = 4.8089*10^-26; 
k_0 = 0.024; 
I_sun = 1000+3.8.*altitude;
% Setting the intensity of the sun to be some arbitrary multiple
I_sun = geom_param(4)*I_sun;
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Derived variables
% ----------------------------------------------------------------------- %
% In the following section, we now proceed to derive other variables also
% used throughout the equations that this code implements, using the
% previously introduced variables and scripts. These are: 

% k_air: conductivity of air, but both pressure and temperature dependent.
% This expression is credited to Wu
% (units of W/m*K)

% k_alumina: conductivity of alumina, which was measured experimentally.
% (units of W/m*K)

% P_star: average pressure between the inner side and outer side of the 3D 
% structure. Taking nanocardboard as an example, P_star would be the 
% average between the pressure on the plate's lower and top sides
% (units of Pa)

% beta_star: variable that denotes the inverse velocity of an air molecule 
% at a particular average temperature in the 3D structure
% (units of s/m)

% lamdba: mean free path of an air molecule at a particular viscosity,
% pressure and temperature
% (units of m)

% Kn: Knudsen number. This is a non-dimensional number used to
% characterizing whether we are operating at the free molecular (high Kn)
% or continuum (low Kn) regimes
% (no units)

% delta: rarefaction coefficient, a dimensionless number 
% (no units)

% In this section, we will also model the temperature difference in the
% nanocardboard sheets. To do so, we will model the conductance across the
% thin alumina channel walls, and the air in and out of the cardboard
% structure itself. We will neglect radiation, and assume that the bottom
% face is illuminate by incident light that is absorbed at 90%. This is in
% order to get more accurate approximations of the actual force that is
% being generated by the structure.

% A1: cross-sectional area of the open channels
% (units of m2)

% R1: thermal resistance against heat conduction of the air column in
% between the alumina walls
% (units of Ohms)

% A2: cross-sectional area occupied by alumina in the cardboard structure
% (units of m2)

% R2: thermal resistance of the vertical alumina walls
% (units of m2)

% A3: cross-sectional area occupied by air trapped inside the cardboard
% (units of m)

% R3: thermal resistance of the air trapped within the walls of the
% structure

% psi: proportion of the absorbed optical flux dissipated upward through
% the nanocardboard. This is assumed to be 0.5.
% (unitless)

% epsilon: apsorption. This is assumed to be 1.
% (unitless)

% t: alumina wall thickness.
% (units of m)

% T_star: average temperature between the inner side and outer side of the 
% 3D structure. With nanocardboard as an example, T_star would be the 
% average between the temperature on the plate's lower and top sides
% (units of K) 
% ----------------------------------------------------------------------- %
lambda = (mu./P).*sqrt(3.14.*kb.*T./(2.*m));
Kn = lambda./A; % Knudsen number
delta = sqrt(3.14)./(2.*Kn); 

k_air = k_0./(1+3.116*lambda/L);
%k_air = k_0./(1+0.000076.*T./(P.*L)) - Wu model;
k_alumina = 1.8;
t = chan_param(6); % ALD thickness, which varies from 25 to 100 nm.
A1 = X*B*A; % channels * width * length
A_cell = (X*B*A + S*B*X);
A2 = X*((B+2*t)*(A+2*t)-B*A); 
A3 = A*B*X/phi - X*(B+2*t)*(A+2*t);% Area that remains after 
% subtracting the channels and the ALD thickness surface
% Secondly, we will write expressions to solve for the various types of
% resistances in the plate:
R1 = L./(k_air.*A1);
R2 = L/(k_alumina.*A2);
R3 = (L-2*t)./(k_air.*A3);
% Now, we can begin to compute the deltaT term across the thickness of the
% nanocardboard:
psi = 0.5;
epsilon = 1;
% Recall that L is the thickness of the structure, while A and B are
% channel parameters depending on the specific geometry
num = psi.*epsilon.*I_sun.*A_cell*(1-phi);
den = ((1./R1)+(1./R2)+(1./R3));
deltaT = num./den;
T_star = (deltaT + 2.*T)./2;
P_star = P;
beta_star = sqrt(m./(2.*kb.*T_star)); 
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Mass flow rate across channel parameters
% ----------------------------------------------------------------------- %
% Now, in this section, we can begin to use the previously defined
% variables to construct equations to find vft and vout for the particular
% 3D geometry in question. We begin with the mass flow rate across the
% structure's channels, whose equation has two important leading
% coefficients. These are:

% alpha 
% (units of kg/s*Pa)

% gamma
% units of (kg/s*K)
% ----------------------------------------------------------------------- %
alpha = (delta./6 + 1).*(1 + 0.25./sqrt(delta)).*(A*A).*B.*beta_star./L;
gamma = (1.1./(1.5+delta)).*(A*A).*B.*P_star.*beta_star./(T_star.*L);
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Solving for vout
% ----------------------------------------------------------------------- %
% Next, in this section, we continue with the derivation of equations to
% compute the outlet velocity of the 3D structure. As was detailed in the
% supplementary report, this resulted in a quadratic, which has the
% following three coefficients:
a = 1-(A_out/A_total)^2;
b = 2*A_out*A*B./(A_total.*alpha.*phi);
c = - 2.*gamma.*deltaT./(rho.*alpha); 
% ----------------------------------------------------------------------- %
% Now, with this information, we can solve the quadratic equation, using 
% matlab's built-in roots command
% First, we initialize a solution vector to store the output roots
sol = zeros(length(altitude),2);
% Then, we loop through all the quadratics corresponding to each altitude,
% and solve individually for each. 
for kk = 1:length(altitude)
    sol(kk,1:2) = roots([a,b(kk),c(kk)]);
end

% Finally, we take the positive roots as the vout
vout = sol(:,2);
% Using the solution to vout, we can now calculate vft
vft = A_out.*vout./A_total;
deltaP = 0.5.*rho.*(vout.^2 - vft.^2);
% ----------------------------------------------------------------------- %


% ----------------------------------------------------------------------- %
% Fitting coefficients
% ----------------------------------------------------------------------- %
% In this section, we use the solutions for vft and vout to calculate the
% forces produced by the 3D structure. C1 and C2 are two coefficients 
% obtained from fitting data from the Ansys Fluent simulations
% ----------------------------------------------------------------------- %
if option == 2
    % Cone Geometry
    C1 = 1.16;
    C2 = 0.86;
    fit = C1 .* 8 .* D .* mu .* vft + C2 .* rho .* A_out .* vout.^2;
    elseif option == 3 
    % Sphere Geometry
    C1 = 1.30;
    C2 = 0.89;
    fit = C1 .* 8 .* D .* mu .* vft + C2 .* rho .* A_out .* vout.^2;
    else
    % Rocket Geometry
    C1 = 1.44;
    C2 = 0.93;
    fit = C1 .* 8 .* D .* mu .* vft + C2 .* rho .* A_out .* vout.^2;
end
% This is the general fitting equation


% Calculating the mass of the walls
A_channel = A*B;
n_channel = floor(A_in/A_channel);
%A_walls = n_channel*(2*A + 2*B)*(L-2*t);
den_alumina = 3950; % kg/m3
vol_walls = n_channel*(L-2*t)*((B+2*t)*(A+2*t)-B*A);
    
% Calculating the net lift that can be generated by the structure. This
% will be equal to the propulsion force minus the structure's own weight
net_lift = fit - 9.8*(A_solid*0.001 + vol_walls*den_alumina);
% ----------------------------------------------------------------------- %

% Writing the simplified temperatures
gam = phi*gamma./(rho*A*B);
vft2 = 0.0275*I_sun./P;
vft3 = gam.*deltaT;

% Writing the aereal densities
aerial = fit./(9.8*(A_total-A_in));

end