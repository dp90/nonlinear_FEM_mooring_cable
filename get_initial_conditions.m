function IC = get_initial_conditions(PRM,dbg)
% This function outputs an intial trial configuration of the mooring cable,
% to prevent numerical instabilities at the first iterations. 
%{

Input
#1: PRM: struct containing all the parameters
#2: dbg: {0,1} 0: nothing
               1: plot the IC

Outputs:
#1: IC:  the initial guess configuration of the system

%}

%% Preparation
L = PRM.L; % [m] length of mooring line
delta = PRM.delta; % [m] depth, should be negative
H = PRM.H; % [m] Water depth, should be positive
vsl_draft = PRM.vsl_draft;
R = PRM.R; % [m] Anchor radius, should be positive
N = PRM.N; % [#] number of elements
g = PRM.g; % [m/s2] gravitational constant; should be negative

% Set dependent parameters
% Get distance between anchor point and connection point
anch_ship_length = sqrt(R^2 + (H-delta+vsl_draft)^2);

%% Checks
if delta > 0
    error('Anchor should be IN the ground (z<0)!')
end
if vsl_draft > 0
    error('Draft should be smaller than zero!')
end
if sign(g) ~= -1
    error('Gravity should point downwards')
end
if L < anch_ship_length
    error('Mooring line length L should be larger than %.0f [m]',anch_ship_length)
end
if L > pi*anch_ship_length/2
    error('Mooring line length L should be smaller than %.0f [m]',pi*anch_ship_length/2)
end
if H <= 0
    error('Water depth should be positive')
end
if R <= 0
    error('Anchor radius should be positive')
end

%% Get IC radius

% Set IC radius equation
% IC_radius_eq = @(delta) asin(anch_ship_length/2/delta) - L/2/delta;
IC_radius_eq = @(delta) sin(L./2./delta) - anch_ship_length./2./delta;

% Get data
IC_radius_guess = linspace(anch_ship_length/2, 10*anch_ship_length, 1E6);
dat = IC_radius_eq(IC_radius_guess);

% Plot
% plot(IC_radius_guess,IC_radius_eq(IC_radius_guess))

% Find when it goes through zero
rel = sign(dat(2:end)) .* sign(dat(1:end-1));
zero_cross_ID = find(rel == -1,1,'first');

% Set radius
IC_radius = IC_radius_guess(zero_cross_ID);

% Get angle of mooring
alpha = L/2/IC_radius;

%% Create IC points

% Set locs
loc_anchor = [0;delta];
loc_ship = [R;H+vsl_draft];

% Get vector from anchor to ship
anchor_ship_vec = loc_ship - loc_anchor;
anchor_ship_n = anchor_ship_vec ./ anch_ship_length;

% Rotate 90 degrees
circle_n_opt1 = [0,-1;1,0] * anchor_ship_n;
circle_n_opt2 = [0,1;-1,0] * anchor_ship_n;

% Select the one which point in the opposite direction as the gravity
% vector
if sign(circle_n_opt1(2)) == -sign(g)
    circle_n = circle_n_opt1;
else
    circle_n = circle_n_opt2;
end

% Get R*
R_star = sqrt(IC_radius^2 - (anch_ship_length/2)^2);

% Get GCS location of origin of circle
circle_org_GCS = loc_anchor + (loc_ship - loc_anchor) / 2 + R_star * circle_n;

% Get vector from circle origin to anchor point
cc_to_anch = loc_anchor - circle_org_GCS;

% Get angle from circle origin to anchor point
angle_cc_ach = atan2(cc_to_anch(2),cc_to_anch(1));

% Get vector of angles for the IC
angle_vec = linspace(angle_cc_ach, angle_cc_ach + 2*alpha, N+1).';

% Create vector of IC
IC_rot_part_x = IC_radius .* cos(angle_vec);
IC_rot_part_z = IC_radius .* sin(angle_vec);
IC = zeros(2*(N+1),1);
IC(1:N+1) = circle_org_GCS(1) + IC_rot_part_x;
IC(N+2:end) = circle_org_GCS(2) + IC_rot_part_z;

%% Plot to validate
if dbg == 1
    hold on
    plot([0,R],[0,0],':k')
    plot([0,R],[H,H],':b')
    plot([0,0],[0,H],':k')
    plot([0,circle_org_GCS(1)],[delta,circle_org_GCS(2)],'--k')
    plot([R,circle_org_GCS(1)],[H+vsl_draft,circle_org_GCS(2)],'--k')
    plot(0,delta,'or')
    plot(R,H+vsl_draft,'or')
    plot(circle_org_GCS(1),circle_org_GCS(2),'ok')
    plot(IC(1:N+1), IC(N+2:end),'--')
    xl = xlim; yl=ylim;
    axis(1.2*[xl(1),xl(2),yl(1),yl(2)])
    axis equal
end

end


