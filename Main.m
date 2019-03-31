clear all; 
clc; 
dbstop if error;
close('all');

%% Model parameters
H.N         = 300;                      % Number of elements

%% System parameters
H.EA        = 1.4544E9;                 % [N]  Young's modulus
H.A         = H.EA/210E9;               % [m^2] Area of cross-section
H.rho       = 7850;                     % [kg/m^3] Mooring line material density
H.rho_w     = 1025;                     % [kg/m^3] Water density
H.delta     = -10;                      % [m] Depth of anchor
H.g         = -9.81;                    % [m/s^2] Gravitational constant
H.vsl_draft = -10;                      % [m] Vessel draft
H.H         = 150;                      % [m] Water level
H.R         = 2*H.H;                    % [m] Anchor radius
H.Msub      = 315.36;                   % [kg/m] Submerged weight
H.L         = 1.10*sqrt(H.H^2 + H.R^2); % [m] Length of mooring line
H.dL        = H.L/H.N;                  % [m] Length of elements
H.k_soil    = 20000;                     % [(N/m)/m] Soil stiffness per unit length
H.Fpoint    = 0;

u0          = get_initial_conditions(H,0);       % Initial conditions vec u0 (1 for plot, 0 for nothing)
% u0          = [[0:H.dL:H.L].'; zeros(H.N,1); 0;]; % base case
u0          = [u0(1:H.N); u0(H.N); u0(H.N+2:end)]; % condition for variation of offset
theta = pi/6;
% u0          = [[0:H.dL:cos(theta)*H.L].'; linspace(H.vsl_draft,H.L*sin(theta)+H.vsl_draft,H.N+1).']; % final case

fn_ID       = [1 H.N+1 H.N+2 2*(H.N+1)];      % Coordinates of fixed nodes
max_it_step_size = 0.2;                  % [m] Maximum step size of u

%% Calculation of stable displacement vector
[u_final,umatrix,CVG,F_CVG] = ...                           % converged displacement vector u
    newton_raphson(u0, fn_ID, max_it_step_size,H);
x_f = u_final(1:H.N+1);
z_f = u_final(H.N+2:end);

figure()
plot(x_f,z_f)
xlabel('x-location [m]')
ylabel('Displacement [m]')
hold on
% figure()
% plot(u0(H.N+2:(end)))

alpha = (1327600/(10*H.A*H.rho));
x_an = linspace(-H.L/2,H.L/2,H.N);
y_an = alpha*cosh(x_an/alpha)-alpha;

% plot(x_an+H.L/2,y_an-max(y_an))


for ii = 1317:CVG-2
    figure(2)
    clf
    plot(umatrix(1:H.N+1,ii),umatrix(H.N+2:end,ii),'k')
    axis equal
    title(['cvg_ID =' num2str(ii)])
    hold on
    xcor = umatrix(225,CVG);
    zcor = umatrix(225+H.N+2,CVG);
    text(xcor-10,zcor,['\leftarrow']);
    text(xcor-120,zcor,['z-coordinate =' num2str(zcor)]);
    pause(0.1)
    
    line([0 300],[0 0])
    hold on
end

