function F = compute_external_forces(u,H)

N = length(u);                                  % Number of dof's
x = u(1 : N/2);                                 % x-positions
z = u( N/2+1 : N );                             % z-positions
Fx = zeros(N/2,1);
Fz = zeros(N/2,1);

k = H.EA/H.dL;

dx = diff(x);
dz = diff(z);
lc = sqrt(dx.^2 + dz.^2);
fx1 = k*(lc - H.dL*ones(H.N,1) ).*(dx./lc);
fx2 = [fx1;0] - [0;fx1];
fz1 = k*(lc - H.dL*ones(H.N,1) ).*(dz./lc);
fz2 = [fz1;0] - [0;fz1];
F_moor = [fx2 ; fz2];

F_soilx = zeros(N/2,1);
a = z<0;
F_soilz = -a.*H.k_soil.*z;
F_soil = [F_soilx ; F_soilz];

% F_soil increases linear for z<0, and is 0 at z>=0
F_g = [zeros(N/2,1) ; ones(N/2,1)*(H.rho-H.rho_w)*H.A*H.dL*H.g];

F_px = zeros(N/2,1);

F_pz = zeros(N/2,1);
    
if mod(H.N,2) ~= 0
    F_pz(N/4) = H.Fpoint/2;
    F_pz(N/4+1) = H.Fpoint/2;
else
    F_pz((N+2)/4) = H.Fpoint;
end

F_point = [F_px ; F_pz];



F = F_moor + F_point + F_g + F_soil;

finished = H.finished;

% cvg_ID = H.cvg_ID;
%     if (cvg_ID == 10) || (cvg_ID == 50) || (cvg_ID == 100) || (cvg_ID == 200) || (finished==1) ==1
%         figure(3)
%         plot(F_moor)
%         hold on
%     end   
end