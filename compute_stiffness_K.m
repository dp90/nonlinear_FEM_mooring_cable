function K = compute_stiffness_K(x, fixed_ID, H)
% Determination of stiffness matrix K
dx = H.dL^(1/3); % finite difference delta
nx = length(x); % degrees of freedom
nf = nx; % number of functions
K = zeros(nf,nx); % matrix of zeros
DOF_vec = 1:nx;
DOF_vec(fixed_ID) = [];
for n = DOF_vec
    % create a vector of deltas, change delta_n by dx
    delta = zeros(nx, 1);
    delta(n) = delta(n)+dx;
    dF = A5_eq(x+delta,H)-A5_eq(x-delta,H); % delta F
    K(:, n) = dF/dx/2; % derivatives dF/d_n
end

% Remove fixed DOF
K(:,fixed_ID) = [];
K(fixed_ID,:) = [];

% Make sparse
K = sparse(K);