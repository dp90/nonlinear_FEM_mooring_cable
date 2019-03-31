function [u_stable,umatrix,CVG,F_CVG] = newton_raphson(IC, fixed_node_ID, max_it_step_size,H )
% This function finds the final converged system configuration (cable
% shape)
%{

Input
#1: eq:             link to the function of your system
#2: IC:             the initial guess configuration
#3: fixed_node_ID:  list of fixed DOFs
#4: max_it_step_size:    the maximum step size of each iteration the solver is allowed to take in [m]. 

Outputs:
#1: u_stable:       The stable configuration of the system 

%}

%% Prep
if isrow(IC)
    IC = IC.';
end

%% Iterative solution using Newton-Raphson method
finished = 0;
cvg_ID = 0;
u = IC;

umatrix1 = zeros(2*H.N+2,1000);

while finished == 0
    
    clc
    cvg_ID = cvg_ID + 1
    
    H.cvg_ID = cvg_ID;
    H.finished = finished;
    
    % Compute
    K = compute_stiffness_K(u, fixed_node_ID, H);
    F = compute_external_forces(u,H);
    
    % Check if singular
    if condest(K) > 1E15
        error('System probably not constrained properly; rigid body motions are possible. Check ''fixed_node_ID''.')
    end
    
    % Remove fixed DOF
    F_CVG = F;
    F(fixed_node_ID) = [];
    
    % Get displacements
    u_incr = K \ F;
    
    % Make sure iteration increment doesn't exceed maximum step size
    step_size = max(abs(u_incr));
    if step_size > max_it_step_size
        % Reduce step size
        u_incr = u_incr * max_it_step_size / step_size;
    end
    
    % Add fixed DOF
    for k = 1:length(fixed_node_ID)
        u_incr = [u_incr(1:fixed_node_ID(k)-1);0;u_incr(fixed_node_ID(k):end)];
    end
    
    % Add the displacement increment
    u_old = u;
    u = u - u_incr;
    
    % Check for convergence
    if max((u_old - u).^2) < 1E-5 % originally at 1E-14
        finished = 1;
        CVG = cvg_ID;
    end
    
%     if (cvg_ID == 10) || (cvg_ID == 50) || (cvg_ID == 100) || (cvg_ID == 200) || (finished==1) ==1
%         figure(2)
%         plot(F)
%         hold on
%     end   

    % Check if stuck
    if cvg_ID > 1999
        plot(u(1:H.N+1), u(H.N+2:H.N*2+2))
        error('1000th iteration. Solver is probably stuck; check!')
    finished = 1;
    end
    
    % Plot iterations
%     if static_iteration_plot == 1
%         plot(u_old(1:q2),u_old(q2+1:q),'r');
%         line(u(1:q2),u(q2+1:q),'color',[.8,.8,.8]*exp(-cvg_ID/20),'linewidth',1);
%     end

%% switch on/off depending on question
umatrix1(:,cvg_ID) = u;

end

% Return the stable config
u_stable = u;

%% switch on/off depending on question
umatrix = umatrix1;

end

