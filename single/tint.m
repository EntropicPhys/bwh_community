function [p, t1] = tint(p, t1, dt, nt, pmod)
% tint.m — Time integration routine using semi-implicit Euler stepping
%
% Inputs:
%   p     - pde2path problem structure (includes solution, matrices, etc.)
%   t1    - initial time
%   dt    - time step size
%   nt    - number of time steps
%   pmod  - plotting interval (plot every pmod steps)
%
% Outputs:
%   p     - updated problem structure (including final solution)
%   t1    - final time after nt steps

%--- Initialization ---%
n = 0;                % time step counter
t = t1;               % current simulation time

% Extract model parameters from the end of the solution vector
par = p.u(p.nu+1:end);
Db  = par(15);        % diffusion coefficient for biomass
Dw  = par(16);        % diffusion coefficient for soil water
Dh  = par(17);        % diffusion coefficient for surface water

%--- Assemble system matrix for semi-implicit stepping ---%
% Construct block-diagonal diffusion stiffness matrix using Kronecker product
K = kron(diag([Db, Dw, Dh]), p.mat.K);

% Construct semi-implicit system matrix for (I + dt*K)
Lam = p.mat.M + dt * K;

% Precompute LU decomposition for efficiency
[L, U, P, Q, R] = lu(Lam);

%--- Time-stepping loop ---%
while (n < nt)
    % Compute the nonlinear right-hand side (explicit part)
    f = -sGdns(p, p.u);  % DNS-specific residual (e.g., nonlinear reaction terms)
    % Compute the RHS of the linear system: M * u_old + dt * f
    g = p.mat.M * p.u(1:p.nu) + dt * f;
    % Solve the linear system using precomputed LU factors
    p.u(1:p.nu) = Q * (U \ (L \ (P * (R \ g))));
    % Update time and step counter
    n = n + 1;
    t = t + dt;
    % Plot current solution every 'pmod' steps
    if mod(n, pmod) == 0
        r = norm(resi(p, p.u), 'inf');
        tits = ['t = ', mat2str(t, 3), ', r = ', mat2str(r, 3)];
        plotsol(p, p.plot.ifig, p.plot.pcmp, p.plot.pstyle);
        title(['u_1, ', tits], 'fontsize', 12);
        set(gca, 'fontsize', 12);
        drawnow;
    end
end
%--- Return final time ---%
t1 = t;
end
