function [p, t1, ts, nc] = tintxs(p, t0, ts, dt, nt, nc, pmod, smod, nffu, varargin)
% tintxs — Time integration with semi-implicit treatment of x-diffusion
%
% Inputs:
%   p     - pde2path problem structure
%   t0    - initial time
%   ts    - time series array (time + residual)
%   dt    - time step
%   nt    - number of steps to integrate
%   nc    - current frame counter
%   pmod  - how often to plot and compute residual
%   smod  - unused (was for saving to disk)
%   nffu  - function handle for the RHS f(u)
%
% Outputs:
%   p     - updated problem structure with solution at final time
%   t1    - final time
%   ts    - updated time series (residual tracking)
%   nc    - updated step counter

% === Extract parameters and dimensions ===
par = p.u(p.nu+1:end);          % Extract parameter vector from solution
Db  = par(15); Dw = par(16); Dh = par(17);  % Diffusion coefficients
N   = p.N;                      % Number of biomass species
n   = p.np;                     % Number of spatial points
ng  = N * n;                    % Total biomass DOFs
u0  = p.u;                      % Store initial state

% === Record initial residual ===
r = norm(resi(p, p.u), 'inf');
ts = [ts, [t0; r]];

% === Build Jacobian for implicit x-diffusion ===
K = p.mat.K;
J = sparse(ng + 2*n, ng + 2*n);  % Full system Jacobian

% Fill diagonal blocks with diffusion matrices
for i = 1:N
    J((i-1)*n+1:i*n, (i-1)*n+1:i*n) = Db * K;  % Biomass components
end
J(N*n+1:(N+1)*n,   N*n+1:(N+1)*n)   = Dw * K;  % Water
J((N+1)*n+1:end,   (N+1)*n+1:end)   = Dh * K;  % Surface water

% === Implicit matrix: (M + dt * J) ===
Lam = p.mat.M + dt * J;

% === LU decomposition (only once, reuse during time loop) ===
[L, U, P, Q, R] = lu(Lam);

% === Time stepping ===
t = t0;
nsteps = 0;

while nsteps < nt
    f = nffu(p, p.u);                  % Evaluate RHS f(u)
    g = p.mat.M * p.u(1:p.nu) + dt * f;
    
    % Solve (M + dt*J) * u^{n+1} = g
    p.u(1:p.nu) = Q * (U \ (L \ (P * (R \ g))));
    
    % Advance time
    t = t + dt;
    nsteps = nsteps + 1;
    
    % Plot and record residual every pmod steps
    if mod(nsteps, pmod) == 0
        r = norm(resi(p, p.u), 'inf');
        ts = [ts, [t; r]];
        tits = ['t = ', mat2str(t, 4), ', r = ', mat2str(r, 3)];
        plotsol(p, p.plot.ifig, p.plot.pcmp, p.plot.pstyle);
        title(['u_1, ', tits], 'fontsize', 12);
        set(gca, 'fontsize', 12);
        drawnow;
    end
    
    % Optional save block (disabled):
    % if mod(nsteps, smod) == 0
    %     ps = p; p = [];
    %     p.t = t; p.ts = ts; p.u = ps.u;
    %     fname = [ps.file.pname, sprintf('%i', nc + nsteps), '.mat'];
    %     save(fname, 'p');
    %     p = ps;
    % end
end

% === Output final time and updated counter ===
t1 = t;
nc = nc + nt;
