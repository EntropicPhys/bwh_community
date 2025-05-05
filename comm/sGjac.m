function J = sGjac(p, u)
global p2pglob
% sGjac — Jacobian of the bwh-community model
% Computes the Jacobian matrix of the system:
%   - Includes finite-difference chi-diffusion
%   - Includes spatial diffusion terms
%   - Uses p2pglob.gu for memory allocation if enabled

% --- Jacobian assembly for the bwh-community model ---
% Input: 
%   p — PDE2PATH problem structure
%   u — state vector [b; w; h; par]
% Output: 
%   J — Jacobian matrix of size (ng × ng), where ng = (N+2)*n

N = p.N;              % Number of traits
n = p.np;             % Number of spatial grid points
ng = (N + 2) * n;     % Total number of unknowns (b, w, h)

% Extract scalar water fields
w = u(N*n+1 : (N+1)*n);
h = u((N+1)*n+1 : (N+2)*n);

% Extract parameters
par = u(p.nu+1:end);
[pp, Lam0, Ga, A, R, L0, f, Q, ...
 Kmin, Kmax, Mmin, Mmax, Ymin, Ymax, ...
 Db, Dw, Dh, Dchi, chimin, chimax] = deal( ...
    par(1), par(2), par(3), par(4), par(5), ...
    par(6), par(7), par(8), par(9), par(10), ...
    par(11), par(12), par(13), par(14), par(15), ...
    par(16), par(17), par(18), par(19), par(20));

% FEM matrices
K = p.mat.K;
M = p.mat.M(1:n, 1:n);

% Trait and biomass terms
bt   = zeros(n, 1);       % Total biomass b̃(x)
btt  = zeros(n, 1);       % Weighted biomass b̃̃(x)
Yi   = zeros(N, 1);       % Trait-dependent Yi
Ki   = zeros(N, 1);       % Trait-dependent Ki
Mi   = zeros(N, 1);       % Trait-dependent Mi

% Trait discretization
dchi = (chimax - chimin) / (N - 1);

% --- Compute b̃, b̃̃, Yi, Ki, Mi ---
for i = 1:N
    chii = chimin + (i - 1) * dchi;
    Bi = u((i - 1)*n + 1 : i*n);
    bt = bt + Bi;
    Yi(i) = Ymax + chii * (Ymin - Ymax);
    Ki(i) = Kmax + chii * (Kmin - Kmax);
    Mi(i) = Mmax + chii * (Mmin - Mmax);
    btt = btt + Yi(i) * Bi;
end

% --- Compute infiltration, evaporation, and their derivatives ---
ov  = ones(n, 1);
I   = A * (btt + f*Q) ./ (btt + Q);               % Infiltration function
L   = L0 ./ (1 + R * bt);                         % Evaporation function
dL  = -L0 * R ./ (1 + R * bt).^2;                 % Derivative of evaporation

% --- Allocate Jacobian ---
if p2pglob.nzi == 0
    J = zeros(ng);
else
    J = p2pglob.gu;
end

% === Block Jacobian assembly ===
for i = 1:N
    Bi    = u((i - 1)*n + 1 : i*n);
    Lami  = Lam0 * Ki(i) ./ (bt + Ki(i));
    dLami = -Lam0 * Ki(i) ./ (bt + Ki(i)).^2;
    
    % Derivative of infiltration with respect to biomass (trait-wise)
    dI = A * Yi(i) ./ (btt + Q) - A * Yi(i) * (btt + f*Q) ./ (btt + Q).^2;

    % --- Loop over all trait interactions ---
    for j = 1:N
        is_diag = (i == j);
        Bj = u((j - 1)*n + 1 : j*n);

        % Diagonal and off-diagonal trait interactions for biomass equation
        djfi = dLami .* w .* Bi + (Lami .* w - Mi(i)) * is_diag;
        J((i - 1)*n + 1 : i*n, (j - 1)*n + 1 : j*n) = ...
            -M * spdiags(djfi, 0, n, n) + is_diag * Db * K;
    end

    % --- Coupling to water equation (∂f_i/∂w) ---
    dwfi = Lami .* Bi;
    J((i - 1)*n + 1 : i*n, ng - 2*n + 1 : ng - n) = -M * spdiags(dwfi, 0, n, n);

    % --- Coupling from biomass to water equations (∂fw/∂b, ∂fh/∂b) ---
    dbfw = dI .* h - dL .* w - Ga * w;
    dbfh = -dI .* h;
    J(N*n + 1 : (N+1)*n, (i - 1)*n + 1 : i*n)      = -M * spdiags(dbfw, 0, n, n);
    J((N+1)*n + 1 : (N+2)*n, (i - 1)*n + 1 : i*n)  = -M * spdiags(dbfh, 0, n, n);

    % --- Trait-diffusion block (finite difference in χ) ---
    row = (i - 1)*n + 1 : i*n;
    if i == 1
        % Left boundary
        J(row, row)           = J(row, row)           - (Dchi / dchi^2) * M;
        J(row, row + n)       = J(row, row + n)       + (Dchi / dchi^2) * M;
    elseif i == N
        % Right boundary
        J(row, row - n)       = J(row, row - n)       - (Dchi / dchi^2) * M;
        J(row, row)           = J(row, row)           + (Dchi / dchi^2) * M;
    else
        % Interior
        J(row, row - n)       = J(row, row - n)       - (Dchi / dchi^2) * M;
        J(row, row)           = J(row, row)           + 2 * (Dchi / dchi^2) * M;
        J(row, row + n)       = J(row, row + n)       - (Dchi / dchi^2) * M;
    end
end

% === Water equation: dwfw = ∂fw/∂w, dhfw = ∂fw/∂h ===
dwfw = -(L + Ga * bt);
J(N*n + 1 : (N+1)*n, N*n + 1 : (N+1)*n)         = Dw * K - M * spdiags(dwfw, 0, n, n);

dhfw = I;
J(N*n + 1 : (N+1)*n, (N+1)*n + 1 : (N+2)*n)     = -M * spdiags(dhfw, 0, n, n);

% === Surface water equation: ∂fh/∂h ===
dhfh = -I;
J((N+1)*n + 1 : (N+2)*n, (N+1)*n + 1 : (N+2)*n) = Dh * K - M * spdiags(dhfh, 0, n, n);

% === Store sparse pattern if first run ===
if p2pglob.nzi == 0
    J = sparse(J);
    p2pglob.gu = J;
    p2pglob.nzi = 1;
end
end
