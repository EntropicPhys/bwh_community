function r = sG(p, u)
% sG — RHS of the bwh-community model with explicit chi-diffusion
% This version treats chi-diffusion explicitly (e.g., for IMEX time integration),
% and omits spatial diffusion (Db = Dw = Dh = 0).

% --- Unpack model dimensions and state variables ---
N = p.N;                 % Number of traits
n = p.np;                % Number of spatial grid points

% Reshape biomass field B(x, χ) into a (space × trait) matrix
b = reshape(u(1:N*n), n, N);

% Extract scalar fields (soil water and surface water)
w = u(N*n+1 : (N+1)*n);
h = u((N+1)*n+1 : (N+2)*n);

% --- Extract parameters from the solution vector ---
par = u(p.nu+1:end);
[pp, Lam0, Ga, A, R, L0, f, Q, ...
 Kmin, Kmax, Mmin, Mmax, Ymin, Ymax, ...
 Db, Dw, Dh, Dchi, chimin, chimax] = deal( ...
    par(1), par(2), par(3), par(4), par(5), ...
    par(6), par(7), par(8), par(9), par(10), ...
    par(11), par(12), par(13), par(14), par(15), ...
    par(16), par(17), par(18), par(19), par(20));

% --- Retrieve FEM matrices ---
K = p.mat.K;               % Spatial stiffness (Laplacian) matrix
M = p.mat.M(1:n,1:n);      % Mass matrix for scalar fields

% --- Initialize variables ---
ov   = ones(n,1);          % All-ones vector
bt   = zeros(n,1);         % Total biomass (b̃)
btt  = zeros(n,1);         % Weighted biomass (b̃̃)
r    = zeros(N*n,1);       % Output for biomass equations

% Trait mesh size
dchi = (chimax - chimin) / (N - 1);

% === Compute total and weighted biomass ===
for i = 1:N
    chii = chimin + (i-1)*dchi;                          % Trait value χ_i
    bt   = bt + b(:,i);                                  % Total biomass
    Yi   = Ymax + chii * (Ymin - Ymax);                  % Trait-dependent water-use efficiency
    btt  = btt + b(:,i) * Yi;                            % Weighted biomass
end

% === Infiltration and evaporation functions ===
I = A * (btt + f*Q) ./ (btt + Q);                        % Infiltration rate I(x)
L = L0 ./ (1 + R * bt);                                  % Evaporation rate L(x)

% === Biomass equations for each trait ===
for i = 1:N
    chii = chimin + (i-1)*dchi;

    % Trait-dependent coefficients
    Ki   = Kmax + chii * (Kmin - Kmax);
    Mi   = Mmax + chii * (Mmin - Mmax);
    Lami = Lam0 * Ki ./ (bt + Ki);                       % Growth efficiency

    bi = b(:,i);                                          % Biomass for trait i

    % Trait diffusion term (finite difference in χ with Neumann BCs)
    switch i
        case 1
            bcc = (b(:,2) - 2*b(:,1)) / dchi^2;           % Left boundary
        case N
            bcc = (-2*b(:,N) + b(:,N-1)) / dchi^2;        % Right boundary
        otherwise
            bcc = (b(:,i+1) - 2*b(:,i) + b(:,i-1)) / dchi^2;
    end

    % Biomass dynamics equation:
    % r1 = -M * [ growth - mortality + trait diffusion ] + spatial diffusion
    r1 = -M * (Lami .* w .* bi - Mi * bi + Dchi * bcc) + Db * K * bi;

    % Store result into appropriate block in r
    r((i-1)*n+1 : i*n) = r1;
end

% === Water equation ===
% r2 = -M * (I*h - L*w - Ga*bt*w) + Dw*K*w
r2 = -M * (I .* h - L .* w - Ga * bt .* w) + Dw * K * w;

% === Surface water equation ===
% r3 = -M * (pp - I*h) + Dh*K*h
r3 = -M * (pp - I .* h) + Dh * K * h;

% === Assemble full residual vector ===
r = [r; r2; r3];

