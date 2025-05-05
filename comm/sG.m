function r = sG(p, u)
% sG — RHS for trait-structured biomass-water-hydrology (bwh) community model
%
% Computes:
%   - Biomass dynamics B(x,chi)
%   - Soil water dynamics W(x)
%   - Surface water dynamics H(x)
%
% The model includes:
%   - Trait diffusion (finite differences)
%   - Space diffusion (finite elements)
%   - Trait-dependent uptake, mortality, and infiltration

% === Unpack field variables ===
N = p.N;            % Number of traits
n = p.np;           % Number of spatial grid points

b = reshape(u(1:N*n), n, N);               % Biomass: B(x,chi)
w = u(N*n+1 : (N+1)*n);                    % Soil water
h = u((N+1)*n+1 : (N+2)*n);                % Surface water

% === Unpack parameters ===
par = u(p.nu+1:end);
[pp, Lam0, Ga, A, R, L0, f, Q, ...
 Kmin, Kmax, Mmin, Mmax, Ymin, Ymax, ...
 Db, Dw, Dh, Dchi, chimin, chimax] = deal(par{:});

% === FEM Matrices ===
K = p.mat.K;             % Stiffness matrix (diffusion operator)
M = p.mat.M(1:n, 1:n);   % Mass matrix

% === Precomputed quantities ===
ov  = ones(n, 1);
bt  = zeros(n, 1);   % Trait-aggregated biomass: sum_i b_i(x)
btt = zeros(n, 1);   % Weighted biomass: sum_i Yi * b_i(x)

chii = linspace(chimin, chimax, N);  % Trait values
dchi = chii(2) - chii(1);            % Trait discretization step

% === Compute trait-integrated quantities ===
for i = 1:N
    bi = b(:, i);
    Yi = Ymax + chii(i) * (Ymin - Ymax);
    bt  = bt + bi;         % Total biomass
    btt = btt + Yi * bi;   % Weighted by Yi
end

% Infiltration and evaporation functions (nonlinear terms)
I = A * (btt + f * Q) ./ (btt + Q);        % Infiltration rate
L = L0 ./ (1 + R * bt);                    % Evaporation rate

% === Initialize residual vector for all components ===
r = zeros(p.nu, 1);

% === Biomass dynamics for each trait ===
for i = 1:N
    bi = b(:, i);

    % Trait-dependent uptake and mortality
    Ki   = Kmax + chii(i) * (Kmin - Kmax);
    Mi   = Mmax + chii(i) * (Mmin - Mmax);
    Lami = Lam0 * Ki ./ (bt + Ki);  % Resource-limited growth

    % Trait diffusion (finite-difference approximation with Neumann BC)
    switch i
        case 1
            bcc = (b(:, 2) - 2 * b(:, 1)) / dchi^2;
        case N
            bcc = (-2 * b(:, N) + b(:, N - 1)) / dchi^2;
        otherwise
            bcc = (b(:, i+1) - 2 * b(:, i) + b(:, i-1)) / dchi^2;
    end

    % RHS for B_i(x)
    growth = Lami .* w .* bi;
    mortality = Mi * bi;
    traitDiff = Dchi * bcc;
    spaceDiff = Db * K * bi;

    r_bi = -M * (growth - mortality + traitDiff) + spaceDiff;

    % Insert into global residual
    r((i - 1) * n + 1 : i * n) = r_bi;
end

% === Water (w) dynamics ===
r_w = -M * (I .* h - L .* w - Ga * bt .* w) + Dw * K * w;

% === Surface water (h) dynamics ===
r_h = -M * (pp - I .* h) + Dh * K * h;

% === Concatenate all components into final residual vector ===
r = [r; r_w; r_h];
