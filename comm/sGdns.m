function r = sGdns(p, u)
% sGdns — RHS of the bwh-community model with explicit chi-diffusion
% This version treats chi-diffusion explicitly (e.g., for IMEX time integration),
% and omits spatial diffusion (Db = Dw = Dh = 0).

% === Unpack dimensions and state variables ===
N = p.N;            % Number of traits
n = p.np;           % Number of spatial points

b = reshape(u(1:N*n), n, N);                     % Biomass B(x,chi)
w = u(N*n+1 : (N+1)*n);                          % Soil water W(x)
h = u((N+1)*n+1 : (N+2)*n);                      % Surface water H(x)

% === Unpack parameter vector ===
par = u(p.nu+1:end);
[pp, Lam0, Ga, A, R, L0, f, Q, ...
 Kmin, Kmax, Mmin, Mmax, Ymin, Ymax, ...
 Db, Dw, Dh, Dchi, chimin, chimax] = deal(par{:});

% === FEM Matrices (used only for mass matrix) ===
K = p.mat.K;                % Stiffness matrix (unused here)
M = p.mat.M(1:n, 1:n);      % Mass matrix
dchi = (chimax - chimin) / (N - 1);   % Trait mesh spacing

% === Initialize RHS components ===
r  = zeros(N*n, 1);         % Biomass part
bt = zeros(n, 1);           % Total biomass
btt = zeros(n, 1);          % Biomass weighted by trait Yi

% === Compute total biomass and weighted biomass ===
for i = 1:N
    chii = chimin + (i-1)*dchi;
    Yi = Ymax + chii * (Ymin - Ymax);
    bi = b(:, i);
    bt = bt + bi;
    btt = btt + Yi * bi;
end

% === Biomass dynamics: loop over traits ===
for i = 1:N
    chii = chimin + (i-1)*dchi;

    % Trait-dependent coefficients
    Mi   = Mmax + chii * (Mmin - Mmax);
    Ki   = Kmax + chii * (Kmin - Kmax);
    Lami = Lam0 * Ki ./ (bt + Ki);
    I    = A * (btt + f * Q) ./ (btt + Q);
    L    = L0 ./ (1 + R * bt);

    bi = b(:, i);

    % Trait diffusion (explicit in chi)
    switch i
        case 1
            bcc = (b(:, 2) - 2 * b(:, 1)) / dchi^2;
        case N
            bcc = (-2 * b(:, N) + b(:, N-1)) / dchi^2;
        otherwise
            bcc = (b(:, i+1) - 2 * b(:, i) + b(:, i-1)) / dchi^2;
    end

    % RHS for biomass B(x,chi_i) without spatial diffusion
    r_bi = -M * (Lami .* w .* bi - Mi * bi + Dchi * bcc);

    % Store into global RHS vector
    r((i-1)*n+1 : i*n) = r_bi;
end

% === Water and Surface Water Equations ===
% Note: Spatial diffusion terms are omitted
I = A * (btt + f * Q) ./ (btt + Q);
L = L0 ./ (1 + R * bt);

r_w = -M * (I .* h - L .* w - Ga * w .* bt);
r_h = -M * (pp - I .* h);

% === Final assembled RHS ===
r = -[r; r_w; r_h];
