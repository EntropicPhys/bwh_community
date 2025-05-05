function J = sGjac(p, u)
% sGjac — Jacobian of the bwh-community model
% Computes the Jacobian matrix of the system:
%   - Includes finite-difference chi-diffusion
%   - Includes spatial diffusion terms
%   - Uses p2pglob.gu for memory allocation if enabled

% === Unpack dimensions and fields ===
N  = p.N;
n  = p.np;
ng = (N + 2) * n;  % total system size

w = u(N*n+1     : (N+1)*n);
h = u((N+1)*n+1 : (N+2)*n);

% === Unpack parameters ===
par = u(p.nu+1:end);
[pp, Lam0, Ga, A, R, L0, f, Q, ...
 Kmin, Kmax, Mmin, Mmax, Ymin, Ymax, ...
 Db, Dw, Dh, Dchi, chimin, chimax] = deal(par{:});

% === FEM matrices ===
K = p.mat.K;
M = p.mat.M(1:n, 1:n);

% === Precompute trait info ===
chii = linspace(chimin, chimax, N);
dchi = chii(2) - chii(1);
bt   = zeros(n, 1);
btt  = zeros(n, 1);
Yi   = zeros(N, 1);
Ki   = zeros(N, 1);
Mi   = zeros(N, 1);

for i = 1:N
    Bi = u((i-1)*n + 1 : i*n);
    bt = bt + Bi;
    Yi(i) = Ymax + chii(i) * (Ymin - Ymax);
    Ki(i) = Kmax + chii(i) * (Kmin - Kmax);
    Mi(i) = Mmax + chii(i) * (Mmin - Mmax);
    btt = btt + Yi(i) * Bi;
end

% === Infiltration and evaporation ===
ov  = ones(n, 1);
I   = A * (btt + f * Q) ./ (btt + Q);
L   = L0 ./ (1 + R * bt);
dL  = -L0 * R ./ (1 + R * bt).^2;

% === Initialize sparse Jacobian ===
if p2pglob.nzi == 0
    J = zeros(ng);  % first allocation as dense
else
    J = p2pglob.gu; % reuse allocated memory
end

% === Main Jacobian assembly ===
for i = 1:N
    Bi    = u((i-1)*n+1 : i*n);
    Lami  = Lam0 * Ki(i) ./ (bt + Ki(i));
    dLami = -Lam0 * Ki(i) ./ (bt + Ki(i)).^2;
    dI    = (A * Yi(i) ./ (btt + Q)) - ...
            (A * Yi(i) * (btt + f * Q) ./ (btt + Q).^2);

    % Derivatives D_j f_i
    for j = 1:N
        djfi = dLami .* w .* Bi + (Lami .* w - Mi(i)) * (i == j);
        J((i-1)*n+1:i*n, (j-1)*n+1:j*n) = ...
            -M * spdiags(djfi, 0, n, n) + (i == j) * Db * K;
    end

    % Derivative w.r.t. w
    dwfi = Lami .* Bi;
    J((i-1)*n+1:i*n, ng-2*n+1:ng-n) = -M * spdiags(dwfi, 0, n, n);

    % Derivative of fw w.r.t. b
    dbfw = dI .* h - dL .* w - Ga * w;
    J(N*n+1:(N+1)*n, (i-1)*n+1:i*n) = -M * spdiags(dbfw, 0, n, n);

    % Derivative of fh w.r.t. b
    dbfh = -dI .* h;
    J((N+1)*n+1:(N+2)*n, (i-1)*n+1:i*n) = -M * spdiags(dbfh, 0, n, n);

    % === Trait diffusion contribution (finite difference in chi) ===
    switch i
        case 1
            % Left Neumann BC
            J(1:n, 1:n) = J(1:n, 1:n) - (Dchi/dchi^2) * M;
            J(1:n, n+1:2*n) = J(1:n, n+1:2*n) + (Dchi/dchi^2) * M;

        case N
            % Right Neumann BC
            J((N-1)*n+1:N*n, (N-2)*n+1:(N-1)*n) = ...
                J((N-1)*n+1:N*n, (N-2)*n+1:(N-1)*n) - (Dchi/dchi^2) * M;
            J((N-1)*n+1:N*n, (N-1)*n+1:N*n) = ...
                J((N-1)*n+1:N*n, (N-1)*n+1:N*n) + (Dchi/dchi^2) * M;

        otherwise
            % Bulk central diff
            J((i-1)*n+1:i*n, (i-2)*n+1:(i-1)*n) = ...
                J((i-1)*n+1:i*n, (i-2)*n+1:(i-1)*n) - (Dchi/dchi^2) * M;
            J((i-1)*n+1:i*n, (i-1)*n+1:i*n) = ...
                J((i-1)*n+1:i*n, (i-1)*n+1:i*n) + 2 * (Dchi/dchi^2) * M;
            J((i-1)*n+1:i*n, i*n+1:(i+1)*n) = ...
                J((i-1)*n+1:i*n, i*n+1:(i+1)*n) - (Dchi/dchi^2) * M;
    end
end

% === Derivatives for water (w) and surface water (h) ===
dwfw = -(L + Ga * bt);
J(N*n+1:(N+1)*n, N*n+1:(N+1)*n) = Dw * K - M * spdiags(dwfw, 0, n, n);

dhfw = I;
J(N*n+1:(N+1)*n, (N+1)*n+1:(N+2)*n) = -M * spdiags(dhfw, 0, n, n);

dhfh = -I;
J((N+1)*n+1:(N+2)*n, (N+1)*n+1:(N+2)*n) = Dh * K - M * spdiags(dhfh, 0, n, n);

% === Finalize as sparse matrix (if first time) ===
if p2pglob.nzi == 0
    J = sparse(J);
    p2pglob.gu = J;
    p2pglob.nzi = 1;
end
