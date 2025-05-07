function Gu = sGjac(p, u)
% sGjac.m — Jacobian of the residual function for the single-species BWH model
% Computes the Jacobian matrix Gu = G_u(p, u) used for Newton steps
% 
% Inputs:
%   p - pde2path problem structure
%   u - current solution vector [b; w; h; par]
%
% Output:
%   Gu - Jacobian matrix
%--- Unpack state variables ---%
n  = p.np;                        % Number of spatial points
b  = u(1:n);                      % Biomass
w  = u(n+1:2*n);                  % Soil water
h  = u(2*n+1:3*n);                % Surface water
par = u(p.nu+1:end);             % Parameter vector
%--- Unpack parameters ---%
pp   = par(1);  Lam0 = par(2);  Ga   = par(3); 
A    = par(4);  R    = par(5);  L0   = par(6);
f    = par(7);  Q    = par(8);  
Kmin = par(9);  Kmax = par(10);
Mmin = par(11); Mmax = par(12);
Ymin = par(13); Ymax = par(14);
Db   = par(15); Dw   = par(16); Dh = par(17);
chi  = par(18); 
l    = par(19);                  % Spatial scaling factor
%--- FEM operators and utilities ---%
M  = p.mat.M;                    % Mass matrix
K  = p.mat.K;                    % Stiffness matrix
ov = ones(n, 1);                 % Vector of ones
%--- Trait-dependent parameters ---%
Yi = Ymax + chi * (Ymin - Ymax);
Ki = Kmax + chi * (Kmin - Kmax);
Mi = Mmax + chi * (Mmin - Mmax);
%--- Functional responses and derivatives ---%
Lam  = Lam0 * Ki ./ (b + Ki);                         % Growth rate
dLam = -Lam0 * Ki ./ (b + Ki).^2;                     % dLam/db
I    = A * (Yi * b + f * Q) ./ (Yi * b + Q);          % Infiltration
dI   = A * Yi ./ (Yi * b + Q) ...
     - A * Yi * (Yi * b + f * Q) ./ (Yi * b + Q).^2;  % dI/db
L    = L0 ./ (1 + R * b);                             % Evaporation
dL   = -L0 * R ./ (1 + R * b).^2;                     % dL/db
%--- Local Jacobian block entries ---%
f1b = dLam .* w .* b + Lam .* w - Mi;   f1w = Lam .* b;   f1h = zeros(n, 1);
f2b = dI .* h - Ga * w - dL .* w;       f2w = -L - Ga * b; f2h = I;
f3b = -dI .* h;                         f3w = zeros(n, 1); f3h = -I;
%--- Assemble local Jacobian matrix (blockwise, sparse) ---%
Fu = [
    spdiags(f1b, 0, n, n), spdiags(f1w, 0, n, n), spdiags(f1h, 0, n, n);
    spdiags(f2b, 0, n, n), spdiags(f2w, 0, n, n), spdiags(f2h, 0, n, n);
    spdiags(f3b, 0, n, n), spdiags(f3w, 0, n, n), spdiags(f3h, 0, n, n)
];
%--- Combine diffusion and reaction terms ---%
Diffusion = kron( ...
    [l^2 * Db, 0, 0;
     0, l^2 * Dw, 0;
     0, 0, l^2 * Dh], K);
Gu = Diffusion - M * Fu;
end
