function r = sG(p, u)
% sG.m — Residual function for the single-species BWH model
% Computes the right-hand side of the PDE system for a given state u.
%
% Inputs:
%   p - pde2path problem structure
%   u - state vector [b; w; h; par]
%
% Output:
%   r - residual vector [r1; r2; r3]

%--- Unpack fields and discretization ---%
n = p.np;                          % Number of spatial points
b = u(1:n);                        % Biomass
w = u(n+1:2*n);                    % Soil water
h = u(2*n+1:3*n);                  % Surface water

%--- Unpack parameters ---%
par     = u(p.nu+1:end);          % Parameter vector
pp      = par(1);                 % Precipitation
Lam0    = par(2); Ga     = par(3); % Growth and mortality rates
A       = par(4); R      = par(5); % Infiltration and evaporation parameters
L0      = par(6); f      = par(7);
Q       = par(8); Kmin   = par(9); Kmax = par(10);
Mmin    = par(11); Mmax  = par(12);
Ymin    = par(13); Ymax  = par(14);
Db      = par(15); Dw    = par(16); Dh   = par(17); % Diffusion coefficients
chi     = par(18); l     = par(19);      % Trait value and domain scaling

%--- FEM matrices ---%
K  = p.mat.K;                     % Stiffness matrix
M  = p.mat.M(1:n,1:n);            % Mass matrix
ov = ones(n,1);                   % Vector of ones (for convenience)

%--- Trait-dependent parameters ---%
Yi  = Ymax + chi * (Ymin - Ymax);     % Infiltration capacity
Mi  = Mmax + chi * (Mmin - Mmax);     % Mortality
Ki  = Kmax + chi * (Kmin - Kmax);     % Growth saturation constant

%--- Model terms ---%
Lam = Lam0 * Ki ./ (b + Ki);              % Growth rate
I   = A * (Yi * b + f * Q) ./ (Yi * b + Q); % Infiltration rate
L   = L0 ./ (1 + R * b);                  % Evaporation rate

%--- PDE residuals ---%
r1 = -M * (Lam .* w .* b - Mi * b)      + l^2 * Db * K * b;
r2 = -M * (I .* h - L .* w - Ga * w .* b) + l^2 * Dw * K * w;
r3 = -M * (pp - I .* h)                  + l^2 * Dh * K * h;

%--- Concatenate residuals ---%
r = [r1; r2; r3];
end
