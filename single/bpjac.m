function duGuph=bpjac(p,u) % second derivative for BP continuation 
% === Unpack state vector and parameters ===
n = p.np;
b = u(1:n);
w = u(n+1:2*n);
h = u(2*n+1:3*n);

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
% Derived trait-dependent quantities
ov  = ones(n,1);
Yi  = Ymax + chi * (Ymin - Ymax);
Ki  = Kmax + chi * (Kmin - Kmax);
Mi  = Mmax + chi * (Mmin - Mmax);

% === Partial derivatives for each equation ===
% Biomass equation (row 1)
f1bb = -2 * (Lam0 * w * Ki^2) ./ (b + Ki).^3;
f1wb =  (Lam0 * Ki^2) * ov ./ (b + Ki).^2;

% Water equation (row 2)
f2bb = -2 * Lam0 * R^2 * w ./ (1 + R*b).^3 ...
       + 2 * A * Q * Yi^2 * (1 - f) * h ./ (Q + Yi*b).^3;
f2bw = -Ga * ov + L0 * R * ov ./ (1 + R*b).^2;
f2bh = -A * Q * Yi * (1 - f) * ov ./ (Q + Yi*b).^2;

% Surface water equation (row 3)
f3bb = -2 * A * (1 - f) * Q * Yi^2 * h ./ (Q + Yi*b).^3;
f3bh =  A * (1 - f) * Q * Yi * ov ./ (Q + Yi*b).^2;

% === Scalar products with perturbation fields ===
ph1 = u(p.nu + 1         : p.nu +     n);
ph2 = u(p.nu + 1 +   n   : p.nu + 2 * n);
ph3 = u(p.nu + 1 + 2 * n : p.nu + 3 * n);

% === Construct diagonal matrices for Frechét derivatives ===
M1 = spdiags(f1bb .* ph1 + f2bb .* ph2 + f3bb .* ph3, 0, n, n);
M2 = spdiags(f1wb .* ph1 + f2bw .* ph2,                0, n, n);
M3 = spdiags(f2bh .* ph2 + f3bh .* ph3,                0, n, n);

% === Block matrix Fu construction ===
Z = sparse(n,n);  % Zero block for structure
Fu = [ M1  M2  M3;
       M2   Z   Z;
       M3   Z   Z ];

% === Compute matrix-vector product ===
M = p.mat.M;
duGuph = - Fu * M;




