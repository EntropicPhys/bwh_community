function out = sgbra(p, u)
% sgbra — Custom branch output function for bwh-community model
%
% Output vector:
%   out(1): L1-norm of B (average biomass)
%   out(2): L2-norm of B (energy-like biomass norm)
%   out(3): chi_center — center of mass of trait distribution
%   out(4): dummy (set to 1, for plotting hacks)
%   out(5): bmax — max value of biomass
%   out(6:end): parameters (for diagnostics)
%
% This function is called by pde2path's continuation framework to record
% custom outputs for bifurcation diagrams via `plotbra`.

% === Extract model structure and parameters ===
par = p.u(p.nu+1:end);       % Parameter vector
u   = u(1:p.nu);             % State vector only
n   = p.np;
N   = p.N;
ng  = n * N;
M   = p.mat.M(1:ng, 1:ng);   % Mass matrix for all B components
vol = p.vol * N;             % Total volume across trait dimension

% === Compute biomass norms ===
l1v = sum(M * u(1:ng)) / vol;                              % L1 norm
l2v = sqrt(u(1:ng)' * (M * u(1:ng)) / vol);                % L2 norm

% === Compute trait-weighted biomass metrics ===
M1  = p.mat.M(1:n, 1:n);         % Mass matrix for 1 trait slice
vol1 = p.vol;
b = p.u(1:ng);
B = reshape(b, n, N);            % B(x, chi)
chi = linspace(par(19), par(20), N);   % Trait grid (chi)
ba = zeros(1, N);                % Integrated biomass over space
baa = zeros(1, N);               % Chi-weighted biomass for CoM

% === Integrate biomass over space for each trait ===
for i = 1:N
    ba(i)  = sum(M1 * B(:, i)) / vol1;
    baa(i) = chi(i) * ba(i);    % Weighted for CoM
end

% === Center of mass in trait space ===
chicent = sum(baa) / sum(ba);

% === Interpolate biomass over finer trait grid to find bmax more smoothly ===
chii = linspace(0, 1, 700);              % Fine trait grid
interbiochi = interp1(chi, ba, chii, 'spline'); % Interpolated biomass
[biomax, ~] = max(interbiochi);          % Max biomass across traits

% === Optional: Trait support detection (bounding chi where biomass > cutoff) ===
try
    frcut = p.frcut;   % optional threshold value
catch
    frcut = 0.001;
end

above_cutoff = find(ba > frcut);
if ~isempty(above_cutoff)
    brmin = chi(min(above_cutoff));
    brmax = chi(max(above_cutoff));
else
    brmin = NaN;
    brmax = NaN;
end

% === Assemble final output vector ===
out = [l1v; l2v; chicent; 1; biomax; par];  % Final output (used by plotbra)
