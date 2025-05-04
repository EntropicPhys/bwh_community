% bwhinit.m — Initialization function for the BWH community model
% This function sets up a pde2path problem structure `p` for a 1D trait-structured
% biomass–water–surface water model, using the OOPDE finite element framework.
% It defines spatial and trait discretizations, boundary conditions, and initial states.

function p = bwhinit(lx, nx, N, par, dir, aux)

% === Create base p-structure and assign function handles ===
p = [];
p = stanparam(p);                     % Initialize standard p-structure
p = setfn(p, dir);                    % Set directory for results
p.fuha.outfu = @sgbra;               % Output function for branch data

% === Define PDE domain (1D) and FEM discretization ===
pde = stanpdeo1Db(0, lx, lx / nx);   % Standard 1D pdeo object over [0, lx] with nx elements
p.np = pde.grid.nPoints;             % Number of spatial grid points
p.pdeo = pde;
n = p.np;                            % Alias for readability

% === Define system and trait discretization ===
p.N = N;                             % Number of trait discretization bins
p.nc.neq = N + 2;                    % Total number of equations (N traits + w + h)
p.ndim = 1;                          % Spatial dimension
p.vol = lx;                          % Domain length
p.nu = p.np * p.nc.neq;              % Total number of degrees of freedom
p.sol.xi = 0.1 / p.nu;               % Arclength continuation parameter

% === Trait discretization and initial condition setup ===
ov = ones(n, 1);                     % Vector of ones for spatial profile
b = zeros(n * N, 1);                 % Initialize biomass component
x = getpte(p);                       % Get spatial node positions

chimin = par(19); chimax = par(20);
delchi = (chimax - chimin) / (N - 1);

switch aux.sw
    case 1  % Localized phenotype near chi = 0.825
        for i = 1:N
            chi = chimin + (i - 1) * delchi;
            b((i - 1) * n + 1 : i * n) = 10 * sech(28 * (chi - 0.825)).^2;
        end
        w = 0.9 * ov;
        h = 3.5 * ov;

    case 2  % Periodic phenotype structure near chi = 0.36
        for i = 1:N
            chi = chimin + (i - 1) * delchi;
            b((i - 1) * n + 1 : i * n) = 8 * sech(28 * (chi - 0.36)) .* cos((2 * pi * 10.8 / p.vol) * x).^2;
        end
        w = 15 * ov;
        h = 28 * ov;

    case 3  % Weak phenotype expression centered at chi = 0.5
        for i = 1:N
            chi = chimin + (i - 1) * delchi;
            b((i - 1) * n + 1 : i * n) = 0.5 * sech(20 * (chi - 0.5)).^2;
        end
        w = 35 * ov;
        h = 125 * ov;
end

% === Assemble full solution vector and finalize FEM setup ===
p.u = [b; w; h; par];                % Initial state: biomass (vector), w, h, and parameters
p.sw.sfem = -1;                      % Use OOPDE FEM engine
p = oosetfemops(p);                  % Assemble FEM matrices

% === Plotting and continuation parameters ===
p.plot.pstyle = -1;                 % Use user-defined plotting style
p.plot.bpcmp = 2;                    % Branch plotting component (e.g. water)
p.plot.pcmp = 2;                     % Solution plotting component
p.nc.ilam = 1;                       % Index of parameter to continue (e.g. precipitation)
end
