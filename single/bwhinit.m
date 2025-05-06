function p = bwhinit(lx, nx, par, b0, w0, h0, dir)
% bwhinit.m — Initialization function for the single-species BWH model
% Inputs:
%   lx  - Length of the spatial domain
%   nx  - Number of discretization points in space
%   par - Vector of model parameters (appended to state vector)
%   b0, w0, h0 - Initial homogeneous values for biomass, soil water, and surface water
%   dir - Directory name for output and problem labeling

%--- Initialize p-structure and set basic fields ---%
p = [];                       
p = stanparam(p);             % Create and fill base p-structure with default settings
p = setfn(p, dir);            % Set working directory

%--- Assign function handles ---%
p.fuha.outfu = @sgbra;        % Output function for bifurcation diagrams (e.g., biomass averages)

%--- Set up 1D finite element problem ---%
pde = stanpdeo1Db(0, lx, lx/nx);  % Generate standard 1D PDE object over [0, lx] with spacing dx
p.pdeo = pde;
p.vol  = lx;
p.np   = pde.grid.nPoints;

%--- Define system size and degrees of freedom ---%
p.ndim   = 1;         % 1D spatial domain
p.nc.neq = 3;         % Three PDE components: B, W, H
p.nu     = p.np * p.nc.neq;   % Total number of unknowns
p.sol.xi = 0.1 / p.nu;        % Initial arc-length predictor step size

%--- Set initial state ---%
b = b0 * ones(p.np, 1);
w = w0 * ones(p.np, 1);
h = h0 * ones(p.np, 1);
p.u = [b; w; h; par];         % Concatenate fields with parameter vector

%--- FEM operators and numerical settings ---%
p.sw.sfem = -1;               % Use OOPDE finite element interface
p = oosetfemops(p);           % Generate finite element matrices

%--- Continuation and bifurcation detection settings ---%
p.sw.bifcheck = 2;            % Enable bifurcation detection (secant method)
p.sw.jac      = 1;            % Use numerical Jacobian
p.nc.ilam     = 1;            % Continuation in first parameter (par(1))

%--- Arc-length continuation parameters ---%
p.sol.ds    = 0.01;           % Initial step size
p.nc.dsmin  = 0.01;           % Minimum step size
p.nc.dsmax  = 3;              % Maximum step size

%--- Bifurcation continuation setting ---%
p.sw.qjac = 0;                % Use numerical Jacobian for continuation from bifurcation points
end
