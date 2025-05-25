function p=bwhinit(lx, nx, par, b0, w0, h0, dir)
% bwhinit.m - Initialization function for the single-species BWH model
%   lx   - domain length
%   nx   - number of spatial discretization points
%   par  - vector of model parameters
%   b0, w0, h0 - initial hom. values for biomass (B), soil water (W), and surface water (H)
%   dir  - directory name for output and problem labeling
p=stanparam;              % Initialize p-structure with default settings
p.fuha.outfu=@sgbra;        % Custom output function for continuation diagnostics
p=setfn(p, dir);            % Set output folder and problem ID
p.file.smod=10;             % save each 10th step to disk 
%--- PDE discretization ---%
pde=stanpdeo1Db(0, lx, lx/nx);  % Create 1D PDE object over [0, lx] with spacing dx=lx/nx
p.pdeo=pde; 
p.vol =lx;                      % Domain volume (length in 1D)
p.np  =pde.grid.nPoints;        % Number of spatial nodes
p.ndim   =1;            % One spatial dimension
p.nc.neq =3;            % Three coupled PDEs: B, W, H
p.nu     =p.np * p.nc.neq;       % Total number of unknowns
p.sol.xi =0.1 / p.nu;            % Predictor step size in continuation
n=p.np;  b=b0 * ones(n, 1);  %--- Initial state vector (homogeneous) ---%
w=w0 * ones(n, 1); 
h=h0 * ones(n, 1); 
p.u=[b; w; h; par];             % Concatenate state & parameter vectors
p.sw.sfem=-1;                  % Use OOPDE-based FEM (standard in pde2path)
p=oosetfemops(p);              % Generate and assign FEM matrices
%--- Bifurcation and continuation settings ---%
p.nc.ilam=1;        % Index of parameter to continue (1=pp)
p.sw.bifcheck=2;               % Use secant-based bifurcation detection
p.sw.jac     =1;               % Use numerical Jacobian for continuation
p.nc.ilam    =1;               % First parameter in 'par' is active for continuation
%--- Arc-length continuation step control ---%
p.sol.ds   =0.01;              % Initial step size
p.nc.dsmin =0.01;              % Minimum step size
p.nc.dsmax =3;                 % Maximum step size
%--- plotting
p.plot.pcmp=[1 2]; p.plot.lw=[2 2]; 
g1=p2pc('g2'); b1=[0 0 0.9]; p.plot.cl=[g1;b1]; 
p.sw.verb=2; % verbosity 
end
