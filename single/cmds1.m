% cmds1.m — Continuation of homogeneous and patterned solutions
% This script computes the homogeneous solution branch and two patterned branches (T1, T2)
% for a 1D spatial vegetation model using pde2path. Assumes bwhinit.m defines the PDE.
% === Trait-based parameter: chi ===
% The parameter `chi` represents a phenotypic trait linked to plant functionality,
% such as water-use efficiency or rooting strategy. In the supplementary material,
% we explore two cases:
%   - chi = 0: corresponds to fast-growing species
%   - chi = 1: corresponds to stress-tolerant species
% These cases help illustrate how different species or functional types can affect
% pattern formation in the vegetation model.
close all; keep p2phome;
global p2pglob; p2pglob.ps = 1; % Set plot style (1 = classic, 2 = fancy)
%% === Parameters ===
% Domain and discretization
pp = 290;         % Precipitation (bifurcation parameter)
lx = 85;          % Length of the spatial domain
nx = 300;         % Number of finite element nodes in 1D
% Model parameters
Lam0 = 8; Ga = 10; A = 3000; L0 = 200; f = 0.01;  % Growth, mortality, infiltration, etc.
Q = 12; R = 0.7; chimin = 0.0; chimax = 1.0;
Kmin = 6.7; Kmax = 28.93;
Mmin = 14.15; Mmax = 20.585;
Ymin = 0.069; Ymax = 0.1041;
% Diffusion coefficients
Db = 1; Dw = 80; Dh = 1800; Dchi = 1e-4;
% Trait and rescaling
chi = 1; l = 1;  % chi = fixed trait value; l = rescaling parameter
Yi = Ymin + chi * (Ymax - Ymin);  % Chi-dependent root-to-shoot ratio
% Parameter vector passed into model definition
par = [pp; Lam0; Ga; A; R; L0; f; Q; Kmin; Kmax; Mmin; Mmax; ...
       Ymin; Ymax; Db; Dw; Dh; chi; l];
% Initial guess for homogeneous steady state (b, w, h)
b0 = 5.25; w0 = 3.15;
h0 = pp * (Yi * b0 + Q) / (A * (Yi * b0 + f * Q));  % Water balance formula
% Output folder
dir0 = 's1'; % usually linked with chi=1
dir = [dir0 '/hom'];
%% === Initialize p-structure and run homogeneous branch ===
p = bwhinit(lx, nx, par, b0, w0, h0, dir);
p.plot.bpcmp = 1;     % Component index to plot during continuation (1 = b)
p.nc.eigref = -0.3;   % Eigenvalue reference for eigs routine
p.nc.neig = 3;        % Number of eigenvalues to compute 
p.sol.ds = -0.1;      % Initial continuation step size (negative = backward in parameter)
p.nc.ilam = 1;        % Index of parameter to continue (1 = pp)
p = cont(p, 60);      % Continue for 60 steps along homogeneous branch
%% === checking user jac === 
[Gu,Gn]=jaccheck(p);Gd=abs(Gu-Gn); em=max(max(Gd)); 
mclf(20); spy(Gd>1e-3);
%% === Branch switching from Turing bifurcation ===
aux = []; aux.m = 3; aux.besw = 0;
p0 = cswibra(dir, 'bpt1', aux);  % Switch at bifurcation point 'bpt1' from s1/hom folder
p0.sw.bifcheck = 2;     % Switch-on bifurcation check type 2
p0.nc.neig = 5;         % Number of eigenvalues to compute 
p0.nc.eigref = -3;      % Eigenvalue reference for eigs routine
p0.sw.spcalc = 1;       % Compute spectrum each step
p0.nc.dsmin = 1e-6;     % Minimum continuation step size
p0.nc.dsmax = 2;       % Maximum continuation step size
p0.nc.tol = 1e-9;      % Tolerance for Newton's method during continuation
%% === Continue patterned branches (T1 and T2) ===
p=continue_branch(p0, [1],    [dir0 '/T1'], 1e-3, 200);  % T1 branch: 1st transverse mode
p=continue_branch(p0, [0, 1], [dir0 '/T2'], 1e-3, 200);  % T2 branch: 2nd transverse mode
%% === Branch switching to Localized structures branch ===
aux = []; aux.m = 3; aux.besw = 0;
p0 = cswibra([dir0 '/T1'], 'bpt2', aux); % Switch at bifurcation point 'bpt2' from s1/hom folder
p0.sw.bifcheck = 2;     % Switch-on bifurcation check type 2
p0.nc.neig = 5;         % Number of eigenvalues to compute 
p0.nc.eigref = -3;      % Eigenvalue reference for eigs routine
p0.sw.spcalc = 1;       % Compute spectrum each step
p0.nc.dsmin = 1e-6;     % Minimum continuation step size
p0.nc.dsmax = 10;       % Maximum continuation step size
p0.nc.tol = 1e-10;      % Tolerance for Newton's method during continuation
%% === Continue Localized patterns branches (T1 and T2) ===
p=continue_branch(p0, [1],[dir0 '/T1s'], 1e-3, 800);  % T1s branch: LPs 
%% === Plot all branches with different green shades ===
% === Colors for each branch ===
colHom = [0.39, 0.83, 0.07];
colT1 = [0.40, 0.53, 0.25];
colT2 = [0.64, 0.81, 0.43];
colT1s = [0.13, 0.7, 0.67];
box on
% Plot branches
clf("reset")
plotbra([dir0 '/T1s'],  'cl', colT1s,  'tyun', '--', 'cmp', 1,'lab',420);
plotbra([dir0 '/hom'], 'cl', colHom, 'tyun', '--', 'cmp', 1,'lab',20);
plotbra([dir0 '/T1'],  'cl', colT1,  'tyun', '--', 'cmp', 1,'lab',180);
plotbra([dir0 '/T2'],  'cl', colT2,  'tyun', '--', 'cmp', 1,'lab',180);
xlabel('Precipitation'); ylabel('||b||');
% Example positions 
text(270,4.9,'Hom','Color',colHom,'FontSize',11,'FontWeight','bold');
text(344,5.5,'Turing 1','Color',colT1,'FontSize',11,'FontWeight','bold');
text(310,4.7,'Turing 2','Color',colT2,'FontSize',11,'FontWeight','bold');
text(298,6.2,'Localized pattern','Color', colT1s,'FontSize',11,'FontWeight','bold');
axis([250 360 3.5 7])
%% === Plotting Solution Snapshots ===
% These correspond to specific solution points along the branches
plotsol('s1/hom','pt20', 'pfig',1);  axis([0 85 0 9]); pause % Homogeneous state
plotsol('s1/T1','pt180','pfig',2);  axis([0 85 0 9]); pause % First patterned state
plotsol('s1/T2','pt180','pfig',3);  axis([0 85 0 9]); pause % Second patterned state
plotsol('s1/T1s','pt420','pfig',4); axis([0 85 0 9]); pause % Localized patter state


%% === Helper Function ===
function p=continue_branch(p0, tauvec, dirname, ds, nsteps)
    % Helper to switch to patterned branch and run continuation
    % Inputs:
    %   p0     = p-structure from cswibra (at bifurcation point)
    %   tauvec = mode selector for gentau (e.g. [1], [0 1], etc.)
    %   dirname= output folder
    %   ds     = step size
    %   nsteps = number of continuation steps
    p = gentau(p0, tauvec);     % Create new branch with transverse mode(s)
    p = setfn(p, dirname);      % Set filename prefix
    p.sol.ds = ds;              % Set step size
    p=cont(p, nsteps);            % Continue
end
