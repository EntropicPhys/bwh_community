close all; keep p2phome;
global p2pglob;
p2pglob.gu = [];   % GUI off
p2pglob.nzi = 0;   % sparse jacobian (no info printing)
p2pglob.ps = 1;    % plot style

%% === Parameters and Initialization ===

% Domain and discretization
lx = 90;           % spatial domain length
nx = 350;          % number of grid points
N  = 33;           % number of trait bins (chi values)

% Initial condition selector
aux.solve = 1;
aux.sw    = 1;      % selects type of initial trait-distribution

% System parameters
pp = 300; Lam0 = 8; Ga = 10; A = 3000; L0 = 200; f = 0.01;
Q = 12; R = 0.7;
chimin = 0; chimax = 1;
Kmin = 6.7; Kmax = 35.599;
Mmin = 14.15; Mmax = 22.515;
Ymin = 0.069; Ymax = 0.11463;
Db = 1; Dw = 80; Dh = 1800; Dchi = 1e-4;

% Parameter vector for the model
par = [pp; Lam0; Ga; A; R; L0; f; Q;
       Kmin; Kmax; Mmin; Mmax; Ymin; Ymax;
       Db; Dw; Dh; Dchi; chimin; chimax];

% Output folders
dir0 = 'comm';
dir  = [dir0 '/hom'];

% Initialize problem structure
p = bwhinit(lx, nx, N, par, dir, aux);

%% === Time Integration to Relax to Steady State ===
% Useful for a good initial guess for Newton
t1 = 0; ts = [];
dt = 2e-3; nt = 4e4; nc = 0;
pmod = nt / 20; smod = pmod;

[p, t1, ts, nc] = tintxs(p, t1, ts, dt, nt, nc, pmod, smod, @sGdns);

%% === Newton Iteration to Find Steady State ===
p.nc.tol = 1e-11;
[p.u, p.r, iter] = nloop(p, p.u);
fprintf('res = %g, iter = %i\n', norm(p.r, Inf), iter);
plotsol(p);

%% === Homogeneous Branch Continuation ===
p.nc.ilam = 1;          % Index of parameter to continue (1 = precipitation)
p.sol.ds = -1;          % Initial continuation step size
p.sw.bifcheck = 2;      % Enable bifurcation detection
p.nc.tol = 1e-10;       % Tolerance for Newton step during continuation
p.sw.verb = 2;          % Verbosity level

tic;
p = cont(p, 30);        % Continue for 30 steps
toc;

%% === Branch Switching at Turing Point (bpt1) ===
aux = []; aux.m = 3; aux.besw = 0;
p0 = cswibra(dir, 'bpt1', aux);

% Configure eigenvalue and numerical tolerances
p0.nc.eigref = -3;
p0.nc.neig   = 3;
p0.nc.dsmin  = 1e-6;
p0.nc.dsmax  = 2;
p0.nc.tol    = 1e-7;

%% === Continue T1 Branch (First Turing Mode) ===
dirT = [dir0 '/T1'];
p = gentau(p0, [1]);     % Projection onto first mode
p = setfn(p, dirT);
p.sol.ds = 1e-1;
p = cont(p, 300);

%% === Continue T2 Branch (Second Turing Mode) ===
dirT = [dir0 '/T2'];
p = gentau(p0, [0 1]);   % Projection onto second mode
p = setfn(p, dirT);
p.sol.ds = 1e-1;
p = cont(p, 300);

%% === Plotting Bifurcation Branches ===
figure;

% Homogeneous solution branch (green)
plotbra('comm/hom', 'cl', [0 0.6 0], ...
        'tyun', '--', 'cmp', 3, 'bplab', 1, 'lab', 20);

% T1 (first pattern) (light blue)
plotbra('comm/T1', 'cl', [0 0.6 1], ...
        'tyun', '--', 'cmp', 3, 'bplab', [1 2], 'lab', 260);

% T2 (second pattern) (red)
plotbra('comm/T2', 'cl', 'r', ...
        'tyun', '--', 'cmp', 3, 'bplab', [1 2], 'lab', 260);

ylabel('$\langle \chi_{max} \rangle$', 'Interpreter', 'latex');
xlabel('$P$', 'Interpreter', 'latex');

%% === Plotting Solution Snapshots ===
% These correspond to specific solution points along the branches
plotsol('comm/hom', 'pt20',  'pfig', 1);  % Homogeneous state
plotsol('comm/T1',  'pt260', 'pfig', 2);  % First patterned state
plotsol('comm/T2',  'pt260', 'pfig', 3);  % Second patterned state
