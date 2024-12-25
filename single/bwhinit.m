function p=bwhinit(lx,nx,par,b0,w0,h0,dir) % bwh-single species init function
p=[]; p=stanparam(p);    % p-structure creation and basic init command
p=setfn(p,dir); % init. dir.
p.fuha.outfu=@sgbra;  % set output quantities from continuation
pde=stanpdeo1Db(0,lx,lx/nx);    % standard PDE objects 1D
p.pdeo=pde;p.vol=lx; p.np=pde.grid.nPoints;  
n=p.np;p.ndim=1;p.nc.neq=3;   % set array-struct dimensions
p.nu=p.np*p.nc.neq;p.sol.xi=0.1/p.nu; 
b=b0*ones(n,1);w=w0*ones(n,1);h=h0*ones(n,1); 
p.u=[b; w;h; par];              % init sol with parameters appended
p.sw.sfem=-1;p=oosetfemops(p); % use OOPDE, generate FEM matrices
p.sw.bifcheck=2;p.sw.jac=1;     % set bp-detection and numerical Jac 
p.nc.ilam=1; % continue in par(1)
p.sol.ds=0.01;p.nc.dsmin=0.01;p.nc.dsmax=3; %  set arc-length cont params
p.sw.qjac=0; % bif. point. cont. pure numeric.
