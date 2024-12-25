close all; keep p2phome; 
global p2pglob; p2pglob.gu=[]; p2pglob.nzi=0;p2pglob.ps=1;
%% parameters and folder
lx=90; nx=350; N=33;  % domain parameters
aux.solve=1;  aux.sw=1; % init solution
pp=300; Lam0=8;Ga=10; A=3000;L0=200; f=0.01; % system params.
Q=12; R=0.7;chimin=0; chimax=1;
Kmin=6.7; Kmax=35.599; Mmin=14.15; Mmax=22.515; Ymin=0.069; Ymax=0.11463;
Dchi=1e-4;Db=1; Dw=80; Dh=1800;
par=[pp; Lam0; Ga; A; R; L0; f; Q; Kmin; Kmax; Mmin; Mmax; Ymin; Ymax; Db; Dw; Dh; Dchi;chimin; chimax]; % system params vector.
dir0='comm'; dir=char([dir0 '/hom']); % principal and homog. sol. directory
p=bwhinit(lx,nx,N,par,dir,aux); % p-struct initialization

%% time integration for initial cont. guess
t1=0;  ts=[]; dt=2e-3; nt=4e4; nc=0; pmod=nt/20; smod=pmod; 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@sGdns); 
%% newton-iteration 
p.nc.tol=1e-11;
[p.u,p.r,iter]=nloop(p,p.u); fprintf('res=%g, iter=%i\n',norm(p.r,Inf),iter); plotsol(p); 
%% homogeneous continuation branch.
p.nc.ilam=1; p.sol.ds=-1; % cont. param. and init. continuation step
p.sw.bifcheck=2;p.nc.tol=1e-10; % bif. detection criteria and residual tol.
p.sw.verb=2;
tic; p=cont(p,30); toc

%% cont. from Turing using cswibra
aux=[]; aux.m=3; aux.besw=0;  p0=cswibra(dir,'bpt1',aux); % basic setting
p0.nc.eigref=-3;p0.nc.neig=3; % ref. and number of eigval to calc.  
p0.nc.dsmin=1e-6; p0.nc.dsmax=2; % min. and max. cont. step size 
p0.nc.tol=1e-7; % residual tolerance.

%% T1
dirT=char([dir0 '/T1']);p=gentau(p0,[1]); p=setfn(p,dirT); % dir and comp. projection to cont.
p.sol.ds=1e-1; p=cont(p,300); % cont. step size and cont.

%% T2
dirT=char([dir0 '/T2']); p=gentau(p0,[0 1]); p=setfn(p,dirT); % dir and comp. projection to cont.
p.sol.ds=1e-1; p=cont(p,300); % cont. step size and cont.

%% plotting solution branches
plotbra('comm/hom','cl',[0 0.6 0],'tyun','--','cmp',3,'bplab',1,'lab',20) % homog. solution branch
plotbra('comm/T1','cl',[0 0.6 1],'tyun','--','cmp',3,'bplab',[1 2],'lab',260) % First Turing pattern
plotbra('comm/T2','cl','r','tyun','--','cmp',3,'bplab',[1 2],'lab',260) % Second Turing pattern
ylabel('$\langle \chi_{max} \rangle$','Interpreter','latex');xlabel('$P$','Interpreter','latex')

%% plotting solution states
plotsol('comm/hom','pt20','pfig',1);
plotsol('comm/T1','pt260','pfig',2);
plotsol('comm/T2','pt260','pfig',3);
