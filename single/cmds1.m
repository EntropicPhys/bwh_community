close all; keep p2phome; % script for 1D 
global p2pglob; p2pglob.ps=1;   % plotstyle 
%% parameters and folder
pp=300;lx=85; nx=300; % domain parameters
Lam0=8;Ga=10; A=3000;L0=200; f=0.01;  % system params.
Q=12; R=0.7;chimin=0.0; chimax=1;
Kmin=6.7; Kmax=28.93; Mmin=14.15; Mmax=20.585; Ymin=0.069; Ymax=0.1041;
Dchi=1e-4;Db=1; Dw=80; Dh=1800;
chi=1;  % \chi value assigned manually
l=1;   % fictitious spatial scale param.
par=[pp; Lam0; Ga; A; R; L0; f; Q; Kmin; Kmax; Mmin; Mmax; Ymin; Ymax; Db; Dw; Dh; chi;l]; % params vector
b0=5.25;w0=3.15;h0=pp*(Yi*b0+Q)/(A*(Yi*b0+f*Q)); % initial homog. sol. guess
dir0='aaa';dir=char([dir0 '/hom']); % principal and homog. sol directories
p=bwhinit(lx,nx,par,b0,w0,h0,dir0); % p-struct init

%% homog. solution continuation
p.plot.bpcmp=1; % component plotted during continuation
p.nc.eigref=-0.3;p.nc.neig=3; % eigenval. reference and number of eigenvals.
p.sol.ds=-0.01; p.nc.ilam=1; % initial cont. step and 1 param. cont. 
p=cont(p,10); % continuation

%% cont. from Turing using cswibra
aux=[]; aux.m=3; aux.besw=0; p0=cswibra(dir,'bpt1',aux); p0.sw.bifcheck=2;
p0.nc.neig=5; p0.nc.eigref=-0.5; p0.sw.spcalc=1;
p0.nc.dsmin=1e-6; p0.nc.dsmax=5; p.nc.tol=1e-11;

%% T1 
dirT=char([dir0 '/T1']);p=gentau(p0,[1]); % 6
p=setfn(p,dirT);p.sol.ds=1e-3;p=cont(p,10);

%% T2 
dirT=char([dir0 '/T2']);p=gentau(p0,[0 1]); % 6
p=setfn(p,dirT);p.sol.ds=1e-3;p=cont(p,10);


%%
plotbra(char([dir0 '/hom']),'cl','b','tyun','--','cmp',1)
plotbra(char([dir0 '/T1']),'cl','b','tyun','--','cmp',1)
plotbra(char([dir0 '/T2']),'cl','b','tyun','--','cmp',1)
