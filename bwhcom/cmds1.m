%% bhwCom with pars from [FPBUM25], 
keep p2phome; global p2pglob; p2pglob.gu=[]; p2pglob.nzi=0; p2pglob.ps=2;
%% folder
dir0='1';  N=35; lx=90; nx=250; 
aux.sw=5; aux.solve=1; % switches to control bwhinit (see there)  
Lam0=8; A=3000; Ga=10; L0=200; f=0.01; pp=350; 
Q=12; R=7e-1;chimin=0.0; chimax=1; Kmin=6.7; Kmax=35.599; 
Mmin=14.15; Mmax=22.515; Ymin=0.069; Ymax=0.11463; 
Db=1; Dw=80; Dh=1800; Dchi=1e-4; 
par=[pp; Lam0; Ga; A; R; L0; f; Q; Kmin; Kmax; Mmin; Mmax; Ymin; Ymax; Db; Dw; Dh; Dchi; ...
    chimin; chimax];
dir=char([dir0 '/hom']);
p=bwhinit(lx,nx,N,par,dir,aux); p=setfn(p,dir);plotsol(p);
%% DNS to get near spatially homogeneous steady state 
t1=0; ts=[]; dt=2e-3; nt=6e4; nc=0; pmod=nt/210; smod=pmod; 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@sGdns); 
%% Newton loop for steady state 
p.nc.almin=0.25;p.nc.tol=1e-10; [p.u,p.r,iter]=nloop(p,p.u); 
fprintf('res=%g, iter=%i\n',norm(p.r,Inf),iter); plotsol(p); 
%% continue hom branch 
p.sw.para=2; p.nc.eigref=-1e-1; p.nc.tol=3e-10; p.nc.neig=5; 
tic; p=cont(p,15); toc % large steps 
p.sol.ds=-0.1; p.nc.dsmax=2; p=cont(p,20); %smaller steps for BP detection
%% swibra to first Turing branch T1
p=swibra('1/hom','bpt1','1/T1',0.001); p.nc.dsmin=1e-6; p=cont(p,10); 
p.sw.bifcheck=0; p.nc.eigref=-4; p.nc.dsmax=10; p.nc.neig=2; p=cont(p,100);
%% swibra to primary snake of localized patterns 
p=swibra('1/T1','bpt1','1/T1-1',0.001); pause; p.nc.tol=1e-8; p=cont(p,10); 
p.file.smod=50; p.sw.bifcheck=0; p.nc.dsmax=5; p.nc.eigref=-4; p=cont(p,1000); 
%% swibra to T2
p=swibra('1/hom','bpt2','1/T2',0.001); p.nc.dsmin=1e-6; p=cont(p,10); 
p.sw.bifcheck=0; p.nc.eigref=-4; p.nc.dsmax=10; p.nc.neig=2; p=cont(p,100);
%% cswibra to T2C, k=0.66 
aux=[]; aux.m=4; aux.besw=0;p0=cswibra('1/hom','bpt2',aux); p0.sw.spcalc=1; 
p0.sw.bifcheck=2; p0.nc.eigref=-3; p0.file.smod=20; p0.nc.neig=4;
p0.nc.dsmin=1e-4; p0.nc.tol=1e-9; p0.nc.bisecmax=10; p0.nc.dsmax=1;
%% gentau
p=gentau(p0,[0 0 1]); p=setfn(p,'1/T2C'); p.sol.ds=0.1; plotsol(p); 
p.sw.bifcheck=0; p.nc.eigref=-4; p.nc.dsmax=10; p.nc.neig=2; p=cont(p,100); 
%% cswibra to T3D, k=0.84
aux=[]; aux.m=5; aux.besw=0;p0=cswibra('1/hom','bpt3',aux); p0.sw.spcalc=1; 
p0.sw.bifcheck=2; p0.nc.eigref=-3; p0.file.smod=20; p0.nc.neig=4;
p0.nc.dsmin=1e-4; p0.nc.tol=1e-9; p0.nc.bisecmax=10; p0.nc.dsmax=1;
%% gentau; inspect kernel vectors here, and choose the one with largest k
p=gentau(p0,[1]); p=setfn(p,'1/T3D'); p.sol.ds=0.1; plotsol(p); 
p.sw.bifcheck=0; p.nc.eigref=-4; p.nc.dsmax=10; p.nc.neig=2; p=cont(p,100); 
%% swibra to T4 (mixed) 
p=swibra('1/hom','bpt4','1/T4',-0.001); p.nc.dsmin=1e-6; p.nc.tol=1e-6;  p=cont(p,10); 
p.sw.bifcheck=0;p.nc.eigref=-4;p.nc.dsmax=10;p.nc.neig=2; p=cont(p,80);
%% cswibra to T4B, k=0.91, but unstable throughout
aux=[]; aux.m=4; aux.besw=0;p0=cswibra('1/hom','bpt4',aux); p0.sw.spcalc=1; 
p0.sw.bifcheck=2; p0.nc.eigref=-3; p0.file.smod=20; p0.nc.neig=4;
p0.nc.dsmin=1e-4; p0.nc.tol=1e-9; p0.nc.bisecmax=10; p0.nc.dsmax=1;
%% gentau 
p=gentau(p0,[0,1]); p=setfn(p,'1/T4B'); p.sol.ds=0.1; plotsol(p); pause 
p.sw.bifcheck=0; p.nc.eigref=-4;p.nc.dsmax=10;p.nc.neig=2;p=cont(p,100); 
%% "delete patches" to go to other k via DNS 
p=loadp('1/T1','pt70','du'); b=p.u(1:p.np*p.N); M=40; % how much to delete
for i=1:p.N; b((i-1)*p.np+p.np-M:i*p.np) = 0; end
p.u(1:p.N*p.np)=b; w=p.u(p.N*p.np+1: (p.N+1)*p.np); wmin=min(w); 
w(p.np-M:end)=wmin; p.u(p.N*p.np+1: (p.N+1)*p.np)=w;
h=p.u((p.N+1)*p.np+1: p.np*(p.N+2)); hmax=max(h); 
h(p.np-M:end)=hmax; p.u((p.N+1)*p.np+1: p.np*(p.N+2))=h; plotsol(p); 
% DNS
t1=0;  ts=[]; dt=0.5e-2; nt=1e5; nc=0; pmod=nt/200; smod=pmod; 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@sGdns); 
%% Newton loop for steady state 
p.nc.almin=0.25;p.nc.tol=1e-10; [p.u,p.r,iter]=nloop(p,p.u); 
fprintf('res=%g, iter=%i\n',norm(p.r,Inf),iter); plotsol(p); 
%% cont steady state 
p=setfn(p,'1/T1b'); p=resetc(p); stansavefu(p); p.sol.ds=-0.01; p=cont(p,41); 
%% reverse direction 
p=loadp('1/T1b','pt0','1/T1bb'); p=resetc(p); p.sol.ds=-p.sol.ds; p=cont(p,41); 
%% once more: delete patches to go to longer wavelength via DNS 
p=loadp('1/T2','pt110','du'); b=p.u(1:p.np*p.N); M=40;
for i=1:p.N; b((i-1)*p.np+p.np-M:i*p.np) = 0; end
p.u(1:p.N*p.np)=b; w=p.u(p.N*p.np+1: (p.N+1)*p.np); 
w(p.np-M:end)=wmin; p.u(p.N*p.np+1: (p.N+1)*p.np)=w;
h=p.u((p.N+1)*p.np+1: p.np*(p.N+2)); hmin=min(h); hmax=max(h); 
h(p.np-M:end)=hmax; p.u((p.N+1)*p.np+1: p.np*(p.N+2))=h; plotsol(p); 
%% DNS
t1=0;  ts=[]; dt=0.5e-2; nt=1e5; nc=0; pmod=nt/200; smod=pmod; 
[p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@sGdns); 
%% Newton loop for steady state 
p.nc.almin=0.25;p.nc.tol=1e-10;
[p.u,p.r,iter]=nloop(p,p.u); fprintf('res=%g, iter=%i\n',norm(p.r,Inf),iter); plotsol(p); pi=p; 
%% cont steady state 
p=setfn(p,'1/T2b'); p=resetc(p); stansavefu(p); p.sol.ds=-0.01; p=cont(p,41); 
%% reverse direction 
p=loadp('1/T2b','pt0','1/T2bb'); p=resetc(p); p.sol.ds=0.5; p=cont(p,3); 
%% plot 
fn=3; mclf(fn); cmp=3; ylab='$\chi_{\max}$'; yl=[0.38 0.85]; % for plotting chi_max 
% cmp=1; ylab='$\|B\|_1$'; yl=[0.02 0.28];   % uncomment for plotting biomass 
plotbra('1/hom','cl',p2pc('g1'),'cmp',cmp,'lab',0,'ln',{'hom'})
plotbra('1/T1','cl',p2pc('g2'),'cmp',cmp,'lab',70,'ln','A');
plotbra('1/T1b','cl',p2pc('v1'),'cmp',cmp,'lab',0,'ln','C');
plotbra('1/T1bb','cl',p2pc('v1'),'cmp',cmp); 
plotbra('1/T2b','cl',p2pc('v3'),'cmp',cmp,'lab',0,'ln','D');
plotbra('1/T2bb','pt20','cl',p2pc('v3'),'cmp',cmp); 
plotbra('1/T2','cl',p2pc('b4'),'cmp',cmp,'lab',110,'ln','B'); 
plotbra('1/T2C','cl',p2pc('r2'),'cmp',cmp,'lab',60,'ln','E'); 
plotbra('1/T3D','cl',p2pc('b2'),'cmp',cmp); 
plotbra('1/T1-1','cl',p2pc('o1'),'cmp',cmp,'lab',[500 600 1000],'ln',{'F','G','H'}); 
ylabel(ylab,'Interpreter','latex'); xlabel('P'); grid on; box on; 
ylim(yl); 
%%
p2pglob.ps=3;  % plot-style 
plotsol('1/hom','pt0'); title('hom'); pause 
plotsol('1/T1','pt70'); title('A'); pause; plotsol('1/T2','pt90');title('B');  pause; 
plotsol('1/T1b','pt0'); title('C'); pause; plotsol('1/T2b','pt0'); title('D');  pause 
plotsol('1/T3D','pt60'); title('E'); pause; plotsol('1/T1-1','pt500'); title('F');pause; 
plotsol('1/T1-1','pt600');  title('G'); pause; plotsol('1/T1-1','pt1400');  title('H'); 