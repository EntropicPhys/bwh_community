% cmds1.m - continuation of homogeneous and patterned solutions
% The parameter `chi` represents a phenotypic trait linked to plant functionality,
% such as water-use efficiency or rooting strategy. In the supplementary material,
% we explore two cases:
%   - chi=0: corresponds to fast-growing species
%   - chi=1: corresponds to stress-tolerant species
% These cases help illustrate how different species or functional types can affect
% pattern formation in the vegetation model.
close all; keep p2phome;
%% === Parameters ===
% Domain and discretization
pp=290;         % Precipitation (bifurcation parameter)
lx=90;  nx=300; % Length and discr of the spatial domain
% Model parameters
Lam0=8; Ga=10; A=3000; L0=200; f=0.01;  % Growth, mortality, infiltration, etc.
Q=12; R=0.7; chimin=0.0; chimax=1.0; Kmin=6.7; Kmax=28.93;
Mmin=14.15; Mmax=20.585; Ymin=0.069; Ymax=0.1041;
Db=1; Dw=80; Dh=1800; Dchi=1e-4; % Diffusion coefficients
% Trait and rescaling
chi=1; l=1;   % chi=fixed trait value; l=rescaling parameter
Yi=Ymin+chi*(Ymax-Ymin); % chi-dependent root-to-shoot ratio
% Parameter vector passed into model definition
par=[pp; Lam0; Ga; A; R; L0; f; Q; Kmin; Kmax; Mmin; Mmax; ...
       Ymin; Ymax; Db; Dw; Dh; chi; l];
b0=5.25; w0=3.15; h0=pp*(Yi*b0+Q)/(A*(Yi*b0+f*Q));  % initial guess 
dir0='s1'; dir=[dir0 '/hom']; % % Output folder, chi=1
%% === Initialize p-structure and run homogeneous branch ===
p=bwhinit(lx, nx, par, b0, w0, h0, dir);
p.plot.bpcmp=1;     % Component index to plot during continuation (1=b)
p.nc.eigref=-0.3;   % Eigenvalue reference for eigs routine
p.nc.neig=10;        % Number of eigenvalues to compute 
p.sol.ds=-0.1;     % Initial cont.step size (negative=backward in parameter)
p=cont(p,60);      % Continue for 60 steps along homogeneous branch
%% === Branch switching to Turing T1
p=swibra('s1/hom','bpt1','s1/T1',0.001); % branch-switching at 'bpt1' in s1/hom 
p.nc.neig=6;    % #eigenvalues to compute; here rather few for speed  
p.nc.eigref=-3;% ref.eigenval for eigs; here rather negative to not miss instab
p.nc.dsmin=1e-6;% min cont. step size; quite small since difficult folds 
p.nc.dsmax=5; p0.nc.dlammax=4;  % max arclength and parameter step sizes 
p.nc.tol=1e-9;      % tolerance for Newton's method during continuation
p=cont(p,200); % go! 
%% === Branch switching to Turing T2, T3 and T4 
p=swibra(dir,'bpt2','s1/T2',0.1); p.nc.neig=5;p.nc.eigref=-3;p.nc.dsmin=1e-6; 
p.nc.dsmax=5; p.nc.dlammax=4; p=cont(p,120); 
p=swibra(dir, 'bpt3','s1/T3',0.1); p.nc.neig=5; p.nc.eigref=-3; p.nc.dsmin=1e-6; 
p.nc.dsmax=5; p.nc.dlammax=4; p=cont(p,120); 
p=swibra(dir, 'bpt4','s1/T4',0.1); p.nc.neig=5; p.nc.eigref=-3; p.nc.dsmin=1e-6; 
p.nc.dsmax=5; p.nc.dlammax=4; p=cont(p,100); 
%% === Branch switching to Localized structures branch ===
p=swibra('s1/T1', 'bpt1','s1/T1s',0.1); p.nc.dsmax=10; p.nc.tol=1e-10; 
p.nc.dlammax=4; p.file.smod=50; % only save each 50st step 
p=cont(p,600); % long branch
%% DNS to k=0.28 branch,preparation 
p=loadp('s1/T1','pt170','du'); p.u(p.nu+1)=180;  % change P 
t1=0; dt=0.01; nt=50000; pmod=100; % #int.steps and plot each pmod-th step 
%% repeat until r=1e-6
[p,t1]=tint(p,t1,dt,nt,pmod); res=pderesi(p,p.u); r=norm(res,'inf') 
%% continue: reset counter and set filename to T2-double 
p=resetc(p); p=setfn(p,'s1/T2d'); p=cont(p,100); 
%% continue in other direction
p=loadp('s1/T2d','pt10','s1/T2db'); p.sol.ds=-p.sol.ds; p=cont(p,100); 
%% DNS to k=0.07 soln (single hump); but cannot continued due to transl.invar
p=loadp('s1/T2db','pt20','single'); p.u(p.nu+1)=155;  
p.u(1:round(p.np/4))=0; plotsol(p); title('one patch deleted'); pause
t1=0; dt=0.01; nt=10000; pmod=100; 
[p,t1]=tint(p,t1,dt,nt,pmod); res=pderesi(p,p.u); r=norm(res,'inf') 
fig(6); title(['t=' mat2str(t1,3)]); 
%% === Plot branches; see plotbra documentation for arguments of plotbra
cHom=[0.39,0.83,0.07]; cT1=[0.40,0.53,0.25]; % colors 
cT2=[0.64,0.81,0.43]; cT1s=[0.13,0.7,0.67]; cT028=[0.13,0.9,0];
c=1; ylab='||B||'; % select compo (see sgbra) and label to plot 
%c=3; ylab='max(B)'; % (uncomment as desired)
mclf(3); plotbra('s1/hom','cl',cHom,'tyun','--','cmp',c,'lab',10,'ln','hom');
plotbra('s1/T1','cl',cT1,'tyun','--','cmp',c,'lab',[100,170],'ln',{'A','B'});
plotbra('s1/T2','cl',cT2,'tyun','--','cmp',c,'lab',110,'ln','C');
plotbra('s1/T1s','cl',cT1s,'tyun','--','cmp',c,'lab',[100,400],'ln',{'F','G'});
plotbra('s1/T3','cl',p2pc('b5'),'tyun','--','cmp',c,'lab',80,'ln','D');
plotbra('s1/T2d','cl',p2pc('o1'),'tyun','--','cmp',c,'lab',[]);
plotbra('s1/T2db','cl',p2pc('o1'),'tyun','--','cmp',c,'lab',20,'ln','E');
xlabel('Precipitation'); ylabel(ylab); grid on
%% === soln plots
plotsol('s1/T1','pt100'); title('A'); pause; 
plotsol('s1/T1','pt170'); title('B'); pause;  
plotsol('s1/T2','pt110'); title('C'); pause; 
plotsol('s1/T3','pt80'); title('D'); pause
plotsol('s1/T2db','pt20'); title('E'); pause
plotsol('s1/T1s','pt100'); title('F'); pause; 
plotsol('s1/T1s','pt400'); title('G'); 
%% === further commands =========================
%% Branch switching to 2nd snake: 
p=swibra('s1/T1', 'bpt2','s1/T1s2',0.1); pause; p.nc.dsmax=10; p.nc.tol=1e-8; 
p.nc.dlammax=4; p.file.smod=50; % only save each 50st step 
p=cont(p,600); % long branch