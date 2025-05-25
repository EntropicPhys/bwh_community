function p=bwhinit(lx,nx,N,par,dir,aux) 
% more or less standard, except maybe l10 where we need to init N b-fields
p=stanparam; p.ndim=1; p=setfn(p,dir); p.fuha.outfu=@sgbra; 
p.N=N; p.nc.tol=1e-8; p.file.smod=10; p.nc.ilam=1; p.sw.verb=2; 
p.nc.neq=N+2; p.plot.bpcmp=3; p.nc.lammin=1; 
p.nc.lammax=500; p.sw.sfem=-1; p.nc.neig=10; 
pde=stanpdeo1Db(0,lx,lx/nx); p.vol=lx; 
p.np=pde.grid.nPoints; p.pdeo=pde; n=p.np; ov=ones(n,1);b=zeros(n*N,1); 
x=getpte(p); chimin=par(19); chimax=par(20); delchi=(chimax-chimin)/(N-1); 
switch aux.sw; 
    case 1; b=0.2*ov; w=20*ov; h=4*ov; b=repmat(b,N,1); % homogen.IC      
    case 2; b=0.3*ov; w=1*ov; h=3*ov; b=repmat(b,N,1); % homogen.IC 
    case 3; iv=1:N; 
        for i=iv; chi=chimin+(i-1)*delchi; 
            b((i-1)*n+1:i*n)=10*sech(28*(chi-0.825)).^2; end 
        w=0.9*ov; h=3.5*ov;     
    case 4; iv=1:N; 
        for i=iv; chi=chimin+(i-1)*delchi; 
            b((i-1)*n+1:i*n)=8*sech(28*(chi-0.36)).*cos((2*pi*10.8/p.vol)*x).^2; end       
        w=15*ov; h=28*ov;
    case 5; iv=1:N; 
        for i=iv; chi=chimin+(i-1)*delchi; 
            b((i-1)*n+1:i*n)=0.5*sech(20*(chi-0.5)).^2; end
        w=35*ov; h=125*ov;    
end 
p.u=[b; w; h; par]; 
p.nu=p.np*p.nc.neq; p.sw.jac=0; p=oosetfemops(p); 
p.nc.dsmin=1e-6; p.nc.dsmax=10; p.nc.dlammax=10; 
p.sol.xi=0.1/p.nu; p.sw.spcalc=0; p.sw.bifcheck=2; p.nc.imax=10; 
p.plot.pstyle=-1; % naturally, bwh requires special plots, via userplot and p2pglob.ps 
p.r=resi(p,p.u); fprintf('inires=%g\n',norm(p.r,Inf)); 
p.plot.pcmp=2;
p.sw.spcalc=1; p.sol.ds=-0.1; p.nc.nsteps=50; p.frcut=0.0001; p.nc.mu1=0.2; 