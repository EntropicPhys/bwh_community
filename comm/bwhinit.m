function p=bwhinit(lx,nx,N,par,dir,aux) % bwh-com init function
p=[]; p=stanparam(p); % p-structure creation and basic init command
p=setfn(p,dir); p.fuha.outfu=@sgbra; % folder creation and output cont. file
pde=stanpdeo1Db(0,lx,lx/nx); % standard PDE object 1D
p.np=pde.grid.nPoints;p.pdeo=pde;n=p.np; % set array-struct dimensions
p.N=N;p.nc.neq=N+2;p.ndim=1;p.vol=lx;p.nu=p.np*p.nc.neq;p.sol.xi=0.1/p.nu;
ov=ones(n,1);b=zeros(n*N,1); % initial conditions switch
x=getpte(p); chimin=par(19); chimax=par(20); delchi=(chimax-chimin)/(N-1); 
switch aux.sw; 
    case 1; iv=1:N; 
        for i=iv; chi=chimin+(i-1)*delchi; 
            b((i-1)*n+1:i*n)=10*sech(28*(chi-0.825)).^2; end 
        w=0.9*ov; h=3.5*ov;     
    case 2; iv=1:N; 
        for i=iv; chi=chimin+(i-1)*delchi; 
            b((i-1)*n+1:i*n)=8*sech(28*(chi-0.36)).*cos((2*pi*10.8/p.vol)*x).^2; end
            %b((i-1)*n+1:i*n)=0.8*sech(27*(chi-0.88)); end
        w=15*ov; h=28*ov;
    case 3; iv=1:N; 
        for i=iv; chi=chimin+(i-1)*delchi; 
            b((i-1)*n+1:i*n)=0.5*sech(20*(chi-0.5)).^2; end
        w=35*ov; h=125*ov;
end 
p.u=[b; w; h; par]; % concatenation of bwh variables and system parameters
p.sw.sfem=-1; p=oosetfemops(p); % use OOPDE, generate FEM matrices
p.plot.pstyle=-1; % naturally, bwh requires special plots, via userplot and p2pglob.ps 
p.plot.bpcmp=2; p.plot.pcmp=2; % component branch and sol plotting
p.nc.ilam=1; % continue in par(1)
