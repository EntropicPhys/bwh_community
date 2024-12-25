function [p,t1,ts,nc]=tintxs(p,t0,ts,dt,nt,nc,pmod,smod,nffu,varargin)
% TINTXS: time integration with time-series output, LU-decomp. of M+dt*K, 
%  Here modified to only treat x-diffusion implicitly, chi-diff expl.
par=p.u(p.nu+1:end); Db=par(15); Dw=par(16); Dh=par(17); N=p.N; n=p.np; ng=N*n; u0=p.u; 
r=norm(resi(p,p.u),'inf'); ts=[ts [t0;r]]; % put time and residual into ts
% if ~exist([p.file.dir '/pt0.mat'],'file'); p.t=0; 
%     r=norm(resi(p,p.u),'inf'); p.ts=[p.t;r]; 
%     p.file.count=0; p.sol.ptype=0;  p.fuha.savefu(p); 
% end
%noa=nargin-9; 
K=p.mat.K; J=sparse(ng); 
for i=1:N; J((i-1)*n+1:i*n,(i-1)*n+1:i*n)=Db*K; end  % diffusion B 
J(N*n+1:(N+1)*n,N*n+1:(N+1)*n)=Dw*K; % diffusion w
J((N+1)*n+1:(N+2)*n,(N+1)*n+1:(N+2)*n)=Dh*K; % diffusion h 
Lam=p.mat.M+dt*J; % mclf(20); spy(Lam); pause 
[L,U,P,Q,R]=lu(Lam); t=t0; n=0; 
while(n<nt) % integration loop
  f=nffu(p,p.u);% f=p.mat.M*f; 
  g=p.mat.M*p.u(1:p.nu)+dt*f; 
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); t=t+dt;  n=n+1;
  if(mod(n,pmod)==0); 
      r=norm(resi(p,p.u),'inf'); ts=[ts [t;r]]; % put time and residual into ts
      tits=['t=' mat2str(t,4) ', r=' mat2str(r,3)];
      plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); 
      title(['u_1, ' tits],'fontsize',12); set(gca,'fontsize',12); 
      drawnow; 
      %max(abs(u0(1:p.nu)-p.u(1:p.nu)));
  end
  % if(mod(n,smod)==0); 
  %   ps=p; p=[]; p.t=t;p.ts=ts; p.u=ps.u; 
  %   fname=[ps.file.pname,sprintf('%i',nc+n),'.mat']; save(fname,'p'); 
  %   p=ps;
  % end
end 
t1=t; nc=nc+nt; 