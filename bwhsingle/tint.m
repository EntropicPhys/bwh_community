function [p,t1]=tint(p,t1,dt,nt,pmod) % mod of library function for DNS for 
% semilinear RD systems. Linearly implicit, with prefactored Lam 
% In : problem struct p as usual, t1=initial time, dt=stepsize, nt=#steps, 
%      pmod=plotting interval (plot every pmod-th step)
% Out: (updated) p and t1
n=0; t=t1; Db=p.u(p.nu+15); Dw=p.u(p.nu+16); Dh=p.u(p.nu+17); 
K=kron(diag([Db,Dw,Dh]),p.mat.K); % diffusion matrix 
Lam=p.mat.M+dt*K; [L,U,P,Q,R]=lu(Lam); % prefactor stepping matrix 
while(n<nt) % integration loop
  f=nodalf(p,p.u);   % the nonlinearity, i.e., everything except diffusion
  g=p.mat.M*p.u(1:p.nu)+dt*f;    
  p.u(1:p.nu)=Q*(U\(L\(P*(R\g)))); n=n+1; t=t+dt; % time stepping:
  if(mod(n,pmod)==0); plotsol(p,p.plot.ifig,p.plot.pcmp,p.plot.pstyle); end
end 
t1=t;