function out=sgbra(p,u)
% branch output for the bwh model; here ||B||_1, ||B||_2, chimax, bmax, fr, 
% then standard stuff, parameters at the end! 
% In plotbra, cmp=1 then yields out(1), etc. 
% Note that the full branch-output has the form 
% [bradat;out],  where bradat contains, e.g., counter, par-value, and 
% #neg.eigenvalues (used for bfbb.m). See bradat.m 
par=p.u(p.nu+1:end); u=u(1:p.nu); n=p.np; N=p.N; ng=n*N; M=p.mat.M(1:ng,1:ng);
vol=p.vol*N; chimin=par(19); chimax=par(20);
l1v=sum(M*u(1:ng))/vol; l2v=sqrt((u(1:ng)'*(M*u(1:ng)))/vol);
M=p.mat.M(1:n,1:n); vol=p.vol; chi=linspace(chimin,chimax,N); ba=0*chi; 
b=p.u(1:N*n); B=reshape(b,n,N);bmax0=max(abs(b)); 
for i=1:N; ba(i)=sum(M*B(:,i))/vol; end
 % chimax=chi(idx); 
% center of mass of biomass distribution
for i=1:N; baa(i)=chi(i)*sum(M*B(:,i))/vol; end
chicent = sum(baa)/sum(ba);
% biomass maxima
chii=linspace(0,1,700);
interbiochi=interp1(chi,ba,chii,'spline');
[biomax,idx]=max(interbiochi);
p.frcut=0.0001; 
try frcut=p.frcut; catch frcut=0.001; end 
[br,ir]=find(ba>frcut); brmin=chi(min(ir)); brmax=chi(max(ir)); 
out=[l1v; l2v; chicent;1;bmax0; par]; 
