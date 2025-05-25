function out=sgbra(p,u)
par=p.u(p.nu+1:end); u=u(1:p.nu); n=p.np; M=p.mat.M(1:n,1:n); 
l1v=sum(M*u(1:n))/p.vol; l2v=sqrt((u(1:n)'*(M*u(1:n)))/p.vol); 
out=[l1v; l2v; max(abs(u(1:n))); min(abs(u(1:n))); par]; 