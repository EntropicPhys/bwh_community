function [kv, pv]=bfbb(kv,pv,dir,k0) % Brute Force Busse Ballon
% scans branch in dir for stable solns and returns the par-values 
% wave-nr k user supplied!  (could be obtained from soln via FFT, 
% but here we keep life simple) 
p=loadpp(dir); pv0=p.branch(4,:); inv=p.branch(3,:); % take data from branch 
lb=length(pv0); 
for i=1:lb;    % if stable, then add point (and wave-nr) to list 
    if inv(i)<1; kv=[kv k0]; pv=[pv pv0(i)]; end 
end 