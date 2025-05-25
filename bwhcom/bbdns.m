function [p,ds]=bbdns(p,ds,nst) % Busse balloon DNS with varying P 
P=p.u(p.nu+1); N=p.N;n=p.np;nx=n; % extracting data from p: prec.P and dims
x=getpte(p); x=x'; lx=max(x); xn=2*pi*x/lx; % space 
dt=0.0075; nt=floor(p.T/dt); pmod=nt/10; smod=0;  
for i=1:nst
 P=P+p.incr; p.u(p.nu+1)=P; t1=0; ts=[]; % change P, prepare DNS  
 pert=0*xn; pwn=[0.5 1:2:5]; npwn=length(pwn); pav=p.pa*[1:npwn]; 
 for l=1:npwn; pert=pert+pav(l)*cos(pwn(l)*xn); end % for pertubing H    
 p.u((N+1)*n+1:(N+2)*n)=p.u((N+1)*n+1:(N+2)*n).*(1+p.pa*pert);
 nc=0; [p,t1,ts,nc]=tintxs(p,t1,ts,dt,nt,nc,pmod,smod,@sGdns); % DNS
 out=sgbra(p,p.u); % output at end of DNS; now extract further data, i.e., 
 % wave number k of patterns, based on H, in different ways
 h=p.u(n*(N+1)+1:n*(N+2)); hdiff=max(h)-min(h); 
 if hdiff<1e-2; hdiff, break; end % stop if soln is bare soil or spatially hom 
 threshold=min(h)+p.threshp*max(h); % extract k by counting patches (via H)
 above_threshold=abs(h)>threshold;
 start_patches=find(diff([false; above_threshold(:)]) == 1);
 km2=2*pi*(length(start_patches)/lx);
 h=h-mean(h); hh=abs(fft(h));hh=movmean(hh,3);hh=hh(1:floor(nx/2)+1); % use FFT 
 [hmax,idx]=max(hh(1:end)); % find argmax(FFT(h))
 fs=nx/p.vol; fr=(0:floor(nx/2))*(fs/nx); fr1=fr(1:end); % FFT scaling 
 try; km=2*pi*fr1(idx); catch km=1; end  % extract k 
 bt=out(1); chi=out(3); P, ds=[ds [P;bt;chi;km;km2]]; % store data
 p.file.count=p.file.count+1; 
 if mod(i,10)==0; p.fuha.savefu(p); end  % save state (for plotting or reload)
end 
