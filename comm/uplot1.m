function uplot1(p,wnr) % special windows plot for continuation
N=p.N;n=p.np; % structure dim. 
 mclf(wnr); % clean and assign window num 'wnr'
p1=[1:4 7:10 13:16 19:22 25:28]; %  windows distribution
p2=37:40; p3=49:52; p4=6:6:28; 
par=p.u(p.nu+1:end);chimin=par(19);chimax=par(20); % problem dims. 
y=linspace(chimin,chimax,N);x=getpte(p);x1=min(x);x2=max(x);
[X,Y]=meshgrid(y,x); % 2D mesh for plot
b=p.u(1:N*n);B=reshape(b,n,N);bmax=3; % spatial-trait biomass re-distr.
w=p.u(N*n+1:(N+1)*n); h=p.u((N+1)*n+1:(N+2)*n); % water variables
subplot(9,6,p1);surf(X,Y,B);caxis([0 bmax]); % spatial-trait biomass plot
view(90,270);shading interp;axis tight;
colormap(flipud(summer));tits=[p.file.pname mat2str(p.file.count-1)]; 
colorbar('southoutside');title(tits);set(gca,'fontsize',12); 
subplot(9,6,p2);plot(x,w,'linewidth',2); % soil-water cont. plot
axis tight;title('w');set(gca,'fontsize',12); 
try; z1=round(min(w)*105)/100; z2=round(max(w)*95)/100; yticks([z1 z2]); 
catch; try; z1=0.9*min(w); z2=1.1*max(w); axis([x1 x2 z1 z2]); zt=round(50*(z1+z2))/100; 
    yticks(zt); catch; end 
end 
subplot(9,6,p3); plot(x,h,'linewidth',2); % surface-water plot 
axis tight; title('h');set(gca,'fontsize',12);  
try; z1=round(min(h)*105)/100; z2=round(max(h)*95)/100; yticks([z1 z2]); 
catch; try; z1=0.9*min(h); z2=1.1*max(h); axis([x1 x2 z1 z2]); zt=round(50*(z1+z2))/100; 
    yticks(zt); catch; end 
end 
M=p.mat.M(1:n,1:n); vol=p.vol; ba=0*y; % spatial biomass integ.
for i=1:N
    ba(i)=sum(M*B(:,i))/vol;
end
subplot(9,6,p4); plot(ba,y,'linewidth',2); % biomass cont. over trait space
xlim([0 bmax])set(gca,'fontsize',12);  
ylabel('\chi'); xlabel('<B(.,\chi)>'); grid on; 