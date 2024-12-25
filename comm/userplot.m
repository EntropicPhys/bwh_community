function userplot(p,wnr) % userplot, called for p.plot.pstyle=-1; 
% here just interface to uplot1,...,  depending on GLOBAL p2pglob.ps 
global p2pglob; 
switch p2pglob.ps 
    case 1; uplot1(p,wnr); 
    case 2; uplot2(p,wnr); 
    case 3; uplot3(p,wnr);       
end 