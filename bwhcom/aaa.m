cmds1: cont of HOM, swibra to Turing branches and first snake, 
       incl.some DNS (to get to HOM, and to get to other branches by 
       "dropping patches")
cmds2: Busse balloon (BB) computations: scan precomputed branches (from cmds1) 
       for stable states (bfbb),  and step through BB by DNS with varying P (bbdns) 
--------
bwhinit:     initialization 
oosetfemops: assemble and store system matrices 
nodalf:    "nonlinearity", i.e., terms without spatial derivatives 
sG,sGjac:  rhs and Jacobian 
sGdns:  for DNS of community model; essentially rhs without x-diffusion    
sgbra:     mod of library function stanbra 
tintxs:    semi-implicit time stepper for DNS
bfbb:      brute force BB: scans branch for stable solns and returns the resp P
bbdns:     Busse balloon DNS with varying P  
---------
userplot:    custom plotting; here just an interface to uplot1,2,3 (see there) 
plotds  :    plot data from bbdns (from stepping through BB) 
