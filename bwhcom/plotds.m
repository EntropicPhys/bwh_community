function h=plotds(ds,kcmp,kc,cc,i1,i2) 
% plotds - Plot output from bbdns.m showing normalized wavenumber trajectories 
% Inputs:
%   ds   - Data from bbdns: ds = [P; ||B||_1; chi; km; km2; km3]
%   kcmp - choose which 'k' to use: 4 = FFT, 5 = patch count, 6 = corrected count
%   kc   - Reference wavenumber for normalization
%   cc   - Column index for color (3 = chi, 1 = ||B||_1)
%   i1   - Number of initial points to skip
%   i2   - Number of final points to skip
len=size(ds,2); n1=1+i1; n2=len-i2; 
k_kc(1,:)=ds(kcmp,:)./kc;                       % Normalize wavenumber
plot(ds(1,n1:n2),k_kc(1,n1:n2),'-k');  hold on; % Plot wavenumber trajectory as line
scatter(ds(1,n1:n2),k_kc(1,n1:n2), 60, ds(cc,n1:n2), 'filled'); % Overlay with colored scatter plot based on chi
if cc==2                                        % Set colormap 
    colormap(flipud(summer)); 
else 
    colormap turbo; 
end 
h=colorbar; 
% Optional: uncomment below to fix color axis range
% caxis([0.2 0.95]);
