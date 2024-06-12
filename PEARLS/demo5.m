%%  demo 5: reproduces PEARLS results shown in Figure 7,8,9 of paper
%
%  J. Bioucas-Dias, V. Katkovnik,  J. Astola, and  K.  Egiazarian, 
%  "Absolute phase estimation: adaptive local denoising and global unwrapping",  
%  vol. 47. no, 29, pp. 5358-5369, Applied Optics, 2008.
%  
% ---
%
%  The unwrapping step is computed with the PUMA algorithm introduced in
%
%  J. Bioucas-Dias  and  G. Valadão, "Phase unwrapping via graph cuts", 
%  IEEE Transactions on Image processing, vol. 16, no. 3, pp. 698-709, 2007.
%
%
%  ---------------------------------------------------------------------
%  data description
%  ---------------------------------------------------------------------
% 
% Original phase: Simulated of SAR Data based on a real surface. Data
% distributed with the book
%
% Two-Dimensional Phase  Unwrapping. Theory, algorithms, and Software, 
% Wiley Inter-Science, 1998.
%
%
%  Authors: J. Bioucas-Dias and V. Katkovnik. Januauy,  2007.
%
%%  begin demo
close all
clear all

tic
I = sqrt(-1);
% set of window sizes
HH=[1 2 3 4]; 
% recommended value for ICI window adaptive selection
gamma=2; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% LOAD ORIGINAL PHASE,  MASKS, AND QUOLITY MAPS      %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat('data\longs.mask.mat'))
load(strcat('data\longs.corr.mat'))
load(strcat('data\longs.phase.mat'))
load(strcat('data\longs.qual.mat'))
load(strcat('data\longs.surf.mat'))


% true surface
varphi = surf; 
% clear surf;
[M,N]=size(varphi);
% contour plot
[c,cl] = contour(varphi,[0  2  4 10:10:120]); 
clabel(c,cl);
colormap gray


% compute the interferogram (bring it to [-pi,pi});
% phi = angle(exp(I*phase));

R = ones(size(surf));
B = surf;
D = 0.8+corr*0.2;
% D = corr;
%%% Generate a noisy pair of SLC SAR images
% [z1, z2] = insarrnd(R, B, D);
[z1, z2] = insarpair(R, D, B, 0);
phi = angle(z1 .* conj(z2));
phi = phi.*mask;
varphi = B;
z = exp(I*phi);


figure(2);
imagesc(flipud(angle(exp(I*phi))))
title('Wrapped noisy phase','FontName', 'TimesNewRoman' , 'FontSize', 14);
colormap gray
axis off
drawnow;

figure(3);
imagesc(flipud(qual.*mask))
title('Quality Map','FontName','TimesNewRoman' , 'FontSize', 14);
colormap gray
axis off
drawnow;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PEARLS Algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[phi_hat, index_h] = PEARLS(phi, 'windows', HH, 'gamma', gamma, ...
                      'unwrapp','no');


figure(4);
imagesc(flipud(index_h)); colormap gray; colorbar; axis off; 
title('Window sizes at each pixel', 'FontSize', 14)
axis off




%fatten mask one pixel
mask = 1-fatten(1-mask);
mask(:,1) = 0;
mask(:,end) = 0;
mask(1,:) = 0;
mask(end,:) = 0;

% compute the ISNR inside the  mask and in the areas of qual == 1
ISNR = var((abs(exp(I*phi(:))-exp(I*varphi(:)))).*mask(:).*qual(:))/...
       var((abs(exp(I*phi_hat(:))-exp(I*varphi(:)))).*mask(:).*qual(:));

ISNR  = 10*log10(ISNR);



figure(5);
imagesc(phi_hat); colormap gray; axis off; 
title(strcat('Denoised wrapped phase; ISNR = ',num2str(ISNR,3),' dB'), ...
     'FontName', 'TimesNewRoman','FontSize', 16)


% isolate mask&qual pixels (those set zero) with discontinuities
[dh dv] = linefg(mask&qual);

%  set qualitymaps to PUMA
qualmap(:,:,1) = 0.99*dh;
qualmap(:,:,2) = 0.99*dv;

figure(6)
% we are  using quadractic potencial - discontinuities are infered beforehand)
potential.quantized = 'yes';
potential.threshold = 0;
phi_hat_puma =puma_ho(phi_hat,2,'potential',potential,'qualitymaps', qualmap);



tt= toc; %PEARLS total time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PEARLS ENDS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% plot error 
figure(7);
imagesc(flipud(((phi_hat_puma-varphi).*mask.*qual))); 
colormap gray; colorbar; axis off; 
title(strcat('(Unwrapped error) *  qualmap '), ...
     'FontName', 'TimesNewRoman','FontSize', 16)

% interpolate in areas where quality is zero
fprintf('Interpolate areas of  qual = 0\n')
phi_hat_recon = reconst(phi_hat_puma,qual,mask,0*dh,0*dv,300);


%total rmse (only in qual == 1 )
index = find(mask(:).*qual(:) == 1);
error = phi_hat_recon(index)-varphi(index);
RMSE = std(error(:));

% find outliers
index_outliers = find(abs(error(:)) > pi);
no_outliers = sum(abs(error(:)) > pi);
index_witout_out = index;
% remove outliers
index_witout_out(index_outliers) = [];

% rmse without ouliers
error = phi_hat_recon(index_witout_out)-varphi(index_witout_out);
RMSE_WOUTLIERS = std(error(:));



fprintf('\nISNR (wrapped) = %2.2f dB\n', ISNR)
fprintf('No. outliers = %d\n', no_outliers)
fprintf('RMSE (with outliers) =  %2.2f\n', RMSE)
fprintf('RMSE (without outliers) =  %2.2f\n', RMSE_WOUTLIERS)
fprintf('TIME = %2.2f sec\n\n', tt)



figure(8)
hist(error(:), 500)
title('Histogram of the phase surface error (without outliers) ', 'FontSize', 14)




figure(9)
surfl((phi_hat_recon).*mask);
shading interp;colormap(gray);
axis off;
title('PEARLS absolute phase estimate', 'FontSize', 14)




figure(10)
surfl((varphi).*mask);
shading interp;colormap(gray);
axis off;
title('True phase', 'FontSize', 14)


