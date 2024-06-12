%%  demo 1: reproduces PEARLS results shown in Figure 2 of paper
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
% Original phase: Gaussian elevation 
%  
% observation model  z = exp(j phi) + n
%
% SNR = 3dB  (sigma = 0.5)
%
%
%  Authors: J. Bioucas-Dias and V. Katkovnik. Januauy,  2007.
%
%%  begin demo

clear all
close all

tic
% stuff
I = sqrt(-1);

% set of window sizes
HH=[1 2 3 4]; 
gamma=2.; % ICI recommended value


% noise variance
sigma = 0.5;
SNR = 10*log10(1/2/sigma^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GENERATE ORIGINAL PHASE,  NOISE, AND OBSERVED DATA %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=150;
N=100;

%build planes
[X,Y] = meshgrid(0:N-1,0:M-1);
%make discontinuous
mask=ones(M,N); mask(:,1:N/2)=0;
varphi=Y.*mask;

clear X Y;

% generate noise
S=2055615866; randn('seed', S);
noise = sigma*(randn(M,N)+I*randn(M,N));
% generate observed data
z = exp(I*varphi)+ noise;
% compute the interferogram;
phi = angle(z);

figure(1);
imagesc(phi); colormap gray; axis off; 
title(strcat('Wrapped phase; SNR = ',num2str(SNR,3),' dB'), ...
     'FontName', 'TimesNewRoman','FontSize', 16)
 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PEARLS Algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[phi_hat, index_h] = PEARLS(phi, 'windows', HH, 'gamma', gamma, ...
                      'unwrapp','no');


figure(2);
imagesc(index_h); colormap gray; colorbar; axis off; 
title('Window size at each pixel', 'FontSize', 16)


% compute ISNR with respect to the wrapped phase, as defined in  the paper
ISNR = var(abs(exp(I*phi(:))-exp(I*varphi(:))))/...
       var(abs(exp(I*phi_hat(:))-exp(I*varphi(:))));
 
ISNR  = 10*log10(ISNR);


figure(3);
imagesc(phi_hat); colormap gray; axis off; 
title(strcat('Denoised wrapped phase; ISNR = ',num2str(ISNR,3),' dB'), ...
     'FontName', 'TimesNewRoman','FontSize', 16)

figure(4)% to show PUMA evolution
% unwrapp  the denoise interferogram (absolute phase reconstruction)
potential.quantized = 'no';
potential.threshold = pi;
% PUMA unwrapping
phi_hat_puma=puma_ho(phi_hat,0.5,'schedule',[1]);


title(strcat('Denoised wrapped phase; SNR = ',num2str(ISNR,3),' dB'), ...
     'FontName', 'TimesNewRoman','FontSize', 16)


tt= toc; %PEARLS total time

figure(5), surfl(phi_hat_puma);shading interp;colormap(gray);
axis off;
title('PEARLS estimate', 'FontSize', 16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PEARLS ENDS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(6)
% unwrapp the noisy interferogram to compair with the PEARLS result
phi_puma = puma_ho(phi,0.5, 'potential',potential, 'schedule', [1]);


figure(6), surfl(phi_puma);shading interp;colormap(gray);
axis off;
title('', 'FontSize', 16)
title('Noisy interferogram unwrapped', 'FontSize', 16)


%compute performance indicatots w.r.t. th unwrapped phase
fprintf('\n\nPerformance indicators w.r.t. the unwrapped phase\n ')
sds_ICI_unwrap1 = function_Errors_NoNormalization(varphi,phi_hat_puma, ...
                                                   phi_puma,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%%%%%%%%%%%%%% SHOW  IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(7); colormap gray; 
subplot(221);
imagesc(phi); colormap gray; axis off;
 title(strcat('Wrapped phase; SNR = ',num2str(SNR,3),' dB'), ...
     'FontName', 'TimesNewRoman' , 'FontSize', 10)

subplot(222);
imagesc(index_h); colormap gray; colorbar; axis off; 
title('Window sizes at each pixel',  ...
       'FontName', 'TimesNewRoman' ,'FontSize', 10)

fprintf('\n\nPerformance indicators w.r.t. the wrapped phase\n ')
fprintf('\nISNR = %2.2f dB\n', ISNR)

fprintf('\nTIME = %2.2f sec\n', tt)



subplot(223);
imagesc(wrap(phi_hat_puma)); colormap gray; axis off;
title(strcat('Denoised Phase; ISNR = ',num2str(ISNR,3),' dB'), ...
    'FontName', 'TimesNewRoman' , 'FontSize', 10)

subplot(224);
surfl(phi_hat_puma);shading interp;colormap(gray);
axis off;
title('Absolute phase estimate', ...
     'FontName', 'TimesNewRoman' ,'FontSize', 10)


