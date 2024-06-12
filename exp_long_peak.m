%%  demo 1: SpInPHASE(Z is for complex data)
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
%%  begin demo

clear all
close all



% stuff
I = sqrt(-1);

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exp -  Long Peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('longs.152x458.mask','r','b');
mask = fread(fid,'uint8');
mask = reshape(mask, [152,458]);
mask = mask/255;

%fatten mask one pixel
mask = 1-fatten(1-mask);
mask(:,1) = 0;
mask(:,end) = 0;
mask(1,:) = 0;
mask(end,:) = 0;


fid=fopen('longs.152x458.corr','r','b');
corr = fread(fid,'uint8');
corr = reshape(corr, [152,458]);
corr = corr/255;
[NA, corr] = BM3D(1, corr, 0.15*255);
% corr = corr.*mask;

fid=fopen('longs.152x458.phase','r','b');
phase = fread(fid,'uint8');
phase = reshape(phase, [152,458]);
phase = wrap(phase/255*(2*pi));

fid=fopen('longs.152x458.surf','r','b');
surf = fread(fid,'float');
surf = reshape(surf, [152,458]);

load longs.qual.mat
% qual = corr>0.45;


% varphi = surf.*mask;% 
% [M,N] = size(varphi);% 
% z = exp(1i*phase);
% % z = z.*mask;
% % compute the interferogram;
% phi = angle(z); 
% % allsigma = 0.3;
% % sigma = allsigma/sqrt(2);

varphi = surf;
[M,N]=size(varphi);
% generate noise
% R = ones(size(surf));
R = corr+1;
B = surf;
D = corr;
%%% Generate a noisy pair of SLC SAR images
[z1, z2] = insarrnd(R, B, D);
% [z1, z2] = insarpair(R, D, B, 0);
phi = angle(z1 .* conj(z2));
z = exp(I*phi.*mask);
% z = z1 .* conj(z2);


pp = insarrndVarToCo(500,10);
sigma = ppval(pp,corr);
sigma = max(sigma*sqrt(2),0.25);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interferometric phase etimation with Complex dictionary learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  set parameter
param.K = 512;
param.T = 500;
param.sigma = sigma;
param.patchnum = round(0.0064*sum(mask(:)));
param.lamda = 0.35;
param.patsize = 10;
param.IstrainDic = 1;
param.IsUnwrap = 0;
param.mask = mask.*qual;
% param.IsSubMean = 1;
% param.originalphase = varphi;

tic;

[phi_hat_SpInPHASE,outputcdl,alphacdl] = SpInPHASE(z, param);%0.5 noise


[dh dv] = linefg(mask&qual);
%  set qualitymaps to PUMA
qualmap(:,:,1) = 0.99*dh;
qualmap(:,:,2) = 0.99*dv;
% we are  using quadractic potencial - discontinuities are infered beforehand)
potential.quantized = 'no';
potential.threshold = 2;
unwrap_phase =puma_ho(phi_hat_SpInPHASE,2,'potential',potential,'qualitymaps', qualmap,'verbose','no');

toc;
%phi_hat_SpInPHASE = angle(exp(I*(phi_pearls + phi_hat_SpInPHASE )));
t = toc;


unwrap_phase = unwrap_phase.*mask;

numberpix = sum(sum(mask.*qual));
ttemp = round(sum(sum((unwrap_phase-varphi).*mask.*qual))/(2*pi*numberpix));
unwrap_phase = unwrap_phase-ttemp*2*pi;

wraperr_norm = norm(wrap(phi_hat_SpInPHASE - varphi).*mask.*qual,'fro')^2;
% RMSE_WFT = sqrt(wraperr_norm/numberpix);
PSNR_SpInPHASE = 10*log10(4*numberpix*pi^2/wraperr_norm);

% err_norm = norm((unwrap_phase - varphi).*mask.*qual,'fro')^2;
% RMSE_SpInPHASE_err = sqrt(err_norm/M/N);
% PSNRa_SpInPHASE = 10*log10(4*numberpix*pi^2/err_norm);

unwrap_err = (unwrap_phase-varphi).*mask.*qual;
errmask = abs(unwrap_err)<=pi;
errmask = errmask.*mask.*qual;
PPSNR_a_SpInPHASE = 10*log10(4*sum(errmask(:))*pi^2/norm(unwrap_err.*errmask,'fro')^2);
NELP_SpInPHASE = numberpix-sum(errmask(:));


fprintf('\n\n------Output Information------\n');
fprintf('PSNR (SpInPHASE) = %2.3f\n', PSNR_SpInPHASE);
fprintf('PSNR_u (SpInPHASE) = %2.3f\n', PPSNR_a_SpInPHASE);
fprintf('NELP (SpInPHASE) = %2.3f\n', NELP_SpInPHASE);
fprintf('time of computation SpInPHASE: %f\n',t);


figure(1);colormap gray;
set (gcf,'position',[232 246 400 600])
fontsize = 10;
subplot(321);
imagesc(wrap(varphi)); axis off;
title('original wrapped phase', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(322);
surfl(varphi);shading interp;axis off; axis([1,M,1,N,min(varphi(:)),max(varphi(:))]);
title('Original surface', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(323);
imagesc(phi);  axis off;
title('Noisy wrapped phase',...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(324);
surfl(unwrap_phase);shading interp;axis off;axis([1,M,1,N,min(varphi(:)),max(unwrap_phase(:))]);
title('Estimate surface', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(325);
imagesc(phi_hat_SpInPHASE.*mask);  axis off;
title('Denoised wrapped phase', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(326);
imagesc(wrap(phi_hat_SpInPHASE - varphi).*mask.*qual);  axis off; 
title('Wrapped phase estimate error', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)

figure(2); axis off;
displayDic(outputcdl);

figure(3); axis off;
hist(sum(alphacdl~=0),1:max(sum(alphacdl~=0)));

%%                    with pre-learned dictionary
%load dic12X12.mat
%save  dic10X10.mat outputcdl
load  dic10X10.mat 

tic;
[phi_hat_SpInPHASE_pre,outputcdl_pre,alphacdl_pre] = SpInPHASE(z, param, D);%0.5 noise


% we are  using quadractic potencial - discontinuities are infered beforehand)
potential.quantized = 'no';
potential.threshold = 2;
unwrap_phase_pre =puma_ho(phi_hat_SpInPHASE_pre,2,'potential',potential,'qualitymaps', qualmap,'verbose','no');

toc;
t_pre = toc;


unwrap_phase_pre = unwrap_phase_pre.*mask;

numberpix = sum(sum(mask.*qual));
ttemp = round(sum(sum((unwrap_phase_pre-varphi).*mask.*qual))/(2*pi*numberpix));
unwrap_phase_pre = unwrap_phase_pre-ttemp*2*pi;

wraperr_norm = norm(wrap(phi_hat_SpInPHASE_pre - varphi).*mask.*qual,'fro')^2;
% RMSE_WFT = sqrt(wraperr_norm/numberpix);
PSNR_SpInPHASE_pre = 10*log10(4*numberpix*pi^2/wraperr_norm);

% err_norm = norm((unwrap_phase_pre - varphi).*mask.*qual,'fro')^2;
% RMSE_SpInPHASE_err = sqrt(err_norm/M/N);
% PSNRa_SpInPHASE_pre = 10*log10(4*numberpix*pi^2/err_norm);


unwrap_err_pre = (unwrap_phase_pre-varphi).*mask.*qual;
errmask_pre = abs(unwrap_err_pre)<=pi;
errmask_pre = errmask_pre.*mask.*qual;
PPSNR_a_SpInPHASE_pre = 10*log10(4*sum(errmask_pre(:))*pi^2/norm(unwrap_err_pre.*errmask_pre,'fro')^2);
NELP_SpInPHASE_pre = numberpix-sum(errmask_pre(:));


fprintf('\n\n------Output Information------\n');
fprintf('PSNR (SpInPHASE_pre) = %2.3f\n', PSNR_SpInPHASE_pre);
fprintf('PSNR_u (SpInPHASE_pre) = %2.3f\n', PPSNR_a_SpInPHASE_pre);
fprintf('NELP (SpInPHASE_pre) = %2.3f\n', NELP_SpInPHASE_pre);
fprintf('time of computation SpInPHASE_pre: %f\n',t_pre);

%%                          WFT
tic;
g=wft2fw('wff',z,4,-pi,0.1,pi,4,-pi,0.1,pi,3*mean(sigma(:)));

phi_hat_WFT = angle(g.filtered);
%phi_hat_WFT = angle(g.r);


% wraperr_norm = norm(wrap(phi_hat_WFT - varphi),'fro')^2;
% % RMSE_WFT = sqrt(wraperr_norm/M/N);
% PSNR_WFT = 10*log10(4*M*N*pi^2/wraperr_norm);

% potential.quantized = 'no';
% potential.threshold = 0.5;
% phi_hat_un = puma_ho(phi_hat_WFT,0.5, 'potential',potential); 


% we are  using quadractic potencial - discontinuities are infered beforehand)
potential.quantized = 'no';
potential.threshold = 2;
phi_hat_un =puma_ho(phi_hat_WFT,2,'potential',potential,'qualitymaps', qualmap,'verbose','no');

t_wft = toc;

phi_hat_un = phi_hat_un.*mask;

numberpix = sum(sum(mask.*qual));
ttemp = round(sum(sum((phi_hat_un-varphi).*mask.*qual))/(2*pi*numberpix));
phi_hat_un = phi_hat_un-ttemp*2*pi;

wraperr_norm = norm(wrap(phi_hat_WFT - varphi).*mask.*qual,'fro')^2;
% RMSE_WFT = sqrt(wraperr_norm/numberpix);
PSNR_WFT = 10*log10(4*numberpix*pi^2/wraperr_norm);

err_norm = norm((phi_hat_un - varphi).*mask.*qual,'fro')^2;
% RMSE_SpInPHASE_err = sqrt(err_norm/M/N);
% PSNR_WFT_un = 10*log10(4*numberpix*pi^2/err_norm);

unwrap_err = (phi_hat_un-varphi).*mask.*qual;
errmask = abs(unwrap_err)<=pi;
errmask = errmask.*mask.*qual;
PPSNR_a_WFT = 10*log10(4*sum(errmask(:))*pi^2/norm(unwrap_err.*errmask,'fro')^2);
NELP_WFT = numberpix-sum(errmask(:));


fprintf('\n\n------Output Information------\n');
fprintf('PSNR (WFF) = %2.3f\n', PSNR_WFT);
fprintf('PSNR_u (WFF) = %2.3f\n', PPSNR_a_WFT);
fprintf('NELP (WFF) = %2.3f\n', NELP_WFT);
fprintf('time of computation WFF: %f\n',t_wft);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NL-InPhase
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters
hW = 10;                % search window of size (2hW+1)^2
hD = 3;                 % patches of size (2hD+1)^2
Lmin = 10;              % Minimum equivalent number of looks
nbit = 10;              % Number of iterations

tic;
R_nlinsar = corr+1;
    B_nlinsar =phi;
    D_nlinsar = corr;
%     R_nlinsar = ones(size(z1));
%     B_nlinsar = zeros(size(z1));
%     D_nlinsar = zeros(size(z1));
    for i = 1:nbit
        fprintf(1, 'Iteration %d/%d\n', i, nbit);
        [R_nlinsar, B_nlinsar, D_nlinsar] = ...
            nlinsar(z1, z2, ...
                    R_nlinsar, B_nlinsar, D_nlinsar, ...
                    hW, hD, ...
                    Lmin);
    end
t_NLInSAR = toc;

phi_hat_NLINSAR = B_nlinsar;
phi_hat_NLINSAR = phi_hat_NLINSAR.*mask;

potential.quantized = 'no';
potential.threshold = 2;
phi_hat_NLun =puma_ho(phi_hat_NLINSAR,2,'potential',potential,'qualitymaps', qualmap,'verbose','no');

numberpix = sum(sum(mask.*qual));
ttemp = round(sum(sum((phi_hat_NLun-varphi).*mask.*qual))/(2*pi*numberpix));
phi_hat_NLun = phi_hat_NLun-ttemp*2*pi;

wraperr_norm = norm(wrap(phi_hat_NLINSAR - varphi).*mask.*qual,'fro')^2;
% RMSE_WFT = sqrt(wraperr_norm/numberpix);
PSNR_NLInSAR = 10*log10(4*numberpix*pi^2/wraperr_norm);

% err_norm = norm((phi_hat_NLun - varphi).*mask.*qual,'fro')^2;
% RMSE_SpInPHASE_err = sqrt(err_norm/M/N);
% PSNR_WFT_un = 10*log10(4*numberpix*pi^2/err_norm);

unwrap_err = (phi_hat_NLun-varphi).*mask.*qual;
errmask = abs(unwrap_err)<=pi;
errmask = errmask.*mask.*qual;
PPSNR_a_NLInSAR = 10*log10(4*sum(errmask(:))*pi^2/norm(unwrap_err.*errmask,'fro')^2);
NELP_NLInSAR = numberpix-sum(errmask(:));


fprintf('\n\n------Output Information------\n');
fprintf('PSNR (NL-InSAR) = %2.3f\n', PSNR_NLInSAR);
fprintf('PSNR_u (NL-InSAR) = %2.3f\n', PPSNR_a_NLInSAR);
fprintf('NELP (NL-InSAR) = %2.3f\n', NELP_NLInSAR);
fprintf('time of computation NL-InSAR: %f\n',t_NLInSAR);