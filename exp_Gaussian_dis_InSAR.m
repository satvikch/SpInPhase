%%  demo 2: SpInPHASE(Z is for complex data)
%
%
%  ---------------------------------------------------------------------
%  data description
%  ---------------------------------------------------------------------
%
% Original phase: Gaussian elevation with discontinuous with InSAR noise
%

%%  begin demo

clear all
close all



% stuff
I = sqrt(-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% GENERATE ORIGINAL PHASE,  NOISE, AND OBSERVED DATA %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=100;
N=100;
varphi=gaussele(M,N,14*pi,10,15);
varphi(1:M/2,1:N/2) = 0;

% noise variance
corr = 0.9; 
% generate noise
R = ones(size(varphi));
B = varphi;
D = corr*ones(size(varphi));
%%% Generate a noisy pair of SLC SAR images
[z1, z2] = insarrnd(R, B, D);
% [z1, z2] = insarpair(R, D, B, 0);
phi = angle(z1 .* conj(z2));
z = exp(I*phi);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interferometric phase etimation with Complex dictionary learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  set parameter

pp = insarrndVarToCo(500,10);
allsigma = ppval(pp,corr);
sigma = allsigma*sqrt(2);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interferometric phase etimation with Complex dictionary learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  set parameter
param.K = 256;
param.T = 500;
param.sigma = sigma;
param.patchnum = 64;
param.lamda = 0.15; %
param.OMPpara = 1.10;
param.patsize = 11;
param.IstrainDic = 1;
param.IsUnwrap = 1;
% param.originalphase = varphi;

tic;
[phi_hat_SpInPHASE,outputcdl,alphacdl,unwrap_phase] = SpInPHASE(z, param);%0.5 noise
toc;
t = toc;


ttemp = round(sum(sum(unwrap_phase-varphi))/(2*pi*M*N));
unwrap_phase = unwrap_phase-ttemp*2*pi;

wraperr_norm = norm(wrap(phi_hat_SpInPHASE - varphi),'fro')^2;
% RMSE_SpInPHASE = sqrt(wraperr_norm/M/N);
PSNR_SpInPHASE = 10*log10(4*M*N*pi^2/wraperr_norm);
% ISNR_SpInPHASE = 10*log10(norm(wrap(phi - varphi),'fro')^2/wraperr_norm);


unwrap_err = unwrap_phase-varphi;
errmask = abs(unwrap_err)<=pi;
% RMSE_a_SpInPHASE = sqrt(norm(unwrap_err.*errmask,'fro')^2/sum(errmask(:)));
% PSNR_a_SpInPHASE = 10*log10(4*M*N*pi^2/norm(unwrap_err,'fro')^2);
PPSNR_a_SpInPHASE = 10*log10(4*sum(errmask(:))*pi^2/norm(unwrap_err.*errmask,'fro')^2);
NELP_SpInPHASE = M*N-sum(errmask(:));

fprintf('\n\n------Output Information------\n');
fprintf('PSNR (SpInPHASE) = %2.3f\n', PSNR_SpInPHASE);
fprintf('PSNR_a (SpInPHASE) = %2.3f\n', PPSNR_a_SpInPHASE);
fprintf('NELP (SpInPHASE) = %2.3f\n', NELP_SpInPHASE);
fprintf('time of computation: %f\n',t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NL-InPhase
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters
hW = 10;                % search window of size (2hW+1)^2
hD = 5;                 % patches of size (2hD+1)^2
Lmin = 10;              % Minimum equivalent number of looks
nbit = 10;              % Number of iterations

tic;
    R_nlinsar = ones(size(z1));
    B_nlinsar = zeros(size(z1));
    D_nlinsar = zeros(size(z1));
    for i = 1:nbit
        fprintf(1, 'Iteration %d/%d\n', i, nbit);
        [R_nlinsar, B_nlinsar, D_nlinsar] = ...
            nlinsar(z1, z2, ...
                    R_nlinsar, B_nlinsar, D_nlinsar, ...
                    hW, hD, ...
                    Lmin);
    end
t1 = toc;
phi_hat = B_nlinsar;


phi_hat_unwrap = puma_ho(phi_hat ,0.5, 'verbose','no' );
ttemp_I = round(sum(sum(phi_hat_unwrap-varphi))/(2*pi*M*N));
phi_hat_unwrap = phi_hat_unwrap-ttemp_I*2*pi;


%% Compare
wraperr_norm = norm(wrap(phi_hat - varphi),'fro')^2;
RMSE_NLInSAR = sqrt(wraperr_norm/M/N);
PSNR_NLInSAR = 10*log10(4*M*N*pi^2/wraperr_norm);

unwrap_err = phi_hat_unwrap-varphi;
errmask = abs(unwrap_err)<=pi;
PPSNR_a_NLInSAR = 10*log10(4*sum(errmask(:))*pi^2/norm(unwrap_err.*errmask,'fro')^2);
NELP_NLInSAR = M*N-sum(errmask(:));


fprintf('\n\n------Output Information------\n');
fprintf('PSNR (NL-InSAR) = %2.3f\n', PSNR_NLInSAR);
fprintf('PSNR_a (NL-InSAR) = %2.3f\n', PPSNR_a_NLInSAR);
fprintf('NELP (NL-InSAR) = %2.3f\n', NELP_NLInSAR);
fprintf('time of computation: %f\n',t1);








figure(1);colormap gray;
set (gcf,'position',[232 246 400 600])
fontsize = 10;
subplot(321);
imagesc(wrap(varphi)); axis off;
title('original wrapped phase', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(322);
surfl(varphi);shading interp;axis off; axis([1,M,1,N,1,max(varphi(:))]);
title('Original surface', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(323);
imagesc(phi);  axis off;
title(strcat('Noisy wrapped phase,sigma = ', num2str(sigma,3)),...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(324);
surfl(unwrap_phase);shading interp;axis off;axis([1,M,1,N,1,max(unwrap_phase(:))]);
title('Estimate surface', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(325);
imagesc(phi_hat_SpInPHASE);  axis off;
title('Denoised wrapped phase', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(326);
imagesc(wrap(phi_hat_SpInPHASE - varphi));  axis off; 
title('Wrapped phase estimate error', ...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
figure(2); axis off;
displayDic(outputcdl);

figure(3); axis off;
hist(sum(alphacdl~=0),1:max(sum(alphacdl~=0)));