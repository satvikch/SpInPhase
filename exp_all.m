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
%%%%%% GENERATE ORIGINAL PHASE,  NOISE, AND OBSERVED DATA %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exp 0 - Gaussian  surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M=100;
% N=100;
% varphi=gaussele(M,N,14*pi,10,15);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exp 1 - Gaussian  truncated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M=100;
N=100;
varphi=gaussele(M,N,14*pi,10,15);
varphi(1:M/2,1:N/2) = 0;
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exp 2 - Sinusoidal surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M=100;
% N=100;
% [X,Y] = meshgrid(-round(N/2): -round(N/2)+N-1, ...
%                     -round(M/2): -round(M/2)+M-1);
%    varphi = 6*sin(2*pi*(X/M+Y/M)*6);
%   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exp 3 - Sinusoidal  discontinuous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M=100;
% N=100;
% [X,Y] = meshgrid(-round(N/2): -round(N/2)+N-1, ...
%                     -round(M/2): -round(M/2)+M-1);
%        
%   mask = Y < 0;
%                 
%  varphi = 6*sin(2*pi*(X/M)*6).*mask -  6*sin(2*pi*(X/M)*6).*(1-mask);
%   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exp 4 - Mountains
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M = 100;
% N = 100;
% 
% S=2055615866; randn('seed', S);
% varphi=randn(M,N);
% fc = 0.01; % define cut frequency
% mask = zeros(M,N);
% mask(round(M/2-fc*M/2):round(M/2+fc*M/2), ...
%     round(N/2-fc*N/2):round(N/2+fc*N/2)) = 1;
% varphi =  2400* real(ifft2(fftshift(fftshift(fft2(varphi)).*mask)));
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exp 5 - Shear Planes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M=100;
% N=100;
% %build planes
% [X,Y] = meshgrid(0:N-1,0:M-1);
% %make discontinuous
% mask=ones(M,N); mask(:,1:N/2)=0;Y = Y-M/5;
% varphi=Y.*mask.*(Y>0);
% 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Exp 6 -  Long Peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load(strcat('data\longs.mask.mat'))
% load(strcat('data\longs.qual.mat'))
% load(strcat('data\longs.surf.mat'))
% 
% surf = surf(:,2:end);
% mask = mask(:,2:end);
% 
% %fatten mask one pixel
% mask = 1-fatten(1-mask);
% varphi = surf.*mask;
% 
% [M,N] = size(varphi);

    
%     [X,Y] = meshgrid(-round(N/2): -round(N/2)+N-1, ...
%                      -round(M/2): -round(M/2)+M-1);
%                  
%     varphi = -(X.^2+Y.^2)/M*pi; 

%         [X,Y] = meshgrid(-round(N/2): -round(N/2)+N-1, ...
%                      -round(M/2): -round(M/2)+M-1);
%                  
%          varphi = -(X.^2+Y.^2)/M*pi - (Y>0)*pi/2; 
%         
%          [X,Y] = meshgrid(-round(N/2): -round(N/2)+N-1, ...
%                      -round(M/2): -round(M/2)+M-1);
%                  
%          varphi = -(X.^2+Y.^2)/M*pi;
%          
%          varphi = varphi -  min(varphi(:));
%          
%          varphi = varphi .* ( Y < 0 | X > 0);
%          



% noise variance
sigma = 0.5;

% generate noise
% S=2055615866; randn('seed', S);
noise = sigma.*(randn(M,N)+I*randn(M,N))/sqrt(2);
% noise = (random('poiss',0,[M,N])+I*random('poiss',0,[M,N]));
% generate observed data
z = exp(I*varphi)+ noise;
% compute the interferogram;
phi = angle(z);

% allsigma = 0.3;
% sigma = allsigma/sqrt(2);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interferometric phase etimation with Complex dictionary learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  set parameter
param.K = 256;         % dict size
param.T = 500;         % iterations
param.sigma = sigma;   % noise value
param.patchnum = 64;   % batch size
param.lamda = 0.11;    %when sigma = 0.9  lamda = 0.15  when sigma < 0.9  lamda = 0.11
param.patsize = 10;    % patch size
param.IstrainDic = 1; 
param.IsUnwrap = 1;
% param.IsSubMean = 1;
% param.originalphase = varphi;

tic;

%[phi_pearls,c1,c2,mag] = denoiseLPAC(phi,1*ones(size(phi)),64);

%z_hat = z.*exp(-I*phi_pearls);

[phi_hat_SpInPHASE,outputcdl,alphacdl,unwrap_phase] = SpInPHASE(z, param);%0.5 noise
toc;

%phi_hat_SpInPHASE = angle(exp(I*(phi_pearls + phi_hat_SpInPHASE )));
t = toc;


% unwrap_phase = unwrap_phase-round((unwrap_phase(1)-varphi(1))/(2*pi))*2*pi;
ttemp = round(sum(sum(unwrap_phase-varphi))/(2*pi*M*N));
unwrap_phase = unwrap_phase-ttemp*2*pi;

%% Compare
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
title(strcat('Noisy wrapped phase,sigma = ', num2str(sigma,3)),...
    'FontName', 'TimesNewRoman' ,'FontSize', fontsize)
subplot(324);
surfl(unwrap_phase);shading interp;axis off;axis([1,M,1,N,min(varphi(:)),max(unwrap_phase(:))]);
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

%%                    with pre-learned dictionary
%load dic12X12.mat
%save  dic10X10.mat outputcdl
load  dic10X10.mat       %pretrained dictionary over 6 surfaces (Check TrainandSaveDict.m)

tic;
[phi_hat_SpInPHASE_pre,outputcdl_pre,alphacdl_pre,unwrap_phase_pre] = SpInPHASE(z, param, D);%0.5 noise
toc;
t_pre = toc;


% unwrap_phase_pre = unwrap_phase_pre-round((unwrap_phase_pre(1)-varphi(1))/(2*pi))*2*pi;
ttemp_pre = round(sum(sum(unwrap_phase_pre-varphi))/(2*pi*M*N));
unwrap_phase_pre = unwrap_phase_pre-ttemp_pre*2*pi;

%% Compare
wraperr_norm_pre = norm(wrap(phi_hat_SpInPHASE_pre - varphi),'fro')^2;
% RMSE_SpInPHASE_pre = sqrt(wraperr_norm_pre/M/N);
PSNR_SpInPHASE_pre = 10*log10(4*M*N*pi^2/wraperr_norm_pre);
% ISNR_SpInPHASE_pre = 10*log10(norm(wrap(phi - varphi),'fro')^2/wraperr_norm_pre);

unwrap_err_pre = unwrap_phase_pre-varphi;
errmask_pre = abs(unwrap_err_pre)<=pi;
% RMSE_a_SpInPHASE_pre = sqrt(norm(unwrap_err_pre.*errmask_pre,'fro')^2/sum(errmask_pre(:)));
% PSNR_a_SpInPHASE_pre = 10*log10(4*M*N*pi^2/norm(unwrap_err_pre,'fro')^2);
PPSNR_a_SpInPHASE_pre = 10*log10(4*sum(errmask_pre(:))*pi^2/norm(unwrap_err_pre.*errmask_pre,'fro')^2);
NELP_SpInPHASE_pre = M*N-sum(errmask_pre(:));



fprintf('\n\n------Output Information------\n');
fprintf('PSNR (SpInPHASE_pre) = %2.3f\n', PSNR_SpInPHASE_pre);
fprintf('PSNR_a (SpInPHASE_pre) = %2.3f\n', PPSNR_a_SpInPHASE_pre);
fprintf('NELP (SpInPHASE_pre) = %2.3f\n', NELP_SpInPHASE_pre);
fprintf('time of computation: %f\n',t_pre);

%%                    PEARLS
% set of window sizes
tic;
HH=[1 2 3 4]; 
gamma= 2*sigma; % ICI recommended value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PEARLS Algorithm
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[phi_hat, index_h] = PEARLS(phi, 'windows', HH, 'gamma', gamma, ...
                      'unwrapp','no');
                  
potential.quantized = 'no';
potential.threshold = 0.5;
figure;
phi_hat_puma = puma_ho(phi_hat,0.5, 'potential',potential); 

tt= toc;
phi_hat_WFT = phi_hat;
[N,M] = size(phi);

wraperr_norm = norm(wrap(phi_hat_WFT - varphi),'fro')^2;
PSNR_PEARLS = 10*log10(4*M*N*pi^2/wraperr_norm);

% unwrap_phase_1 = phi_hat_puma-round((phi_hat_puma(1)-varphi(1))/(2*pi))*2*pi;
ttemp_pe = round(sum(sum(phi_hat_puma-varphi))/(2*pi*M*N));
unwrap_phase_1 = phi_hat_puma-ttemp_pe*2*pi;

unwrap_err_pe = unwrap_phase_1-varphi;
errmask = abs(unwrap_err_pe)<=pi;
% PSNR_a_PEARLS = 10*log10(4*M*N*pi^2/norm(unwrap_err_pe,'fro')^2);
PPSNR_a_PEARLS = 10*log10(4*sum(errmask(:))*pi^2/norm(unwrap_err_pe.*errmask,'fro')^2);
NELP_PEARLS = M*N-sum(errmask(:));

fprintf('\n\n------Output Information------\n');
fprintf('PSNR (PEARLS) = %2.3f\n', PSNR_PEARLS);
fprintf('PSNR_a (PEARLS) = %2.3f\n', PPSNR_a_PEARLS);
fprintf('NELP (PEARLS) = %2.3f\n', NELP_PEARLS);
fprintf('time of computation: %f\n',tt);

%%                          WFT
tic;
g=wft2fw('wff',z,4,-pi,0.1,pi,4,-pi,0.1,pi,3*sigma);

phi_hat_WFT = angle(g.filtered);
%phi_hat_WFT = angle(g.r);


% wraperr_norm = norm(wrap(phi_hat_WFT - varphi),'fro')^2;
% % RMSE_WFT = sqrt(wraperr_norm/M/N);
% PSNR_WFT = 10*log10(4*M*N*pi^2/wraperr_norm);

potential.quantized = 'no';
potential.threshold = 0.5;
phi_hat_un = puma_ho(phi_hat_WFT,0.5, 'potential',potential); 

t_wft = toc;

ttemp_wft = round(sum(sum(phi_hat_un-varphi))/(2*pi*M*N));
unwrap_phase_wft = phi_hat_un-ttemp_wft*2*pi;


%% Compare
wraperr_norm_wft = norm(wrap(phi_hat_WFT - varphi),'fro')^2;
PSNR_WFT = 10*log10(4*M*N*pi^2/wraperr_norm_wft);

unwrap_err_wft = unwrap_phase_wft-varphi;
errmask_pre = abs(unwrap_err_wft)<=pi;
% PSNR_a_WFT = 10*log10(4*M*N*pi^2/norm(unwrap_err_wft,'fro')^2);
PPSNR_a_WFT = 10*log10(4*sum(errmask_pre(:))*pi^2/norm(unwrap_err_wft.*errmask_pre,'fro')^2);
NELP_WFT = M*N-sum(errmask_pre(:));



fprintf('\n\n------Output Information------\n');
fprintf('PSNR (WFT) = %2.3f\n', PSNR_WFT);
fprintf('PSNR_a (WFT) = %2.3f\n', PPSNR_a_WFT);
fprintf('NELP (WFT) = %2.3f\n', NELP_WFT);
fprintf('time of computation: %f\n',t_wft);


% [M,N] = size(phi_hat_un);
% ttemp = round(sum(sum(phi_hat_un-varphi))/(2*pi*M*N));
% unwrap_phase = phi_hat_un-ttemp*2*pi;
% 
% err_norm = norm(unwrap_phase - varphi,'fro')^2;
% % RMSE_SpInPHASE_err = sqrt(err_norm/M/N);
% PSNR_WFT_un = 10*log10(4*M*N*pi^2/err_norm);
% 
% unwrap_err = unwrap_phase-varphi;
% errmask = abs(unwrap_err)<=pi;
% NELP_WFT = M*N-sum(errmask(:));




% 
% fprintf('\n\n------Output Information------\n');
% fprintf('PSNR (WFT) = %2.3f\n', PSNR_WFT);
% fprintf('PSNR_u (WFT) = %2.3f\n', PSNR_WFT_un);
% fprintf('NELP (WFT) = %2.3f\n', NELP_WFT);
