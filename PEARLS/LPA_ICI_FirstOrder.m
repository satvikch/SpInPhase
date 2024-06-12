%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LPA FIRST ORDER USED FOR THE WINDOW SIZE SELECTION  %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Wrapped phase: z_est
%% window sizes: HH
%% mask matrix: mask
%% gamma-parameter of ICI rule: gamma

%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%
%% absolute phase estimate :phi_estimate
%% adaptive window sizes: h_opt
%% estimates : PHI_estimate_h

% V. Katkovnik, Tampere University of Technology, Tampere, Finland,
% October 2007
% Modified by J. Bioucas-Dias February 2008

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi_estimate, index, std_opt]=LPA_ICI_FirstOrder(z_est,HH,mask,gamma)
I=sqrt(-1);
maxwinsize = max(HH);

[xN,yN]=size(z_est);
lenh=length(HH);
% arrays
phi_estimate_h=zeros(xN,yN,lenh);

std_estimate_h=zeros(xN,yN, lenh);
maxmagnitu = zeros(xN,yN, lenh);

sigma11=function_stdEst2D_phase(z_est,1); % estimation of the standard
% deviation of the noise in the absolute phase
% Loops over the data grid:
for s_h=1:lenh
    h=HH(s_h);

    % Phase estimates
    
    [phi_estimate_h(:,:,s_h),c1_FFT,c2_FFT, maxmag] = denoiseLPAC(z_est,h*ones(xN,yN),16);
    [result,c1,c2,maxmagnitu(:,:,s_h)] = denoiGradC(z_est,h*ones(xN,yN),c1_FFT,c2_FFT,1e-4,16);
    phi_estimate_h(:,:,s_h) = result;
    std_estimate_h(:,:,s_h)=1/((2*h+1));
%     std_estimate_h(:,:,s_h)=1./sqrt(maxmagnitu(:,:,s_h));
end

[phi_estimate, index, std_opt]=function_ICI_ZeroOrderPhaseFiltering(phi_estimate_h,std_estimate_h*sigma11,gamma);

