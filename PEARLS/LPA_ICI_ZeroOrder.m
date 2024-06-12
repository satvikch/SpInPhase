%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% LPA ZERO ORDER USED FOR THE WINDOW SIZE SELECTION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% true signal: y_true
%% true signals derivatives on x and y : y_x,y_y 
%% cos and sin of the observed wrapped phase: z_x1,z_x2
%% Wrapped phase: z_est
%% derivatives of the wrapped phase: z_x_est,z_y_est
%% window sizes: HH
%% indicator of the used initialization: true_initialization
%% mask matrix: mask
%% gamma-parameter of ICI rule: gamma

%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%
%% absolute phase estimate :phi_estimate
%% absolute phase estiamte derivatives: phi_x_estimate,phi_y_estimate
%% adaptive window sizes: h_opt
%% estimates (approximate) with constant window sizes: PHI_estimate_h

% V. Katkovnik, Tampere University of Technology, Tampere, Finland,
% October 2007
%
% Modified by J Bioucas, September, 2012 
%             denoiseLPACZero -> implemented in C
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phi_estimate, index, std_opt]=LPA_ICI_ZeroOrder(z_est,HH,mask,gamma);

%ddelta=1;
% disp('------------------------')
% disp('LPA ZERO ORDER filtering')
% disp('------------------------')

[yN,xN]=size(z_est);
lenh=length(HH);
% arrays
phi_estimate_h=zeros(yN,xN,lenh);
std_estimate_h=zeros(yN,xN, lenh);
% estimation of the standard 
sigma11=function_stdEst2D_phase(z_est,1); 


%%z_x1=cos(z_est); z_x2=sin(z_est);
% deviation of the noise in the absolute phase

% % 
% % %Loops over the data grid: 
% % for s_h=1:lenh
% %     h=HH(s_h);
% % 
% %  for s_y=1:yN
% % %   %  s_y;
% % % 
% % % %for s_x=h+1:xN-h-1
% %  for s_x=1:xN
% %      %%%% PLACE FOR FUSING WITH Zstep algorithm %%%%%%%%%%%%%%%%%%%%%%%%
% %     
% % % Zero Oeder LPA
% % [phi_estimate1, std_C1]=PHASE_LA_ESTIMATION_ZeroOrder(z_x1,z_x2,s_y,s_x,h,yN,xN,mask);
% % 
% % phi_estimate_h(s_y,s_x,s_h)=phi_estimate1;  % Phase estimates 
% % std_estimate_h(s_y,s_x,s_h)=std_C1;
% % 
% % end
% % end
% % end
% % 

% C implementation of zero order  LPA 
for i=1:lenh
    [phi_estimate_h(:,:,i),std_estimate_h(:,:,i)] = denoiseLPACZero(z_est,HH(i)*ones(yN,xN));
    std_estimate_h(:,:,i) = 1./sqrt(std_estimate_h(:,:,i));
end


[phi_estimate, index, std_opt]=function_ICI_ZeroOrderPhaseFiltering(phi_estimate_h,std_estimate_h*sigma11,gamma);




% keyboard