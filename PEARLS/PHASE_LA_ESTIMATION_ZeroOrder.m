% PHASE_LA_ESTIMATION.m

% Point-wise adaptive window size estimation for masked and full neighbourhood
% Inputs:
% C - Estimate 
% C_initial - Initial Estimate
% Cosine and Sine observations data, z_x1=cos(z_est); z_x2=sin(z_est);
% Data size yN,xN
% Estimation point coordinates s_y,s_x,h
% HH a set of the window size parameters
% std of the noise - sigma
% mask matrix of the size yN,xN
% Outputs:
% Estimates with window sizes corresponding to the elements of HH- YY_h
% Adaptive window size estimates for the phase and the derivatives are the
% elements of the vector C_opt
% Indexes (numbers of HH) for adaptive window sizes h_opt_IN
%% V. Katkovnik, Tampere University of Technology, Tampere, Finland, October 2007


function [YY,std]=PHASE_LA_ESTIMATION_ZeroOrder(z_x1,z_x2,s_y,s_x,h,yN,xN,mask)


ddx=s_x+max(1-s_x,-h):s_x+min(h,xN-s_x); %% Grids of the function values taken boundary conditions
ddy=s_y+max(1-s_y,-h):s_y+min(h,yN-s_y);
  
    
 local_mask=mask(ddy,ddx);
cos1=z_x1(ddy,ddx).*local_mask;
sin1=z_x2(ddy,ddx).*local_mask;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HESSIAN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


C=atan2(sum(sin1(:)),sum(cos1(:))); %atan(sum(sin1(:))/sum(cos1(:)));



%C=atan2(sum(cos1(:)),sum(sin1(:)));
 
 std=1/sqrt(length(cos1(:)));

   YY=C;
