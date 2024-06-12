% [phi_hat,c1,c2] = denoiseLPAC(phi,h_opt,FFTWinsize)
% 
% C implementation of function denoiseLPA_Varying_Window_Sizes
%
%  Description:
%
%   Denoise  a wrapped interferogram by adjusting local first order polinomials  
%   (planes) by solving in tfe frequency domain the equations (4) and (5)
%   shown in 
%
%  J. Bioucas-Dias, V. Katkovnik,  J. Astola, and  K.  Egiazarian, 
%  "Absolute phase estimation: adaptive local denoising and global unwrapping",  
%  vol. 47. no, 29, pp. 5358-5369, Applied Optics, 2008.
%  
%
% 
%%
%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% phi           --    the noisy wrapped phase   
% h_opt         --    adaptive window sizes of each pixel 
% FFTWinsize    --    the fixed FFT window size. It determines the precsion
%                     of the solution. with is about  pi / FFTWinsize
%
%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% phi_hat    --   phase estimate       
% c2,c2      --  first order polinomial coefficients
% mag        --
%
% Author: Hongxing Hao, July, 2012

