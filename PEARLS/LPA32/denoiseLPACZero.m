% [phi_hat,wout_size] =  denoiseLPACZero(phi,wsize);
% 
% C implementation of zero oeder LPA
%
%  Description:
%
%   Denoise  a wrapped interferogram by adjusting local zero order polinomials  
%   (constants). See
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
% wsize         --    window sizes matrix,  with the same size as phi. 
%                     The element wsize(i,j) contains the size h of the 
%                     square window centered at (i,j), withsize (1+2h)^2
%
%
%
%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% phi_hat    --  zero order LPA  phase estimate       
% wout_size  --  matrix containing the number of pixel of   intersection 
%               of the window with the image grid. 
%
% Author: Hongxing Hao, July, 2012

