%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Denoise  a wraped interferogram by the gradient method given the fixed FFT window estimated parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%[result,c1,c2] = denoiGradC(phi,h_opt,c1_FFT,c2_FFT,delta,FFTWinsize);
%%                 INPUTS
% the noisy wrapped phase:                                                     phi
% adaptive window sizes of each pixel:                                         h_opt
% the first order polinomial parameters we get from the small windowed FFT:    c1_FFT,c2_FFT
% the Increment to compute the gradient:                                       delta
% the fixed FFT window size:                                                   FFTWinsize
%%                 OUTPUTS
% phase estimate:                                              result
% the first order polinomial parameter:                        c1,c2
% 
% Author: Hongxing Hao, July, 2012

