function  [x1, x2] = insarpair(power, cohe, phase, npower)
%   [x1 x2] = insarpair(power, cohe, phase, npower);
%
%   Generates a  pair of SAR images according to the model (1,2)  
%   presented in IGARSS98
% 
%   MATLAB 5.x
%
%   x1 = z1*exp(-j*phi1)+n1
%   x2 = z3*exp(-j*phi2)+n2
%   power  - E[|z1|^2]=E[|z2|^2]
%   cohe   - E[z1*z2^{*}]/E[|z1|^2]
%   phase  - phase = phi2-phi1	(interferometric phase)
%   npower - additive noise power
%
%   Author J.M. Bioucas Dias 1998
%   Topic - SAR Interferometry

[M N] =size(power);
%-------------------------------------
% % generate w1 and w2 
% w1 = 1/sqrt(2)*(random('norm',0,1,M,N) + i*random('norm',0,1,M,N) ).* sqrt(power);
% w2 = 1/sqrt(2)*(random('norm',0,1,M,N) + i*random('norm',0,1,M,N) ).* sqrt(power);
% %-------------------------------------
% % generate n1 and n2 
% n1 = sqrt(npower/2)*(random('norm',0,1,M,N) + i*random('norm',0,1,M,N) );
% n2 = sqrt(npower/2)*(random('norm',0,1,M,N) + i*random('norm',0,1,M,N) );
% generate w1 and w2 
w1 = 1/sqrt(2)*(randn(M,N) + i*randn(M,N) ).* sqrt(power);
w2 = 1/sqrt(2)*(randn(M,N) + i*randn(M,N) ).* sqrt(power);
%-------------------------------------
% generate n1 and n2 
n1 = sqrt(npower/2)*(randn(M,N) + i*randn(M,N) );
n2 = sqrt(npower/2)*(randn(M,N) + i*randn(M,N) );

%-------------------------------------
% generate x1 and x2 
x1 = cohe.* w1+sqrt(1-cohe.^2).*w2+n1;
x2 = w1.* exp(-i*phase)+n2;








