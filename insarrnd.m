function [z1 z2] = insarrnd(R, B, D)

% INSARRND
%   [Z1 Z2] = INSARRND(R, B, D) generates InSAR noisy data
%
%       ARGUMENT DESCRIPTION:
%               R - Reflectivity image (greater than zero)
%               B - Phase image
%               D - Coherence image (between 0 and 1)
%
%       OUTPUT DESCRIPTION:
%               [z1 z2] - noisy pair of SLC image s.t z1.*conj(z2)
%                         provides a noisy interferogram
%
%       AUTHOR:
%               Charles Deledalle
%               email: deledalle@telecom-paristech.fr
%
%       VERSION: 23 August 2010
    
    l1 = sqrt(R);
    l2 = sqrt(R) .* D .* exp(-1i * B);
    l3 = sqrt(R) .* sqrt(1 - D.^2);
    
%     S=2055615866; randn('seed', S);
    
    x1 = (randn(size(D)) + 1i * randn(size(D))) / sqrt(2);
    x2 = (randn(size(D)) + 1i * randn(size(D))) / sqrt(2);
    z1 = l1 .* x1;
    z2 = l2 .* x1 + l3 .* x2;
