% ICI rule for adaptive image estimation
% ----------------------------------------------------------------------
%
% Complete version of the ICI algorithm, including weighted
% order-statistics filtering (WOS) and adaptive variance (AV).
%
%
% [YY, IN,std_opt]=function_ICI(B,std_B,gam,N,O);
%
% INPUTS:
% B      : an array of image estimates obtained for set of increasing window sizes, size(B)=N1*N2*N_h
% std_B  : an array with the corresponding estimates' standard deviations
% gam    : the ICI Gamma threshold (a scalar)
%
% OPTIONAL INPUTS:
% N      : the cell array providing the weight matrices for the weighted
%          order-statistics (WOS) filtering for each ICI iteration.
% O      : vector of the weighted statistics offsets.
%           ( offset=0 corresponds to standard weighted median
%             e.g. N={ones(3),ones(3),ones(1)} and O=(0 0 0)
%                  will result in 2 recursive 3x3 median filters )
%        weights have to be odd1 x odd2 matrices and sum of weights + offset
%        should be odd
%
% OUTPUTS:
% YY      : is a N1xN2 matrix of the adaptive estimates
% IN      : is a N1xN2 matrix of the adaptive estimate/scale indexes (1 corresponds to the smallest scale)
% std_opt : is a N1xN2 matrix of of the adaptive estimates' standard deviations
%
%
% Alessandro Foi - Tampere University of Technology - 2004-2005
% V. Katkovnik, Tampere University of Technology, Tampere, Finland, October
% 2007, modified 

% ---------------------------------------------------------------
function [YY, IN,std_opt]=function_ICI_ZeroOrderPhaseFiltering(B,std_B,gam,N,O);
[N1,N2,N_h]=size(B);

% Upper and Lower bounds
D=gam*std_B;
U=B+D;
L=B-D;

IN=ones(N1,N2);   % indexes
YY=B(:,:,1);
TT=cell(N_h,1);   % truth table
TT{1}=~zeros(N1,N2);

% ICI rule into practice...
for s=2:N_h
    L(:,:,s)=max(L(:,:,s-1),L(:,:,s));         % updates L and U
    U(:,:,s)=min(U(:,:,s-1),U(:,:,s));
    TT{s}=U(:,:,s)>=L(:,:,s);                  %checks for non-empty intersection
    if nargin==5 % WOS filtering
        TT{s}=function_WOSFilters(TT{s},N{s-1},O(s-1));      %filtering of the test
        if s==N_h
            TT{s}=function_WOSFilters(TT{s},N{s-1},O(s-1));   %again for the last scale
        end
        TT{s}=TT{s}&TT{s-1};
    end
    YY=TT{s}.*B(:,:,s)+~TT{s}.*YY;  % updates estimates
  %  IN=TT{s}.*repmat(s,N1,N2)+~TT{s}.*IN;  % updates indexes
    IN=TT{s}*s+~TT{s}.*IN;  % updates indexes
%    keyboard
end
% ... end of ICI
if nargin==5 % WOS filtering
    IN=function_WOSFilters(IN,N{N_h},O(N_h));   %filtering of the adaptive scales
    YY=zeros(N1,N2);  % after filtering the scales the estimates need to be recombined
end

IN=medfilt2(IN);  IN(IN==0)=N_h;
std_opt=zeros(N1,N2);
for s=1:N_h
   if nargin==5 % WOS filtering
        YY=YY+(IN==s).*B(:,:,s);   % updates estimates
   end
    std_opt=std_opt+std_B(:,:,s).*(s==IN);
end

%keyboard