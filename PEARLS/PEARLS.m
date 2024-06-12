function [phi_hat, index_h] = PEARLS(phi, varargin)
%        [wphi_hat, index_h] = PEARLS(wphi, windows, gama)
%
%  PEARLS -- phase estimation using adaptive regularization
%            based on local smoothing
%
%
%  This function implements the absolute phase estimation algorithm
%  introduced in
%
%  J. Bioucas-Dias, V. Katkovnik,  J. Astola, and  K.  Egiazarian,
%  "Absolute phase estimation: adaptive local denoising and global unwrapping",
%  vol. 47. no, 29, pp. 5358-5369, Applied Optics, 2008.
%
%  The unwrapping step is computed with the PUMA algorithm introduced in
%
%  J. Bioucas-Dias  and  G. Valadão, "Phase unwrapping via graph cuts",
%  IEEE Transactions on Image processing, vol. 16, no. 3, pp. 698-709, 2007.
%
%
%
%%  ===== Required inputs =============
%
%   phi  -- noisy  wrapped phase  inaage in the interval [-pi,pi)
%
%%  ===== Optional inputs =============
%
%  'windows'  --  set of local window sizes to implement local polynomial
%                 approximations. The size of the windows is to be selected
%                 by the ICI method.
%                 Default: [1 2 3 4]
%
%  'gamma'    -- ICI parameter influencing the adaptive window  selection
%                mechanism.  Higher values of gamma yields larger windows.
%                Default: 2.0;
%
%  'unwrapp'  -- {'no', 'yes'} weather or not the unwrapping step is
%                applied
%                Default: 'no'
%
%  'PUMA_pot' -- (struct) type of PUMA potential
%                Default: potential.quantized = 'no'
%                         potential.threshold = 2;
%
%
%  'PUMA_exp' -- Exponent of the potential in PUMA
%                Default: 1;
%                
%
%
%  'verbose'   = {'yes', 'no'};
%                 'no' - work silently
%                 'yes' - display warnings
%                  Default 'no'
%
%%  ===== Output =============
%
%  phi_hat  -- estimate phase
%
%  index_h  -- indexes for the estimated windows
%
%
%
% Authors: Jose M. Bioucas Dias  and Vladimir Katkovnil, 2012

%%  begin


%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 1
    error('Wrong number of required parameters');
end



%%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
% window sizes
HH = [1 2 3 4];
% ICI parameter
gamma = 2;
%unwrapp
unwrapp = 'no';
% verbose
verbose = 'no';
% PUMA parameters
potential.quantized = 'no';
potential.threshold = 2;
p_exponent = 1;



%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'WINDOWS'
                HH = varargin{i+1};
            case 'GAMMA'
                gamma = varargin{i+1};
            case 'PUMA_POT'
                potential = varargin{i+1};
            case 'PUMA_EXP'
                p_exponent = varargin{i+1};
            case 'UNWRAPP'
                unwrapp = varargin{i+1};
            case 'VERBOSE'
                verbose = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end


% stuff
I = sqrt(-1);

% ensures tha the interferogram is in [-pi,pi]
phi = angle(exp(I*phi));

%%
%select windows using ICA based on zero-order LPA
if verbose
    disp('------------------------')
    disp('LPA ZERO ORDER filtering')
    disp('------------------------')
end
mask = ones(size(phi));
[phi_hat, index_h, std_h]=LPA_ICI_ZeroOrder(phi,HH,mask,gamma);

h_opt = HH(index_h); % Adaptive window sizes
% index_h - adaptive indexes


% LPA first order approximation
if verbose
    disp('------------------------')
    disp('LPA FIRST ORDER filtering')
    disp('------------------------')
end
%[phi_hat] = denoiseLPA_Varying_Window_Sizes(phi,h_opt,mask);
% C implementation
[phi_hat,c1_hat,c2_hat,mag_hat] = denoiseLPAC(phi,h_opt,16);

%%
% if verbose
%     disp('------------------------')
%     disp('LPA First ORDER filtering')
%     disp('------------------------')
% end
% mask = ones(size(phi));
% [phi_hat, index_h, std_h]=LPA_ICI_FirstOrder(phi,HH,mask,gamma);

%%


if strcmp(unwrapp,'yes')
    if verbose
        disp('------------------------')
        disp('UNWRAPPING')
        disp('------------------------')
    end
    % unwrapp  the denoise interferogram (absolute phase reconstruction)
    potential.quantized = 'no';
    potential.threshold = 2;
    phi_hat = puma_ho(phi_hat,p_exponent, 'potential',potential );
end





