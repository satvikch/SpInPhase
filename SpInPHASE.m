function [image_esti,outputD,alpha, phi_hat_unwrap] = SpInPHASE(Z, param, D)
%        [image_esti,output,alpha, phi_hat_unwrap] = SpInPHASE(Z, param, D)
%
% -- INon-Convex Optimization For Sparse Interferometric Phase Estimation 
%
% Satvik Chemudupati; Praveen Kumar Pokala; Chandra Sekhar Seelamantula, 
% "Non-Convex Optimization For Sparse Interferometric Phase Estimation", 
% IEEE International Conference on Image Processing (ICIP), (published) 2020.
%
%
%
% -- input ---------------------------------------------------------------
%
% Z                            noise data (image of complex numbers)
%
% param                       parameters (stucure)
%  param.patsize              width or height of the patch(usually width=height)
%                             defaul:12
%
%  param.IstrainDic           train the dictionary
%                             (1 - train 0 - do not train)
%                             defaul:1
%
%  param.sigma                noise standard deviation variance of the 
%                             complex noise image:  
%                             either a matrix with the pixels stds   
%                             or a real number with the same std for all
%                             pixels.
%                             default: estimate automatically
%
%  param.mask                 image defining the usefull pixels to be 
%                             estimated 
%                             default: no mask
%
%  param.IsUnwrap             perform unwrapping  with PUMA
%                             (1 - unwrapping 0 - do not unwrapping)
%                             defaul: 0
%
%  param.IsOutputC            0 - output the  estimated wrapped 
%                                 interferometrc phase
%                             1 - output the   complex exponential of the 
%                                 estimated wrapped interferometrc phase
%                             defaul: 0
%
%  param.IsSubMean            Substract the mean
%                             (1 - substract mean of every patch)
%                            (0 - do not substract mean)
%                             defaul: 0
%
%  param.originalphase        original unwrapped phase  used to compute  the 
%                             Mean Square Error and the ISNR 
%                             default: not used
%
% ---  parameters for dictionary learning algorithm: DicLearningM ---------
%
%  param.K                    number of atoms in the trained dictionary
%                             defaul: 256
%
%  param.T                    number of iterations
%                             defaul: 500
%
%  param.patchnum             number of patches sampled on each iteration
%                             defaul: 64
%
%  param.lamda                SparSAL regularization parameter to solve
%                             the BPDN  l2-l1 problems
%                             defaul: 0.11
%
%  param.err                  tolerance error in the SpaSAL algorithm
%                             defaul: 1e-3
%
% ------- optional input dictionary ---------------------------------------
%
%  D                          use a pre-learned dictionary. If D is input
%                             then param.IstrainDic is meaningless 
%
% --------- output --------------------------------------------------------
%
%  image_esti                   either  estimated interferometric phase or
%                               the equivalent unit circle  image of
%                               complexes
%
%  outputD                       learned dictionary
%
%  alpha                        code for the patches based on the dictionary
%
%  phi_hat_unwrap               unwrapped surface if param.IsUnwrap==1
%
% --------------------------------------------------------
%  
% Authors: Hao Hongxing (hongxinghao87@gmail.com) 
%         J. M. Bioucas-Dias (bioucas@lx.it.pt) 
%
% August, 3013
% 





% begin SpInPHASE
% Z is real if it is the input wrapped noisy phase
if imag(Z)==0
    Z = exp(sqrt(-1)*Z);
end

if isfield(param, 'patsize')
    patsize = param.patsize;
else
    patsize = 12;
end

if isfield(param, 'IstrainDic')
    IstrainDic = param.IstrainDic;
else
    IstrainDic = 1;
end

if isfield(param, 'sigma')
    sigmaO = max(param.sigma,0.1);
else
    sigmaO = max(function_stdEst2D_phase(real(Z),1),0.1); 
end
if isfield(param, 'mask')
    mask = param.mask;
else
    mask = ones(size(Z));
end

if size(sigmaO) == size(Z)
    Z = Z./sigmaO;
    sigma = 1;
else
    sigma = sigmaO;
end

if isfield(param, 'IsUnwrap')
    IsUnwrap = param.IsUnwrap;
else
    IsUnwrap = 0; 
end
if isfield(param, 'IsOutputC')
    IsOutputC = param.IsOutputC;
else
    IsOutputC = 0; 
end
if isfield(param, 'IsSubMean')
    IsSubMean = param.IsSubMean;
else
    IsSubMean = 0; 
end
%% initial parameters for dictionary learning
if isfield(param, 'K')
    paramDL.K = param.K;
else
    paramDL.K = 256; 
end
if isfield(param, 'T')
    paramDL.T = param.T;
else
    paramDL.T = 500; 
end
if isfield(param, 'patchnum')
    paramDL.patchnum = param.patchnum;
else
    paramDL.patchnum = 64; 
end
if isfield(param, 'lamda')
    paramDL.lamda = param.lamda;
else
    paramDL.lamda = 0.11; 
end
if isfield(param, 'err')
    paramDL.err = param.err;
else
    paramDL.err = 1e-3; 
end
if isfield(param, 'Dinitial')
    paramDL.Dinitial = param.Dinitial;
end
if isfield(param, 'originalphase')
    originalphase = param.originalphase;
    haveoriginal = 1;
else
    haveoriginal = 0; 
end

if isfield(param, 'L')
    paramOMP.L = param.L;
end

isinputdic = 0;
if(nargin==3)
    if size(D,1)~=patsize*patsize
        error(['The size of the dictionary do not match the patsize']);
    else
        isinputdic = 1;
    end
end


X = im2col(Z,[patsize patsize],'sliding');
maskC = im2col(mask,[patsize patsize],'sliding');
indexm = find(sum(abs(maskC))~=0);
clear maskC;
Xmask = X(:,indexm);
if IsSubMean == 1
    meancolx = angle(sum(Xmask));
    Xmask = Xmask.*(ones(size(Xmask,1),1)*exp(-1i*meancolx));

end

%train the dictionary from the noisy image
if isinputdic==0
    if IstrainDic==1
        [D] = DicLearningM(Xmask,paramDL);
    else
        error(['Please input the dictionary!!!']);
    end
end

outputD = D;

fprintf('Evaluating cost function...\n');
paramOMP.err = chi2inv(0.96,size(Xmask,1)*2)*sigma^2/2;
alpha=OMP_C(outputD,Xmask,paramOMP);



xx_hatmask = outputD*alpha;
if IsSubMean == 1
    xx_hatmask=xx_hatmask.*(ones(size(xx_hatmask,1),1)*exp(1i*meancolx));
end


X(:,indexm) = xx_hatmask;

blockx = patsize;
blocky = patsize;
final_numestimate = zeros(size(Z));
final_extentestimate = zeros(size(Z));
for indexi = 1:blocky
    for indexj = 1:blockx
        tempesti = reshape(X((indexi-1)*blockx+indexj,:),size(Z)-[blockx,blocky]+1);
        numestimate = zeros(size(Z));
        extentestimate = zeros(size(Z));
        extentestimate(1:size(tempesti,1),1:size(tempesti,2)) = tempesti;
        numestimate(1:size(tempesti,1),1:size(tempesti,2)) = 1;
        
        extentestimate = circshift(extentestimate, [indexj,indexi]-1);
        numestimate = circshift(numestimate, [indexj,indexi]-1);
        
        final_numestimate = final_numestimate+numestimate;
        final_extentestimate = final_extentestimate+extentestimate;
    end
end
image_esti_c = final_extentestimate./final_numestimate;
%%
image_esti = angle(image_esti_c);
phi_hat_unwrap = zeros(size(image_esti));
if IsUnwrap==1
    phi_hat_unwrap = puma_ho(image_esti ,0.5, 'verbose','no' );
end

if haveoriginal==1
    wraperr_norm = norm(wrap(image_esti - originalphase),'fro')^2;
    RMSE_SpInPHASE = sqrt(wraperr_norm/M/N);
    PSNR_SpInPHASE = 10*log10(4*M*N*pi^2/wraperr_norm);
    fprintf('\n\n Performance indices\n')
    fprintf('RMSE (SpInPHASE) = %2.3f\n', RMSE_SpInPHASE)
    fprintf('PSNR (SpInPHASE) = %2.3f\n', PSNR_SpInPHASE)
end
if IsOutputC==1
    image_esti = image_esti_c;
end
