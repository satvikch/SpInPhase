function A = OMP_C(D,X,param)
%        A = OMP_C(D,X,param)
%
% -- orthogonal matching pursuit to solve an l2-l0 parallel and independent 
%    optimization problems in the complex field.
% 
% -- the original methods was  introduces in 
%

% Y. Pati, R. Rezaiifar, and P. Krishnaprasad, “Orthogonal matching
% pursuit: Recursive function approximation with applications to wavelet
% decomposition,” in Conference Record of The Twenty-Seventh Asilomar
% Conference on Signals, Systems and Computers,, vol. 1, 1993, pp. 40–44.
%
% OMP_C solves approximately  one of the following optimization set of 
%  problems:
%
%
%   a) min_{a_i} ||a_i||  s.t  ||x_i-D*a_i||_2^2 <= err  for i=1,2, ....
%  
%         or
%
%   b) min_{a_i} ||x_i-D*a_i||_2^2  s.t. ||a_i||_0 <= L, for i=1,2,...
%
%
% --- Input parameters ----------------------------------------------------
%
% D   ([m,k] dictionary - complex valued matrix) 
%
% X   ([m,n] Data vectors to code (one per column) - complex valued matrix) 
%      X = [x_1,...,x_n];
%
% ---  input parameters ---------------------------------------------------
%
%     param.L     maximum number of active atoms per vector 
%
%     param.err   tolerance error parameter in ||x_i-D*a_i||_2^2 <= err
%
% --- output --------------------------------------------------------------
%
% A   Sparse codes of each row in X according to dictionary 
%     D: X \simeq D*A + error
%   
% Authors: Hao Hongxing (hongxinghao87@gmail.com) 
%         J. M. Bioucas-Dias (bioucas@lx.it.pt) 
%
% August, 3013


[m,k] = size(D);
n = size(X,2);

if isfield(param, 'L')
    L = param.L;
else
    L = m;
end

if isfield(param, 'err')
    err = param.err;
else
    err = 0;
end


A = zeros(k,n);
A = sparse(A);
Dtx  = @(x) D'*x;
for index = 1:n
    x_i = X(:,index);
    % -- Intitialize --
    % start at a_i = 0, so r_i = x_i - D*a_i = x_i
    a_i = [];
    r_i = x_i;
    err_iter = sum(abs(r_i).^2);
    numberL = 0;
    indexselect = [];
    while (numberL<L && err_iter>err)||(isempty(indexselect))%&&~IsSubMean
        Derr = Dtx(r_i);
        [maxvalue,maxindex] = max(abs(Derr));
        indexselect = [indexselect,maxindex];
        numberL = numberL+1;
        a_i = pinv(D(:,indexselect))*x_i;
        r_i = x_i-D(:,indexselect)*a_i;
        err_iter = sum(abs(r_i).^2);
    end
    if (~isempty(indexselect))
       A(indexselect,index)=a_i;
    end
end
