function [D,A,B] = DicLearningM(X,param)
%       [D,A,B] = DicLearningM(X,param)
%
% Online dictionary learning method for the complex numbers by solving 
% the optimization problem
%
%   min_{D,alpha} ||x-D*alpha||_2^2+lamda*||alpha||_1
%
% where D and alpha are matrices
%
% -- the original methods was introduced in 
%
% J. Mairal, F. Bach, J. Ponce, and G. Sapiro, “Online dictionary learning
% for sparse coding,” in Proceedings of the 26th Annual International
% Conference on Machine Learning, ser. ICML ’09. New York, USA:
% ACM, 2009, pp. 689–696.
% 
% ------------------------input parameters---------------------------------
% 
% X sample vectors (columns) used to train the dictionary
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
%------------------------output--------------------------------------------
% 
%  outputD                   learned dictionary
%
%  A, B                      information record along 
%                            the iterations (see ref above)
%
% --------------------------------------------------------
%  
% Authors: Hao Hongxing (hongxinghao87@gmail.com) 
%         J. M. Bioucas-Dias (bioucas@lx.it.pt) 
%
% August, 3013
% 




if isfield(param, 'K')
    k = param.K;
else
    k = 256; 
end
if isfield(param, 'T')
    T = param.T;
else
    T = 500; 
end
if isfield(param, 'patchnum')
    patchnum = param.patchnum;
else
    patchnum = 64; 
end
if isfield(param, 'lamda')
    lambda = param.lamda;
else
    lambda = 0.11; 
end
if isfield(param, 'err')
    err = param.err;
else
    err = 1e-3; 
end


t0 = 1e-3;
m = size(X,1);
% A = zeros(k,k);
% B = zeros(m,k);
A = t0*eye(k);

if isfield(param, 'Dinitial')
    D = param.Dinitial;
else
    DCT=zeros(m,k);
    for ii=0:1:k-1,
        V=cos([0:1:m-1]'*ii*pi/k);
        if ii>0, V=V-mean(V); end;
        DCT(:,ii+1)=V/norm(V);
    end;
    D = DCT;
    D=D./repmat(sqrt(sum(abs(D).^2)),[size(D,1) 1]);
end


B = t0*D;


A = sparse(A);
B = sparse(B);

X = X./repmat(sqrt(sum(X.*conj(X))),size(X,1),1);
Sele = randi(size(X,2),1,T*patchnum);

for ite = 1:T
    
    index = Sele(1,(ite-1)*patchnum+1:ite*patchnum);
    Xi = X(:,index);  
    alphat = SpaRSAL(D,Xi,lambda,err);
    alphat(abs(alphat)<1e-4) = 0;
    alphat = sparse(alphat);

    rho = 2;
    beta = (1-1/ite)^rho;
    
    
    Atemp = zeros(k,k);
    Btemp = zeros(m,k);
    Atemp = sparse(Atemp);
    Btemp = sparse(Btemp);
    for tt = 1:patchnum
        Atemp = Atemp+alphat(:,tt)*alphat(:,tt)';
    end
    A = beta*A+Atemp;
    for tt = 1:patchnum
        Btemp = Btemp+Xi(:,tt)*alphat(:,tt)';
    end
    B = beta*B+Btemp;

    for j = 1:k
        if A(j,j)~=0
            uj = (B(:,j)-D*A(:,j))/A(j,j)+D(:,j);
%             dj = uj/norm(uj);
            dj = uj/max(norm(uj),1);
            D(:,j) = dj;
        end
    end
    fprintf('Iteration: %d\n',ite);
end

D=D./repmat(sqrt(sum(abs(D).^2)),[size(D,1) 1]);

