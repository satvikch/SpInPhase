%% train dictionary by solving the problem
%
%   min_{D,alpha} ||x-D*alpha||_2^2+lamda*||alpha||_1
%

%%                                    input
% X samples organized in column to train the dictionary

%
%%%%%%%%%%%%%%%%%%%%%%%%input parameters for dictionary learning%%%%%%%%%%%
%  param.K                    number of atoms in the trained dictionary
%                             defaul: 256
%  param.T                    number of iterations in the trained dictionary
%                             defaul: 500
%  param.patchnum             number of patches use in one ieration in the trained dictionary
%                             defaul: 64
%  param.lamda                parameter in SUNSAL to solve the 
%                             defaul: 1.8
%  param.err                  tolerant error in the SUNSAL algorithm
%                             defaul: 1e-2

%%%%%%%%%%%%%%%%%%%%%%%output%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  D                         trained dictionary 
%  A, B                      information record in the iteration of training dictionary.


function [D,energy,Etime,alphat,A,B] = DicLearningM_batch(X,param)

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
% if isfield(param, 'TT')
%     TT = param.TT;
% else
%     TT = 100; 
% end
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
A = zeros(k,k);
B = zeros(m,k);
A = t0*eye(k);

% Aenergy = zeros(k,k);
% Benergy = zeros(m,k);
% Xenergy = zeros(m,m);
% sumalphan1 = 0;


DCT=zeros(m,k);
for ii=0:1:k-1,
    V=cos([0:1:m-1]'*ii*pi/k);
    if ii>0, V=V-mean(V); end;
    DCT(:,ii+1)=V/norm(V);
end;
D = DCT;
% D = rand(m,k);
D=D./repmat(sqrt(sum(abs(D).^2)),[size(D,1) 1]);
% load('Dini.mat');
% load('D_nonoise.mat');

B = t0*D;


A = sparse(A);
B = sparse(B);


X = X./repmat(sqrt(sum(X.*conj(X))),size(X,1),1);

% S=2055615866; randn('seed', S);
% D=randn(m,k)+sqrt(-1)*randn(m,k);
% D=D./repmat(sqrt(sum(abs(D).^2)),[size(D,1) 1]);

% rng('default');
% Sele = randi(size(X,2),1,param.T);
% patchnum = param.patchnum;
% rng('default');
% rng('shuffle');
% Sele = randi(size(X,2),1,T*patchnum);
% Alpha_all = zeros(k,size(X,2));

energy = zeros(T,1);
Etime = zeros(T,1);
Xi = X;

% index = Sele(1,1:T*patchnum);
% Xi = X(:,index);
    
for ite = 1:T
%     tic;    
%     indeend = min(Sele(1,ite)+patchsize,size(X,2));
%     Xi = X(:,Sele(1,ite):indeend);


%     Xi = X(:,Sele(1,ite));
    
    %creat initial dictionary
    
%     alpha0 = zeros(k,size(Xi,2));    
%     alphat=mexL1L2BCD(Xi,D,alpha0,listgroups,parambcd);

%     lambda = param.lamda;
%     alphat = sunsal(D,Xi,'lambda',lambda,'verbose','yes');
%     alphat = sunsal(D,Xi,'lambda',lambda);
    alphat = SpaRSAL(D,Xi,lambda,err);
%     [alphat, Xvar, paramm] = EMBGAMP(Xi, D) ;
    alphat(abs(alphat)<1e-4) = 0;
    alphat = sparse(alphat);
%     energy(ite) = norm(D*alphat-Xi,'fro')^2+lambda*norm(alphat,1);
    
%     delta=sqrt(kais)*param.sigma;
%     q = 1;
%     [alphat,res1,res2,res3] = csunsal_p(D,Xi,'Q',q,'delta',delta,'ADDONE','no',...
%         'POSITIVITY','no','AL_iters',100,'verbose','no');
%     alphat(abs(alphat)<1e-4) = 0;
    
%     for tempindex = 1:patchsize
%         A = A+alphat(:,tempindex)*alphat(:,tempindex)';
%         B = B+Xi(:,tempindex)*alphat(:,tempindex)';
%         for j = 1:k
%             if A(j,j)~=0
%                 uj = (B(:,j)-D*A(:,j))/A(j,j)+D(:,j);
%                 dj = uj/norm(uj);
%                 D(:,j) = dj;
%             end
%         end
%     end    
    
    %%
%     if ite<patchnum
%         theta = ite*patchnum;
%     else
%         theta = patchnum*patchnum+ite-patchnum;
%     end
%     beta = (theta+1-patchnum)/(theta+1);
%%                  
%     rho = 2;
%     beta = (1-1/ite)^rho;
%     beta = 1;


    
    Atemp = zeros(k,k);
    Btemp = zeros(m,k);
    Atemp = sparse(Atemp);
    Btemp = sparse(Btemp);
    for tt = 1:size(Xi,2)
        Atemp = Atemp+alphat(:,tt)*alphat(:,tt)';
    end
    A = Atemp;
%     A = beta*A+Atemp;
    for tt = 1:size(Xi,2)
        Btemp = Btemp+Xi(:,tt)*alphat(:,tt)';
    end
    B = Btemp;
%     B = beta*B+Btemp;
    
    %%                     variable for compute the energy
    Xtemp = zeros(m,m);
    for tt = 1:size(Xi,2)
        Xtemp = Xtemp+Xi(:,tt)*Xi(:,tt)';
    end
    
%     Xenergy = Xtemp;
%     Aenergy = Atemp;
%     Benergy = Btemp;
%     sumalphan1 = lambda*sum(abs(alphat(:)));
    
%     Xenergy = beta*Xenergy+Xtemp;
%     Aenergy = beta*Aenergy+Atemp;
%     Benergy = beta*Benergy+Btemp;
%     sumalphan1 = beta*sumalphan1+lambda*sum(abs(alphat(:)));

%     Xenergy = Xtemp;
%     Aenergy = Atemp;
%     Benergy = Btemp;
%     sumalphan1 = sum(abs(alphat(:)));
    %%
    
%     A = A+alphat*alphat';
%     B = B+Xi*alphat';
    
%     A = A+alphat*alphat'-Alpha_all(:,Sele(1,ite))*Alpha_all(:,Sele(1,ite))';
%     B = B+Xi*alphat'-Xi*Alpha_all(:,Sele(1,ite))';    
%     Alpha_all(:,Sele(1,ite)) = alphat;
    
    
    for j = 1:k
        if A(j,j)~=0
            uj = (B(:,j)-D*A(:,j))/A(j,j)+D(:,j);
            dj = uj/norm(uj);
%             dj = uj/max(norm(uj),1);
            D(:,j) = dj;
        end
    end
    
%     for ttt = 1:5
%         for j = 1:k
%             if A(j,j)~=0
%                 uj = (B(:,j)-D*A(:,j))/A(j,j)+D(:,j);
%                 dj = uj/max(norm(uj),1);
%                 D(:,j) = dj;
%             end
%         end
%     end
%     for j = 1:k        
%         D(:,j) = D(:,j)/norm(D(:,j));        
%     end
        fprintf('Iteration: %d\n',ite);
% %          energy(ite) = 1/(ite*patchnum)*(1/2*trace(D*Aenergy*D')-1/2*trace(Benergy*D')-1/2*trace(D*Benergy'));
%         energy(ite) = 1*(1/2*trace(D*Aenergy*D')-1/2*trace(Benergy*D')-1/2*trace(D*Benergy')+1/2*trace(Xenergy)+full(sumalphan1))/(size(Xi,2));
% %         eenergy(ite) = 0.5*norm(Xi-D*alphat,'fro').^2+lambda*(sum(sum(abs(alphat))));
%         Etime(ite) = toc;
        
%         energy(ite) = 1/ite*(1/2*trace(D*A*D')-1/2*trace(B*D')-1/2*trace(D*B'));
%         if ite == 1
%                 energy(ite) = 1/(patchnum)*(1/2*trace(D*Aenergy*D')-1/2*trace(Benergy*D')-1/2*trace(D*Benergy')+1/2*trace(Xenergy)+lambda*sumalphan1);
%         else
%             energy(ite) = 1/(ite*patchnum)*(energy(ite-1)*(ite-1)*patchnum+1/2*trace(D*Aenergy*D')-1/2*trace(Benergy*D')-1/2*trace(D*Benergy')+1/2*trace(Xenergy)+lambda*sumalphan1);
%         end
%         energy(ite) = 1/(ite*patchnum)*(1/2*trace(D'*D*Aenergy)-trace(D'*Benergy)+1/2*trace(Xenergy)+lambda*sumalphan1);
%         energy(ite) = 1/ite*(1/2*trace(D'*D*A)-1/2*trace(B*D')-1/2*trace(D*B'));
end

D=D./repmat(sqrt(sum(abs(D).^2)),[size(D,1) 1]);

