function [z,iternum, res_p,res_d] = SpaRSAL(M,y,lambda,err)


%% --------------- Description --------------------------------------------
%
%  SpaRSAL solves the following l2-l1 optimization  problem 
%  [size(M) = (L,p); size(X) = (p,N)]; size(Y) = (L,N)]
%
%         min  (1/2) ||M X-y||^2_F + lambda ||X||_1
%          X              
%
%  where ||X||_1 = sum(sum(abs(X)).
% 




% mixing matrixsize
[LM,p] = size(M);
% data set size
[L,N] = size(y);
if (LM ~= L)
    error('mixing matrix M and data set y are inconsistent');
end
% if (L<p)
%     error('Insufficient number of columns in y');
% end


%%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
% maximum number of AL iteration
AL_iters = 100;
% tolerance for the primal and dual residues
tol = err;
%%
%--------------------------------------------------------------
% Local variables
%--------------------------------------------------------------

%---------------------------------------------
%  If lambda is scalar convert it into vector
%---------------------------------------------
lambda = lambda*ones(p,N);

% compute mean norm
norm_y = sqrt(mean(mean(abs(y).^2)));
% rescale M and Y and lambda
M = M/norm_y;
y = y/norm_y;
lambda = lambda/norm_y^2;


%%
%---------------------------------------------
%  Constants and initializations
%---------------------------------------------
mu_AL = 0.01;
mu = 10*mean(lambda(:)) + mu_AL;

% F = M'*M+mu*eye(p);
% [UF,SF] = svd(M'*M);
[UF,SF,VF] = svdecon(M'*M);
sF = diag(SF);
IF = UF*diag(1./(sF+mu))*UF';
% IF = eye(p)/mu-M'*(eye(LM)+M*M'/mu)^(-1)*M/(mu*mu);
% IF = inv(F);
% Aux = IF*B'*inv(B*IF*B');


yy = M'*y;

%%
%---------------------------------------------
%  Initializations
%---------------------------------------------

% no intial solution supplied
x= IF*M'*y;

z = x;
% scaled Lagrange Multipliers
d  = 0*z;

iternum = 0;

%%
%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
tol1 = sqrt(N*p)*tol;
tol2 = sqrt(N*p)*tol;
i=1;
res_p = inf;
res_d = inf;
mu_changed = 0;
while (i <= AL_iters) && ((abs (res_p) > tol1) || (abs (res_d) > tol2)) 
    
    iternum = iternum+1;
    % save z to be used later
    if mod(i,10) == 1
        z0 = z;
    end

%   minimize with respect to z
    %%
    % SpInPhase
%     z =  soft(x-d,lambda/mu); % L1

    %%
    % MC penalty
%     n = 50; 
%     z =  thresh(x-d,lambda/mu, n*(lambda/mu)); % MC

    %%
    % SCAD penalty
    n = 50; 
    z = Sthresh(x-d,lambda/mu, n); % SCAD

    %% 
    % WMCPM
%     z =  soft(x-d,lambda/mu); % L1
%     lambda = max(lambda - (abs(x-d)/n), 0);  % WMCPM

    %%
    % MCpNM
%     p = 0.8; % p value 
%     n = 50;  % n value
%     z =  solve_Lp_w(x-d,lambda/mu, p); % Lp
%     lambda = max(20 - ((abs(x-d).^p)/n), 0); % MCpNM % fix a labmda value based on application/surface.

    %%
    % WMCpNM
%     p = 0.8; % p value
%     n = 50;  % n value
%     z =  solve_Lp_w(x-d,lambda/mu, p); % Lp 
%     lambda = max(lambda - ((abs(x-d).^p)/n), 0);  % WMCpNM

    %%
    x = IF*(yy+mu*(z+d));
    
    % Lagrange multipliers update
    d = d -(x-z);

    % update mu so to keep primal and dual residuals whithin a factor of 10
    if mod(i,10) == 1
        % primal residue
%         res_p = norm(x-z,'fro');
        res_p = norm(abs(x-z),'fro');
        % dual residue
        %res_d = mu*norm(z-z0,'fro');
        res_d = mu*norm(abs(z-z0),'fro');
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            d = d/2;
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;
            d = d*2;
            mu_changed = 1;
        end
        if  mu_changed
           IF = UF*diag(1./(sF+mu))*UF';
%            IF = eye(p)/mu-M'*(eye(LM)+M*M'/mu)^(-1)*M/(mu*mu);
%            Aux = IF*B'*inv(B*IF*B');
            mu_changed = 0;
        end
        
        
    end
    
    i=i+1;
        
   
       
end


    
 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
