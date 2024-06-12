% True Phase models
% Input parameters: test - number of the test function
% coherence parameter for the model 14
% sigma - additive gaussian noise variance
% SNR for the model 15

% Outs: true phase signal - y_true,
% noisy wrapped phase for the models 14 and 15 - z_est;
% otherwises z_est=0;

function [y_true, z_est]=function_TruePhaseModel(test,coherence, sigma,SNR)

      z_est=0;
%%%%%%%%% Quadratic Model %%%%%%%%%%%%%%%%%%%%%%%   
if test==1
% ampl=pi*5; ddelta=pi*5/100;
%         %ampl=2; ddelta=.005*2;
%         [X,Y] = meshgrid(-ampl:ddelta:ampl-ddelta, -ampl:ddelta:ampl-ddelta);
% V2=(X.^2+Y.^2); %0*randn(size(V));
% y_true=V2/2;  y_x=X; y_y=Y; % Quadratic Phase

ddelta=pi*5/10;
        %ampl=2; ddelta=.005*2;
        AMPL=1/10/100/1.2534;
        %[X,Y] = meshgrid(0:255, 0:255);
        [X,Y] = meshgrid(-128:2:127, -128:2:127);  
X=X*ddelta; Y=Y*ddelta;
ddelta=1;
V2=(X.^2+Y.^2); %0*randn(size(V));
y_true=AMPL*V2/30;  
y_true=max(y_true(:))-y_true; y_x=-AMPL*(X-150); y_y=-AMPL*(Y-150); % Quadratic Phase
%[y_x, y_y]=function_Differentiation(y_true,ddelta);

end
%%%%%%%%% Sinc Model %%%%%%%%%%%%%%%%%%%%%%%   
if test==2
         ddelta=2/200;
        %ampl=2; ddelta=.005*2;
       % [X,Y] = meshgrid(-ampl:ddelta:ampl-ddelta, -ampl:ddelta:ampl-ddelta);
         [X,Y] = meshgrid(-128:2:127, -128:2:127);  
         
X=X*ddelta; Y=Y*ddelta;
ddelta=1;
        X=X+0.00001;  Y=Y+0.00001;
        V2=(X.^2+Y.^2);
        frequency=pi*4; 
        ampl=75/5*2/2.7908;
y_true=ampl*sin(frequency*V2)./V2/frequency/2; % 0*randn(size(V));
% y_x=ampl*(frequency*cos(frequency*V2)./V2-sin(frequency*V2)./(V2.^2)).*2.*X/frequency;
% y_y=ampl*(frequency*cos(frequency*V2)./V2-sin(frequency*V2)./(V2.^2)).*2.*Y/frequency;
%[y_x, y_y]=function_Differentiation(y_true,ddelta);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Linear Model %%%%%%%%%%%%%%%%
if test==3
% ampl=1.5; 
ddelta=.01;

        
 %[X,Y] = meshgrid(-128:127, -128:127);  
[X,Y] = meshgrid(0:2:255, 0:2:255);  
X=X*ddelta; Y=Y*ddelta;
ddelta=1;
alpha=25; % [0 .25 .5 .75 1];
y_true=X*alpha; y_x=alpha*ones(size(X)); y_y=0*ones(size(Y)); % Linear phase
%y_true=Y*alpha; %y_x=alpha*ones(size(X)); y_y=0*ones(size(Y)); % Linear phase



%[y_x, y_y]=function_Differentiation(y_true,ddelta);
end

%%%%%%%% Mexican Hat Model %%%%%%%%%%%
if test==4
    
 ddelta=.03;
            
[X,Y] = meshgrid(-128:2:127, -128:2:127);  
X=X*ddelta; Y=Y*ddelta;
     
SIGMA=.8; SIGMA2=SIGMA^2; V2=(X.^2+Y.^2);
AMPL=2*(4*pi/2/pi/SIGMA2^2)/.2383/4/1.4987;            % with '/10' works very well
y_true=AMPL*(2-V2/SIGMA2).*exp(-V2/SIGMA2/2);

%[y_x, y_y]=function_Differentiation(y_true,ddelta);
ddelta=1;
end



%%%%%%%% Gaussian Model %%%%%%%%%%%%%%
if test==6

% SIGMA=.35;
% ddelta=.01;
%  [X,Y] = meshgrid(-128:2:127, -128:2:127);  
% X=X*ddelta; Y=Y*ddelta;
%       ddelta=1;
% 
% 
% SIGMA2=SIGMA^2; V2=(X.^2+Y.^2);
% AMPL1=25/.19305;
% AMPL=AMPL1/2/2; %*pi; %% .25 .5 .75 1 
% 
% y_true=AMPL.*exp(-V2/SIGMA2/2)/2/pi/SIGMA2;

%y_true(1:64,65:end)=0; % Discontinuous model

% 
%y_x=AMPL*(-X/SIGMA2*2).*exp(-V2/SIGMA2/2)+y_true.*(-X/SIGMA2);
% y_y=(AMPL*(-Y/SIGMA2*2).*exp(-V2/SIGMA2/2)+y_true.*(-Y/SIGMA2));

%[y_x, y_y]=function_Differentiation(y_true,ddelta);
M=100;
N=100;
y_true=gaussele(M,N,14*pi,10,15); % large amplitude

y_true(1:50,1:50)=0; % Discontinuous model
end

%%%%%%%% Cosine Model %%%%%%%%%%%%%%
if test==7
ddelta=pi/255*5;

 
%[X,Y] = meshgrid(-ampl:ddelta:ampl, -ampl:ddelta:ampl);  
 [X,Y] = meshgrid(-128:2:127, -128:2:127);  
X=X*ddelta; Y=Y*ddelta;
      ddelta=1;
         
[yN,xN]=size(X);


 V2=sqrt(X.^2+Y.^2);

% y_true=42*(cos(X).*cos(Y));
% y_true=15*(cos(X).*cos(Y));
y_true=5*cos(V2)/2;
image_shift=0;
y_true(:,ceil(xN/2):end)=y_true(:,ceil(xN/2):end)+image_shift;

%[y_x, y_y]=function_Differentiation(y_true,ddelta);



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PIECE-WISE LINEAR MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test==8
%ampl=30;
ddelta=.15;
 
%[X,Y] = meshgrid(0:ddelta:ampl, 0:ddelta:ampl);  
 [X,Y] = meshgrid(-128:2:127, -128:2:127);  
X=X*ddelta; Y=Y*ddelta;
      ddelta=1;
         

         
a1=10/4; a2=-4*2;
y_true=zeros(size(X));
y1=(a1*X+20); y1=y1.*(y1>0);
y2=a2*X+20*2;
y2=y2.*(y2>0);

y12=y1.*(y1<y2)+y2.*(y2<=y1);
y_true=y12;



[y_x, y_y]=function_Differentiation(y_true,ddelta);

end
% 


%%%%%%%%%%%%%%%%%% IMAGES MODELS FOR PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if test==9
yy=im2double(imread('image_Cameraman256.png'));
% yy=im2double(imread('image_Lena256.png'));
%yy=im2double(imread('image_Lena512.png'));
% yy=im2double(imread('image_Cheese128.png'));
% yy=0.8*im2double(imread('image_Cheese128.png'))+0.2;
%yy=im2double(imread('image_Boats512.png'));
%yy=im2double(imread('barbara_port.png'));
%yy=im2double(imread('kodim23.png')); % parrots
% yy=im2double(imread('kodim19.png')); % fence
% yy=im2double(imread('kodim04.png')); % fence
% yy=rgb2Oppon(yy);
%yy=im2double(imread('IMAGE_boats512.png'));
%yy=im2double(imread('testpat1.png'));
%yy=0.8*im2double(imread('testpat1.png'))+0.2;
[n,m]=size(yy);
for s=1:n
    yyy(s,:)=yy(n-s+1,:);
    
end


%yy=im2double(imread('image_Cheese128.png'));
%yy=yy(127:end-127,127:end-127);
y_true=(yyy); 
[X,Y] = meshgrid(0:255, 0:255);  
X=X/255*20; Y=Y/255*20;
      ddelta=1;
         y_true=Y+X+y_true;

ddelta=1; X=yy;

%[y_x, y_y]=function_Differentiation(y_true,ddelta);
end
%%%%%%%%%%%%%%%%%%% Square %%%%%%%%%%%%%%%%%%%%%%%%%%%
if test==10
    y_true(1:2:256,1:2:256)=0;
    
       [X,Y] = meshgrid(-128:2:127, -128:2:127);  

      ddelta=1;
  ampl=.3;       
y_true(84:172,84:172)=1*ampl;  ;
y_true=y_true+image_shift;
[y_x, y_y]=function_Differentiation(y_true,ddelta);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Pyramid %%%%%%%%%%%%%%%%%%%%%%%%%%%

if test==11

alpha=10*pi/126/.9973;
 alpha=alpha*4/2;
[X,Y] = meshgrid(0:2:255, 0:2:255);  
y1=X*alpha;
y2=Y*alpha;
y3=-(X-255)*alpha;
y4=-(Y-255)*alpha;

ddelta=1;
for s1=1:256/2
for s2=1:256/2

y_true(s1,s2)=min([y1(s1,s2) y2(s1,s2) y3(s1,s2) y4(s1,s2)]);


end
end

%[y_x, y_y]=function_Differentiation(y_true,ddelta);


end



%%%%%%%%%% J.M.B.Dias and J.M.N.Leitao models %%%%%%%%%%%%%%%%%%%%%%%%

if test==14 
    
%init=2055615866; randn('seed', init);

%M=100;
M=100;
N=100;
z1=gaussele(M,N,14*pi,10,15); % large amplitude
%[z1, z_est, ddelta]=function_TruePhaseModel(4, coherence, sigma,SNR);
%z1=z1/1.7759; % small amplitude

%z1(1:50,1:50)=0; % Discontinuous model
%z1=ones(size(z1))*3;
[x1 x2] = insarpair(ones(M), coherence*ones(M), z1, sigma^2);	
% estimate interferogram and lambda 
eta  = angle(x1.*conj(x2));	
xx = 1:100;
yy = xx;
ddelta=1;
y_true=z1;
z_est=eta;
ddelta=1;
 ssigma_std=std(wrap(z1(:)-eta(:)))
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if test==15 
    
   % SNR=10;
    
init=2055615866; randn('state', init);

[x,y]=meshgrid(1:128,1:128);
a=128; b=128;
origphase=3-(x-64).^2/80-(y-64).^2/80;
origphase(find(origphase<0))=0;  % formular for phase is truncated parabola


sigma=sqrt(10^(-SNR/10));


complexsig=exp(j*2*pi*(origphase))+sigma*randn(a,b)+i*sigma*randn(a,b);

wrapphase=angle(complexsig)/2/pi;

%%%%%%%%%%  MY variables %%%%%%%%%%%%%%%%%%%%%
y_true=origphase*2*pi;
z_est=wrapphase*2*pi;

ddelta=1;    
end

if test==15 

% demoz3 (z-step) - Ramp with and without discontinuity information
% 

M=150;
N=100;

%build planes
[X,Y] = meshgrid(0:N-1,0:M-1);
%make discontinuous
mask=ones(M,N); mask(1:M/2,:)=0;
X=X.*mask;
% co = 1.0;	%coherence
%generate insarpair
[x1 x2] = insarpair(ones(M,N), co*ones(M,N), X, 0);
eta=angle(x1.*conj(x2));
y_true=X;

z_est=eta;


% set discontinuities
discv = zeros(M,N);
disch = zeros(M,N);
disch(M/2+1,5:N)=1;

%Z-step with know discontinuities
% apha = zstepd(eta,disch,discv); 


end