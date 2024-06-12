% function_Differentiation.m 
function [y_x, y_y]=function_Differentiation(x,ddelta)

% v1= [0 0 0;1/2 0 -1/2;0 0 0];
% v2=-[0 1/2 0; 0 0 0;0 -1/2 0];
v1= [0 0 0;0 1 -1;0 0 0];
v2= [0 1 0; 0 -1 0;0 0 0];

% v1=zeros(11,11);
% v1(6,1)=1/10; v1(6,11)=-1/10;
% v2=zeros(11,11);
% v2(1,6)=-1/10; v2(11,6)=-1/10;
% v_support=3;
% [v1]=function_KernelsForFFT(v1, v_support);
% [v2]=function_KernelsForFFT(v2, v_support);

border0=0;
%x=padarray(x,[border0,border0],'replicate','both');

% v1=padarray(v1,[border0/2,border0/2]);
% v2=padarray(v2,[border0/2,border0/2]);

[yN,xN]=size(x);
size_z_1=yN; size_z_2=xN;
[ghy1,ghx1]=size(v1);

big_v1=zeros(yN,xN); big_v1(1:ghy1,1:ghx1)=v1; 

big_v1=circshift(big_v1,-round([(ghy1-1)/2 (ghx1-1)/2])); % pad PSF with zeros to whole image domain, and centers it.

V1=fft2(big_v1); % Frequency responde of the PSF

big_v2=zeros(yN,xN); big_v2(1:ghy1,1:ghx1)=v2; 
big_v2=circshift(big_v2,-round([(ghy1-1)/2 (ghx1-1)/2])); % pad PSF with zeros to whole image domain, and centers it.

V2=fft2(big_v2); % Frequency responde of the PSF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y_x=fft2(circshift(y_x,-round([(ghy1-1)/2 (ghx1-1)/2])));
% y_y=fft2(circshift(y_y,-round([(ghy1-1)/2 (ghx1-1)/2])));
X=fft2(x);
y_x=real(ifft2(V1.*X));
y_x=y_x(border0+1:end-border0,border0+1:end-border0);
y_y=real(ifft2(V2.*X));
y_y=y_y(border0+1:end-border0,border0+1:end-border0);
%keyboard