clear all
close all

I = sqrt(-1);
%%
% % size of the image
% M=100;
% N=100;
% 
% % horizintal ramp
% 
% [X,Y] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
% varphi = (X);
% 
% % interferogram
% z = exp(I*varphi);
% % compute the interferogram;
% phi = angle(z);
%%
M=150;
N=100;

%build planes
[X,Y] = meshgrid(1:N,1:M);
%make discontinuous
mask=ones(M,N); mask(:,1:N/2)=0;
varphi=1*X;
% varphi=X.*mask;

clear X Y;
z = exp(I*varphi);
phi = angle(z);

% zero padding
% window size (n = (2*ws+1)^2)
ws = 1;
bloc_size = 16; % ti be incremente by zero padding

z_pad = zeros(M+2*ws,N+2*ws);
z_pad (ws+1:end-ws,ws+1:end-ws ) = z;

phi_hat = zeros(M,N);
mag = zeros(M,N);
for i=ws+1:(M+2*ws)-ws
    for j=ws+1:(N+2*ws)-ws
        block = z_pad(i-ws:i+ws, j-ws:j+ws);
        block(bloc_size,bloc_size) = 0;
        Fz = fft2(block);
        [max_F, index] = max(abs(Fz(:)));
        % conver linear index to array subscript
        [k,l]=ind2sub([bloc_size,bloc_size],index);
        ky = 2*pi/bloc_size*(k-1);
        kx = 2*pi/bloc_size*(l-1);
        % correct for the center of the window
        phi_hat(i-ws,j-ws) = angle( Fz(index)*exp(I*ky*ws+I*kx*ws));
        mag(i-ws,j-ws) = max_F;
    end
end

figure(1)
imagesc(wrap(phi_hat-varphi))


[phi_hat_hao,c1_hat,c2_hat, mag_hat] = denoiseLPAC(phi,ws*ones(M,N),bloc_size);
figure(2)
imagesc(wrap(phi_hat_hao-varphi))


