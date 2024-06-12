function [wphf] = denoiseLPA_Varying_Window_Sizes(wph,h_opt,mask)
%        wphf = denoiseLPA(wph,n1,n2)
% 
%  Denoise  a wraped interferogram by adjusting local
%  first order polinomials  (planes). See 
%
%  J. Bioucas-Dias, V. Katkovnik,  J. Astola, and  K.  Egiazarian, 
%  "Absolute phase estimation: adaptive local denoising and global unwrapping",  
%  vol. 47. no, 29, pp. 5358-5369, Applied Optics, 2008.
%
%
%  %%%%%%% Input parameters %%%%%%%%%%%%%%%%
%
%  wph     -- wrapped noisy phase
%  ny * nx -- window size support for  local polinomials
%  
%   
%  %%%%%%% Output  parameters %%%%%%%%%%%%%%%%
%
%  wphf -- denoised phase
%
%
% Authors: 
% V. Katkovnik, Tampere University of Technology, Tampere, Finland,
% J. Bioucas-Dias, Instituto Superior Tecnico, Lisbon, Portugal
% October 2007

disp('-----------------------------------------')
disp('Varying_Window Size First Order Filtering')
disp('-----------------------------------------')

I= sqrt(-1);
fft_padd = 64;
% data size
[M,N] = size(wph);

wphf = zeros(M,N);
for i=1:M
    for j= 1:N  
        ny=2*h_opt(i,j)+1;  %%% ADAPTIVE WINDOW SIZE SUPPORTS
        nx=ny;
% define central points
o_y = ceil((ny-1)/2)+1;
o_x = ceil((nx-1)/2)+1;

Mp= M + 2*(o_y -1);
Np= N + 2*(o_x -1);

% padded data size

% insert a frame on wph
wph_padded = zeros(Mp,Np);
mask_padded=zeros(Mp,Np);
%keyboard
% copy wrapped phase
wph_padded(o_y:(Mp-o_y+1), o_x:(Np-o_x+1)) = wph;
mask_padded(o_y:(Mp-o_y+1), o_x:(Np-o_x+1)) = mask;
data_window = exp(I*wph_padded(i:i+ny-1,j:j+nx-1)).*mask_padded(i:i+ny-1,j:j+nx-1); %% Original size data
               
       % padd with zeros 
       data_window(fft_padd,fft_padd) =0;  %% padded data
       % compure 2 FFT
       %%%%%%%%     MATRIX FFT2 %%%%%%%%%%%%%%
       % fdata=conj(tran_fft_matrix')*data_window*tran_fft_matrix; %% fdata novel matrix FFT2
        fdata = fft2(data_window);
       %compute maximum  of abs(fdata)
       bb=abs(fdata);
       %bb=sqrt(cos(fdata).^2+sin(fdata).^2);
       [dummy colm] = max(max(bb));
       [mse(i,j) linem] = max(bb(:,colm));
%%%%%%%%%%%%%% Alternative maximization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%        dummy=abs(fdata);
%        [linem colm] = find(dummy==max(dummy(:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % get co and correct for  originoffset
       wphf(i,j) = angle(fdata(linem,colm)*exp(I*(2*pi)/fft_padd*((o_x-1)*(colm-1)+(o_y-1)*(linem-1))));
       % here we may compute c1 and c2 from colm and linem
       % mse per pixel
       %mse(i,j) = 2-2*mse(i,j)/nx/ny;
       if h_opt(i,j)==0
           wphf(i,j)=wph(i,j);
       end
     
    end
end
