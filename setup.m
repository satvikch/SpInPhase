% setup 

%
addpath(strcat(pwd,'\svdandpca'));
addpath(strcat(pwd,'\SpaRSAL'));
addpath(strcat(pwd,'\nlinsar'));
addpath(strcat(pwd,'\data'));
% addpath(strcat(pwd,'\PUMA'));
addpath(strcat(pwd,'\BM3D'));
addpath(strcat(pwd,'\bookdata'));
addpath(strcat(pwd,'\SLC'));
% addpath(strcat(pwd,'\PUMA'));
addpath(strcat(pwd,'\PEARLS'));
addpath(strcat(pwd,'\WFF'));
addpath(strcat(pwd,'\batch_method'));

% add PEARLS to the path.
addpath(strcat(pwd,'\PEARLS\PUMA'));
addpath(strcat(pwd,'\PEARLS\data'));
st = computer; %teste machine
if strcmp(st,'PCWIN64' )
   addpath(strcat(pwd,'\PEARLS\LPA64'));
elseif strcmp(st,'PCWIN' )
   addpath(strcat(pwd,'\PEARLS\LPA32'));
else
   fprintf('\nThis code run only on PC\n')
end
clear st;









    


