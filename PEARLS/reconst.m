function  apha = reconst(epha,qual,mask,disch,discv,iter)
%   apha = reconst(epha,mask,iter);
%   Reconstruct phase in ~qual&mask areas by minimizing the 
%   least square error differences using ICM, given the estimated 
%   phase at the boundary
%   eapha    - estimated absolute  aphase in mask pixels
%   iter     - number of ICM iterations 
%   qual     - float matrix with mask (1 - don't visit; 0-visit)
%   mask     - pixls outside mask have no meaning
%   discv    - vertical discontinuity matrix 
%              disch(i,j) = 1 means a horizontal discontinuity 
%                           between site (i,j) and site (i,j-1);
%                           disch(i,j) \in [0,1] 
%   disch    - horizontal discontinuity matrix 
%              discv(i,j) = 1 means a horizontal discontinuity 
%                           between site (i,j) and site (i-1,j)
%                           discv \in [0 1]
%
%   Author J.M. Bioucas Dias 2000
%   Topic - Interferometry




[M N] = size(epha);
apha = epha;

for m=1:iter
   for nextsite=1:M*N  
      if ((~qual(nextsite))& mask(nextsite))
         % find nextsite neighbors 
         c= floor((nextsite-1)/M)+1;  % column
         l= nextsite - (c-1)*M;       % line         
         % set 1  moreneighs(i) to 1 if there is no discontinuity
         % to the correspondent neighbor
         neighs = [0 0 0 0]; %[(l,c+1) (l-1,c)(l,c-1)(l+1,c)]
         discweight = [0 0 0 0]; 
         if ((c+1) <= N) & mask(l,c+1)
            neighs(1) = nextsite+M;
            discweight(1) = 1 - discv(l,c+1);
         end
         if ((l-1) >= 1) & mask(l-1,c)
            neighs(2) = nextsite-1;
            discweight(2) = 1 - disch(l,c);
         end
         if ((c-1) >= 1) & mask(l,c-1)
            neighs(3) = nextsite-M;
            discweight(3) = 1 - discv(l,c);
         end
         if ((l+1) <= M) & mask(l+1,c)
            neighs(4) = nextsite+1;         
            discweight(4) = 1 - disch(l+1,c);
         end
         discsum = sum(discweight);
         meanpha = 0;
         for j=1:4
            if neighs(j) ~= 0
               meanpha  = meanpha  + apha(neighs(j))*discweight(j);
            end
         end
         if discsum ~= 0.0
             %icm steep
            apha(nextsite) = meanpha /discsum;
         end
      end % end if
   end  % end next site
end % iter
