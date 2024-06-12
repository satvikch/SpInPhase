function  mf = fatten(mask)
%         mf = fatten(mask);
%
%         Fatten mask to first order neifgbors
%         Note: fatten site if mask(site) == 1

[M N] = size(mask);
mf=mask;

for nextsite=1:M*N  
   if (mask(nextsite) == 1)
      % find nextsite neighbors 
      c= floor((nextsite-1)/M)+1;  % column
      l= nextsite - (c-1)*M;       % line         
      if ((c+1) <= N)    %[(l,c+1)
         mf(l,c+1) = 1;
      end
      if ((l-1) >= 1)    %(l-1,c)
         mf(l-1,c) = 1;
      end
      if ((c-1) >= 1)    %(l,c-1)
            mf(l,c-1) = 1;
      end
      if ((l+1) <= M)    %(l+1,c)
         mf(l+1,c) = 1;
      end
   end
end
