function  [dh, dv] = linefg(mask);
%         [dh, dh] = linefg(mask);
%
%         Generate  line field discontinuities from mask
%         if mask(i) == 0 set discontinuities

[M N] = size(mask);
dh=zeros(M,N);
dv=zeros(M,N);

for nextsite=1:M*N  
   if (mask(nextsite) == 0)
      % find nextsite neighbors 
      c= floor((nextsite-1)/M)+1;  % column
      l= nextsite - (c-1)*M;       % line         
      if ((c+1) <= N)    %[(l,c+1)
         dv(l,c+1) = 1;
      end
      if ((l-1) >= 1)    %(l-1,c)
         dh(l,c) = 1;
      end
      if ((c-1) >= 1)    %(l,c-1)
            dv(l,c) = 1;
      end
      if ((l+1) <= M)    %(l+1,c)
         dh(l+1,c) = 1;
      end
   end
end
