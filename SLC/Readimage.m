function re = Readimage(filename,small,large)
A = imread(filename);
AA=A(:,:,1);
AA = double(AA);
re = AA/255*(large-small)+small;