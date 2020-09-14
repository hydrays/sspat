for i = 0 : 1 : 2000
    i
    ss = sprintf('%s%05d%s', 'c', i, '.dat');
    ssout = sprintf('%s%05d%s', 'fig', i, '.png');
    A = readmatrix(ss);
    A = 0.5*(A + 1);
    imwrite(A, ssout);
end
