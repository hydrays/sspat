L = 512;
A = imread('fibroblast4.png');
A = rgb2gray(A);
%A = A(200:1800-200, 200:2700-200);
B = double(imresize(A, [L, L]));
B = B/max(B(:));
imshow(B);
dlmwrite('matrix_512.txt', B, ' ');