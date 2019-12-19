L = 1024;
A = imread('init.png');
A = imresize(A, [L, L]);
B = A(300:700, :, :);
imshow(B)
imwrite(B, 'init_cut.png')
