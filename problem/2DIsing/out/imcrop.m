clearvars;
dH = 10;
dW = 10;
pad = 7;
figureid = 0;

for fileid = 1:10
    filename_x = sprintf('%s%05d%s', 'fig', fileid, '.png');
    filename_y = sprintf('%s%05d%s', 'fig', fileid+1, '.png');
    BW_x = imread(filename_x);
    BW_y = imread(filename_y);
    %BW = imadjust(BW1);
    fileid
    
    LH = size(BW_x);
xlist = (pad+1) : dW : (LH(1)-dW-pad+1);
ylist = (pad+1) : dH : (LH(2)-dH-pad+1);    
    
    for i = 1:length(xlist)
        for j = 1:length(ylist)
            clf;
            figureid = figureid + 1;
            B_x = BW_x( (xlist(i)-pad) : (xlist(i)+dW+pad-1), (ylist(j)-pad) : (ylist(j)+dH+pad-1));
            B_y = BW_y( xlist(i) : (xlist(i)+dW-1), ylist(j) : (ylist(j)+dH-1));
            montage({B_x, B_y});
            outname_x = sprintf('%s%08d%s', 'output/x_', figureid, '.png');
            outname_y = sprintf('%s%08d%s', 'output/y_', figureid, '.png');
            imwrite(B_x, outname_x);
            imwrite(B_y, outname_y);
         end
    end
end
