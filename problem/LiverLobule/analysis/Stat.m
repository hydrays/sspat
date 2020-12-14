Nsmall = 4;
for timepoint = 1:3
    timepoint
    statfilename = sprintf('%s%02d%s', 'stat/stat_size', timepoint, '.csv');  
    fileID = fopen(statfilename,'w');
    fprintf(fileID,'%6s, %12s\n','id','size');
    for rep = 1:20
        rep
        filename = sprintf('%s%02d%s%04d%s', 'stat/timepoint', timepoint, '/brick_', rep, '.tif');
        A = imread(filename);
        imshow(A);
        CC = bwconncomp(A);
        for i = 1:CC.NumObjects
            if (length(CC.PixelIdxList{i}) > Nsmall)
                fprintf(fileID,'%6d, %12.1f\n', i, ...
                        length(CC.PixelIdxList{i}));
            end
        end
    end
    fclose(fileID);
end
