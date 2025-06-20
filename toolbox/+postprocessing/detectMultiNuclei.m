function [maxObjCount, objCountSeries] = detectMultiNuclei(maskMat)
%% test mode
saveMultiNuc = false;
savePath = 'E:\MATAB_NC_project\exports\test_sperm/multiNuc';
if ~exist(savePath, 'dir')
    mkdir(savePath);
end
%%
    thresframe = 10;
    cnt = 0;
    maxObjCount = 1;
    objCountSeries = [];
    if contains(maskMat, 'dna')
        mask = load(maskMat).("dnaMask");
   
    elseif contains(maskMat, 'nuclear')
        mask = load(maskMat).("nuclearMask");

    else
        fprintf('invalid file type.')
        return
    end
    
    for i = 1:size(mask, 3)
        im = mask(:, :, i);
        % bw = imbinarize(im);
        % D = bwdist(~bw);
        % D_neg = -D;
        % L = watershed(D_neg);
        % bw_separated = bw;
        % bw_separated(L == 0) = 0;
        % cc = bwconncomp(bw_separated);
        cc = bwconncomp(im);
        numObjects = cc.NumObjects;
        objCountSeries = [objCountSeries, numObjects];
        %{
        if saveMultiNuc && numObjects > 0
            disp(numObjects);
            fn = sprintf("multiNuc_%d.png",i);
            imwrite(bw_separated, fullfile(savePath,fn));
        end
        %}
 
        if numObjects > maxObjCount 
           %disp(numObjects);
           maxObjCount = numObjects;
        end
        %{
        if cc.NumObjects == 1
            retv = 1;
            
        elseif cc.NumObjects > 1
            lens = sort(cellfun('length', cc.PixelIdxList));
            lcclen = lens(end);
            secondlen = lens(end - 1);
            if secondlen > lcclen * 0.75
                cnt = cnt + 1;
                if cnt >= thresframe
                    retv = 2;
                    return
                end
            end
        end
        %}
    end

end

