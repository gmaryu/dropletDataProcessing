function retv = detectMultiNuclei(maskMat)
%%
    thresframe = 10;
    cnt = 0;
    retv = nan;
    mask = load(maskMat).("nuclearMask");
    for i = 1:size(mask, 3)
        im = mask(:, :, i);
        cc = bwconncomp(im);

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
    end
end

