function [hoechstsum, npts, smooththres, smoothbg, idxToFrame] = sumHoechstIntwNucMask(files, nucSegOutputMat, output)
% sumHoechstIntwNucMask Quantifies Hoechst intensity within nuclear masks.
%
%   [hoechstsum, npts, smooththres, smoothbg, idxToFrame] = sumHoechstIntwNucMask(files, nucSegOutputMat, output)
%
% This function reads a series of Hoechst-stained TIFF images defined by 'files' and uses
% nuclear segmentation results from 'nucSegOutputMat' to compute, for each frame, the total
% Hoechst intensity (above background) and the number of pixels above threshold.
% Results are saved to the specified output file.
%
% Inputs:
%   files           - (1,1) string. File pattern for Hoechst images, e.g., "Pos6_DAPI_???.tif".
%   nucSegOutputMat - (1,1) string. Filename of the .mat file containing nuclearMask and idxToFrame.
%   output          - (1,1) string. Output filename to save the results.
%
% Outputs:
%   hoechstsum - 1×N vector of total Hoechst intensity per frame.
%   npts       - 1×N vector of number of pixels above the threshold per frame.
%   smooththres- 1×N vector of threshold values computed from the nuclear region.
%   smoothbg   - 1×N vector of background intensity values.
%   idxToFrame - 1×N vector mapping image index to frame number.
%
% Example:
%   [hsum, npts, sthres, sbg, idx] = postprocessing.sumHoechstIntwNucMask("Pos6_DAPI_???.tif", "nuclearSeg.mat", "output.mat");

    arguments
        files (1,1) string
        nucSegOutputMat (1,1) string
        output (1,1) string
    end

    radiusMargin = 0.15 + 0.05;
    fs = dir(files);
    N = length(fs);
    if N == 0
        error("No files found matching pattern: %s", files);
    end
    root = fs(1).folder;
    
    % Load nuclear segmentation results.
    nucData = load(nucSegOutputMat);
    nuclearMask = nucData.nuclearMask;
    idxToFrame = nucData.idxToFrame;
    
    if size(nuclearMask,3) ~= N
        error("Frame count mismatch: nuclearMask and image files differ.");
    end
    
    rawImages = cell(1, N);
    bgs = zeros(1,N);
    ts = zeros(1,N);
    for i = 1:N
        rawImages{i} = imread(fullfile(root, fs(i).name));
        if i == 1
            nPixels = size(rawImages{1}, 1);
            [xx, yy] = meshgrid(1:nPixels, 1:nPixels);
            radius = nPixels / 2 / (1 + radiusMargin);
            mask = uint16(hypot(xx - (nPixels+1)/2, yy - (nPixels+1)/2) < radius);
            nanmask = double(mask);
            nanmask(nanmask==0) = NaN;
        end
        
        rm = rawImages{i};
        mask_current = bitand(imdilate(nuclearMask(:,:,i) > 0, strel("disk", round(radius*0.1))), ~isnan(nanmask));
        mask_bg = bitand(nuclearMask(:,:,i) == 0, ~isnan(nanmask));
        
        if sum(nuclearMask(:,:,i), "all") == 0
            thres_val = nan;
            bg = nan;
        else
            thres_val = quantile(rm(mask_current), 0.5);
            bg = quantile(rm(mask_bg), 0.5);
        end
        
        ts(i) = thres_val;
        bgs(i) = bg;
    end
    
    qtl = zeros(1,N);
    npts = zeros(1,N);
    for i = 1:N
        rm = rawImages{i};
        rm = rm(nuclearMask(:,:,i) > 0);
        if isnan(ts(i))
            qtl(i) = nan;
            npts(i) = nan;
        else
            qtl(i) = sum(rm(rm > ts(i)) - bgs(i));
            npts(i) = sum(rm > ts(i));
        end
    end
    
    hoechstsum = qtl;
    smooththres = ts;
    smoothbg = bgs;
    
    save(output, "hoechstsum", "npts", "smooththres", "smoothbg", "idxToFrame");
end
