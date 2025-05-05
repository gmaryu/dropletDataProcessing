function [nuclearArea, idxToFrame] = cropBrightChunk(files, labels, output)
% cropBrightChunk Segments nuclei from a series of bright field images.
%
%   [nuclearArea, idxToFrame] = cropBrightChunk(files, output)
%
% This function reads a set of TIFF images specified by the file pattern in 'files', applies a 
% Gaussian filter and outlier detection to identify the nucleus, fits a two‐dimensional Gaussian 
% model, and produces a binary mask of the nucleus. It returns the nuclear area for each frame 
% and a mapping from image index to frame number. Results are saved to the specified output file.
%
% Inputs:
%   imgfiles  - (1,1) string. File pattern for the images, e.g., "Pos6_RFP-T_???.tif".
%   labelfiles    - (1,1) string. File pattern for the labels, e.g., "Pos0_label_???.tif"
%   output - (1,1) string. Output .mat filename to save the results.
%
% Outputs:
%   nuclearArea - 1×N vector of nuclear area (in pixels) for each frame.
%   idxToFrame  - 1×N vector mapping each image to its frame index.
%
% Example:
%   [area, idxFrame] = postprocessing.cropBrightChunk("test/droplet_000/Pos6_RFP-T_???.tif", "Pos6_000.mat");
    
    arguments
        files (1,1) string
        labels (1,1) string
        output (1,1) string
    end

    %% parameters for nuclear segmentation
    gfilterPixels = 1;
    outlierThresFactor = 1.5;
    refQuantileLower = 0.25;
    refQuantileUpper = 0.75;
    outlierThres = 0.01;
    smallCcThres = 0.0005;
    gfitMaxStdFactor = 10;
    intensityThresFactor = 1.5;
    
    %% fileIO and image2stack
    fs = dir(files);
    N = length(fs);
    if N == 0
        error("No files found matching pattern: %s", files);
    end

    ls = dir(labels);
    N2 = length(ls);
    if N2 == 0
        error("No label images: %s", labels);
    end

    imgroot = fs(1).folder;
    lblroot = ls(1).folder;
    rawImages = cell(1, N);
    for i = 1:N
        % raw image stack
        rawImages{i} = imread(fullfile(imgroot, fs(i).name));
        nPixels = size(rawImages{i}, 1);

        % label image stack
        tmplabel = imbinarize(imread(fullfile(lblroot, ls(i).name)));
        labelImage{i} = tmplabel;
        %{
        % artificial droplet mask 
        if i == 1
            [xx, yy] = meshgrid(1:nPixels, 1:nPixels);
            radius = nPixels / 2 / (1 + radiusMargin);
            mask = uint16(hypot(xx - (nPixels+1)/2, yy - (nPixels+1)/2) < radius);
            nanmask = double(mask);
            nanmask(nanmask==0) = NaN;
        end
        %}
    end
    
    try
        % multi-cell to 3-dimensional matrix 
        rawImagesCat = cat(3, rawImages{:});
        labelImageCat = cat(3, labelImage{:});
    catch
        nuclearArea = [];
        idxToFrame = [];
        return;
    end
    
    % replace pixels of droplet outside to 0
    maskedImages = uint16(labelImageCat) .* rawImagesCat; % xyt stack
    % matrix allocation for out
    nuclearMask = double(maskedImages);
    
    %% Process each frame to segment nucleus.
    for i = 1:N
        currentImage = double(maskedImages(:,:,i));
        nanmask = labelImage{i};
        mask = nanmask;

        % gaussian fitting of to detect nuclear area (brightest area) and 
        % estimate cytoplasm area intensity as ref
        gfilt = imgaussfilt(currentImage .* nanmask, gfilterPixels);
        gfilt(gfilt==0)=NaN;
        [mint, idx] = max(gfilt, [], "all", "linear");  % find brightest spot
        [ix, iy] = ind2sub(size(gfilt), idx);           % position of brightest sport
        ref = median(gfilt(gfilt > quantile(gfilt(:), refQuantileLower) & gfilt < quantile(gfilt(:), refQuantileUpper))); % collect intensity profile lower and upper quantile (possibliy cytoplasm area)
        t_val = outlierThresFactor * -1/(sqrt(2)*erfcinv(3/2)) * median(abs(gfilt(:)-ref), "omitnan");
        outliers = imgaussfilt(double(gfilt - ref > t_val), gfilterPixels);
        
        if outliers(ix, iy) > outlierThres
            [my, mx] = meshgrid(1:size(mask,1), 1:size(mask,2));
            idx_fit = ~isnan(gfilt);
            g2d_residuals = @(p, x, y, z) p(1)*exp(-((x-ix).^2+(y-iy).^2)/p(2)) + p(3) - z;
            opts = optimset('Display','off');
            initialGuess = [mint - mean(gfilt(:), 'omitnan'), 1, mean(gfilt(:), 'omitnan')];
            lb = [0, 0, 0];
            ub = [inf, (size(gfilt,1)/gfitMaxStdFactor)^2, inf];
            gfit = lsqnonlin(@(p) g2d_residuals(p, mx(idx_fit), my(idx_fit), gfilt(idx_fit)), initialGuess, lb, ub, opts);
            mask_cyto = gfit(1)*exp(-((mx-ix).^2+(my-iy).^2)/gfit(2)) + gfit(3) < gfit(1)*exp(-3) + gfit(3);
            masked_cyto_images = gfilt .* double(mask_cyto);
            masked_cyto_images(masked_cyto_images == 0) = NaN;
            mu_val = mean(masked_cyto_images(:), 'omitnan');
            sigma_val = std(masked_cyto_images(:), 'omitnan');
            standardized = double(nanmask) .* ((gfilt - mu_val) / sigma_val);
            %standardized = double(mask) .* ((gfilt - mu_val) / sigma_val);
            nuclearMask(:,:,i) = bwareaopen(standardized > intensityThresFactor, floor(nPixels^2 * smallCcThres));
        else
            % no nucleus-like structure
            nuclearMask(:,:,i) = 0;
        end
    end
    
    nuclearArea = sum(reshape(nuclearMask, nPixels*nPixels, N), 1);
    
    %% (Optional) Save overlay images.
    for i = 1:N
        overlay = imoverlay(rawImagesCat(:,:,i), bwperim(nuclearMask(:,:,i)), [0, 1, 0]);
        name = sprintf("NucMask_Overlay_%03d.tif",i);

        % Uncomment the line below to save overlay images.
        %imwrite(overlay, fullfile(imgroot, name));
    end
    
    %% Determine frame indices from file names.
    [~, basename, ext] = fileparts(files);
    nameref = convertStringsToChars(basename + ext);
    idxToFrame = zeros(1,N);
    for i = 1:N
        currentName = fs(i).name;
        idx = "";
        for j = 1:strlength(currentName)
            if nameref(j)=='?'
                idx = idx + currentName(j);
            end
        end
        idxToFrame(i) = str2double(idx);
    end
    
    save(output, "nuclearMask", "nuclearArea", "idxToFrame", "fs");
end
