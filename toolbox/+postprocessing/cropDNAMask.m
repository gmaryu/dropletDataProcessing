function [hoechstArea, idxToFrame] = cropDNAMask(files, labels, nucSegOutputMat, output)
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
        labels (1,1) string
        nucSegOutputMat (1,1) string
        output (1,1) string
    end

    radiusMargin = 0.15 + 0.05;
    gfilterPixels = 1;
    outlierThresFactor = 3;
    refQuantileLower = 0.5;
    refQuantileUpper = 0.95;
    outlierThres = 0.01;
    smallCcThres = 0.0005;
    gfitMaxStdFactor = 10;
    intensityThresFactor = 1.5;

    % Load nuclear segmentation results.
    nucData = load(nucSegOutputMat);
    nucMask = nucData.nuclearMask;

    % output allocation
    dnaMask = zeros(size(nucMask));
    NucDNAMask = zeros(size(nucMask));

    % Nucler ojbect check if no nucleus in entire frames, detection process
    % is passed. Save empty variables.

    if max(nucData.nuclearArea) ~= 0

        % DAPI images
        fs = dir(files);
        N = length(fs);
        if N == 0
            error("No files found matching pattern: %s", files);
        end
        imgroot = fs(1).folder;

        % Droplet labels
        ls = dir(labels);
        N2 = length(ls);
        if N2 == 0
            error("No label images: %s", labels);
        end
        lblroot = ls(1).folder;

        if size(nucMask,3) ~= N
            error("Frame count mismatch: nuclearMask and image files differ.");
        end

        rawImages = cell(1, N);
        bgs = zeros(1,N);
        ts = zeros(1,N);
        for i = 1:N
            rawImages{i} = imread(fullfile(imgroot, fs(i).name));
            nPixels = size(rawImages{i}, 1);
            tmplbl = imread(fullfile(lblroot, ls(i).name));
            %{
        if i == 1
            nPixels = size(rawImages{1}, 1);
            [xx, yy] = meshgrid(1:nPixels, 1:nPixels);
            radius = nPixels / 2 / (1 + radiusMargin);
            mask = uint16(hypot(xx - (nPixels+1)/2, yy - (nPixels+1)/2) < radius);
            nanmask = double(mask);
            nanmask(nanmask==0) = NaN;
        end
            %}

            % label image stack
            tmplabel = imbinarize(imread(fullfile(lblroot, ls(i).name)));
            labelImage{i} = tmplabel;
        end

        try
            % multi-cell to 3-dimensional matrix
            rawImagesCat = cat(3, rawImages{:});
            labelImageCat = cat(3, labelImage{:});
        catch
            %nuclearArea = [];
            %idxToFrame = [];
            %return;
        end

        maskedImages = uint16(labelImageCat) .* rawImagesCat;
        

        %% Process each frame to segment nucleus.
        for i = 1:N
            currentImage = double(maskedImages(:,:,i));
            nanmask = labelImage{i};
            mask = nanmask;
            
            gfilt = imgaussfilt(currentImage .* nanmask, gfilterPixels);
            gfilt(gfilt==0)=NaN;
            [mint, idx] = max(gfilt, [], "all", "linear");
            [ix, iy] = ind2sub(size(gfilt), idx);
            ref = median(gfilt(gfilt > quantile(gfilt(:), refQuantileLower) & gfilt < quantile(gfilt(:), refQuantileUpper)));
            t_val = outlierThresFactor * -1/(sqrt(2)*erfcinv(3/2)) * median(abs(gfilt(:)-ref), "omitnan");
            outliers = imgaussfilt(double(gfilt - ref > t_val), gfilterPixels);
            outliers(outliers < 0.3) =0;
            
            
            currentImage2 = currentImage .* nanmask;
            ref2 = median(currentImage2(currentImage2 > quantile(currentImage2(:), refQuantileLower) & currentImage2 < quantile(currentImage2(:), refQuantileUpper)));
            t_val2 = outlierThresFactor * -1/(sqrt(2)*erfcinv(3/2)) * median(abs(currentImage2(:)-ref2), "omitnan");
            outliers2 = double(currentImage2 - ref2 > t_val2);
            dnaMask(:,:,i) = outliers2;

            %% too large autofluorescence junk
            dropletArea = sum(double(nanmask(:)));
            binaryDnaMask = double(outliers2 > 0);
            brightArea = sum(binaryDnaMask(:));
            if brightArea/dropletArea > 0.2
                fprintf(' - There might be a large junk. Skip this droplet -');
                hoechstArea = NaN*ones(1,N);
                idxToFrame = NaN*ones(1,N);
                return
            end

            %% compare nuclear amsk and DNA mask upgrade the accuracy of nuclear segmentation result
            tmpNucMask = double(nucMask(:,:,i));
            diffMask = tmpNucMask - binaryDnaMask;
            if any(diffMask(:) < 0)
                %fprintf(' - Nuc Area updated based on Hoechst data -');
                NucDNAMask(:,:,i) = tmpNucMask + binaryDnaMask;

            end

            %{
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
                dnaMask(:,:,i) = bwareaopen(standardized > intensityThresFactor, floor(nPixels^2 * smallCcThres));
            else
                dnaMask(:,:,i) = 0;
            end
            %}
        end
        hoechstsum =[];
        npts = [];
        smooththres = [];
        smoothbg = [];
    else
        fprintf(' - No nuclei detection in FP channel -');
        hoechstsum =[];
        npts = [];
        smooththres = [];
        smoothbg = [];
        hoechstArea = NaN*ones(1,size(nucMask,3));
        idxToFrame = NaN*ones(1,size(nucMask,3));
        return
    end

    positiveHoechst = dnaMask;
    positiveHoechst(dnaMask > 0) = 1;
    hoechstArea = sum(reshape(positiveHoechst, nPixels*nPixels, N), 1);


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

    %% (Optional) Save overlay images.
    for i = 1:N
        overlay = imoverlay(rawImagesCat(:,:,i), bwperim(dnaMask(:,:,i)), [0, 1, 0]);
        name = strrep(fs(i).name, ".tif", "_segmented_DNA.tif");
        % Uncomment the line below to save overlay images.
        %imwrite(overlay, fullfile(imgroot, name));
    end
    %{
        %% using segmented nucleus data
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
    
    
    %}
    %save(output, "dnaMask","hoechstArea", "npts", "smooththres", "smoothbg", "idxToFrame");
    save(output, "dnaMask","hoechstArea", "idxToFrame","fs", "NucDNAMask");
end
