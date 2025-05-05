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
%   labels          - (1,1) string. File pattern for label images
%   nucSegOutputMat - (1,1) string. Filename of the .mat file containing nuclearMask and idxToFrame.
%   output          - (1,1) string. Output filename to save the results.
%
% Outputs:
%   dnaMask    - 3-D matrix for dna mask result
%   hoechstArea- 1×N vector of total number of positive Hoechst intensity pixels
%   idxToFrame - 1×N vector mapping image index to frame number.
%   fs         - string of reffered image file path
%   NucDNAMask - 3-D matrix for nuclear mask result. When nuclear area is
%                   empty in orignal result, DNA mask is transferred as
%                   nuclear tentative mask.%   
%
% Example:
%   [hoechstArea, idxToFrame] = cropDNAMask(files, labels, nucSegOutputMat, output)

    arguments
        files (1,1) string
        labels (1,1) string
        nucSegOutputMat (1,1) string
        output (1,1) string
    end

    %% parameters
    bwareaopenthresh = 15;

    %% fileIO
    % Load nuclear segmentation results. All droplets have nuclear mask.
    nucData = load(nucSegOutputMat);
    nucMask = nucData.nuclearMask;

    % output allocation
    dnaMask = zeros(size(nucMask));
    NucDNAMask = zeros(size(nucMask));

    %% DNA area segmentation
    % Nucler ojbect check: if no nucleus in entire frames, detection process
    % is passed. Save empty variables.
    if max(nucData.nuclearArea) ~= 0

        % collect DAPI images information
        fs = dir(files);
        N = length(fs);
        if N == 0
            error("No files found matching pattern: %s", files);
        end
        imgroot = fs(1).folder;

        % collect Droplet labels information
        ls = dir(labels);
        N2 = length(ls);
        if N2 == 0
            error("No label images: %s", labels);
        end
        lblroot = ls(1).folder;

        % sanity check
        if size(nucMask,3) ~= N
            error("Frame count mismatch: nuclearMask and image files different.");
        end

        % separate image to image stack
        rawImages = cell(1, N);
        labelImage = cell(1, N);
        for i = 1:N
            rawImages{i} = imread(fullfile(imgroot, fs(i).name));
            nPixels = size(rawImages{i}, 1);
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
        
        %% skewness values of intensity histogram for each frame
        % The skewness of the histogram changes with the degree of DNA condensation.
        skewVec = zeros(N,1);
        for t = 1:N
            I = double(maskedImages(:,:,t));

            se = strel('disk', 20);
            I_tophat = imtophat(I, se);
            
            mask = I_tophat > 0;
            pix = I(mask);

            skewVec(t) = skewness(double(pix));
        end
        smoothSkew = movmean(skewVec, [1 1]);


        %% Process each frame to segment nucleus.
        for i = 1:N
            currentImage = double(maskedImages(:,:,i));
            % background cleaning 
            se = strel('disk', 20);
            I_tophat = imtophat(currentImage, se);
            I_tophat_NaN = I_tophat;
            I_tophat_NaN(I_tophat==0) = NaN;

            % pixel information of NaN
            pix = I_tophat_NaN(~isnan(I_tophat_NaN));

            if smoothSkew(i) > median(skewVec) || smoothSkew(i) > 4
                % detect outlier pixels in very skewed intensity histogram
                tf = isoutlier(pix,'quartiles');
                outlierMask = false(size(I));
                outlierMask(~isnan(I_tophat_NaN)) = tf;
                BW = outlierMask;

            else
                % percentile scaling in modarate skewed intensity histogam
                p1 = prctile(I_tophat_NaN(:), 1);
                p90= prctile(I_tophat_NaN(:), 90);
                %p1_series = [p1_series, p1];
                %p90_series = [p90_series, p90];
                I_norm = imadjust(I_tophat/max(I_tophat(:)), [p1 p90]/max(I_tophat(:)), [0 1]);
                BW = imbinarize(I_norm);
                
            end

            % data cleaning
            % 1) labeling
            CC = bwconncomp(BW);

            % 2) calculation of detected object stats
            stats = regionprops(CC, 'Area', 'Eccentricity', 'Solidity', 'PixelIdxList');

            % 3) filtering parameters
            eccThresh = 0.95;
            solThresh = 0.5;

            
            toRemove = false(size(stats));
            for h = 1:numel(stats)
                badShape = stats(h).Eccentricity > eccThresh || stats(h).Solidity < solThresh;

                % if detected area is on the periferi
                regionMask = false(size(BW));
                regionMask(stats(h).PixelIdxList) = true;
                % 1-pixel dilation
                regionDilated = imdilate(regionMask, strel('disk',1));
                % if dialated area is overllaped with outer region of dropletMask -> true
                touchesEdge = any( regionDilated(:) & ~I(:) );

                if badShape && touchesEdge
                    toRemove(h) = true;
                end

                % detected reagion is too large ignore the area (0.2 based on
                % maximum area nucleus is ~10% of droplet area)
                d_area = sum(I(:)>0);
                if stats(h).Area > d_area * 0.33
                    toRemove(h) = true;
                end
            end
            
            for h = find(toRemove)'
                BW(stats(h).PixelIdxList) = false;
            end

            BW = imfill(BW, 'holes');
            BW = bwareaopen(BW, bwareaopenthresh);

            dnaMask(:,:,i) = BW;


            %% compare nuclear amsk and DNA mask upgrade the accuracy of nuclear segmentation result
            tmpNucMask = double(nucMask(:,:,i));
            binaryDnaMask = dnaMask(:,:,i);
            diffMask = tmpNucMask - binaryDnaMask;
            if any(diffMask(:) < 0) && ~any(tmpNucMask, 'all')
                %fprintf(' - Nuc Area updated based on Hoechst data -');
                NucDNAMask(:,:,i) = tmpNucMask + binaryDnaMask;
            end

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
        name = sprintf("DnaMask_Overlay_%03d.tif",i);
        
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
    %save(output, "dnaMask","hoechstArea", "npts", "smooththres", "smoothbg", "idxToFrame");
    %}
    
    save(output, "dnaMask","hoechstArea", "idxToFrame","fs", "NucDNAMask");
end
