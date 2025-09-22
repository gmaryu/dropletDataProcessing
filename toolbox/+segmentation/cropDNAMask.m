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
    opts = struct;
    opts.tophatR     = 8;
    opts.stdThresh   = 0.004;
    opts.ratioThresh = 2;
    opts.minArea     = 5;
    opts.maxSpots    = 5;

    R      = getfield(opts, 'tophatR');
    stdT   = getfield(opts, 'stdThresh');
    ratT   = getfield(opts, 'ratioThresh');
    minA   = getfield(opts, 'minArea');
    maxN   = getfield(opts, 'maxSpots');
    %dbg    = getfield(opts, 'showDebug');
    
    bwareaopenthresh = minA;

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

        % sanity check
        if size(nucMask,3) ~= N
            error("Frame count mismatch: nuclearMask and image files different.");
        end

        % separate image to image stack
        rawImages = cell(1, N);
        for i = 1:N
            rawImages{i} = imread(fullfile(imgroot, fs(i).name));
        end

        try
            % multi-cell to 3-dimensional matrix
            rawImagesCat = cat(3, rawImages{:});
        catch
            % padding with 0
            T = numel(rawImages);
            H = cellfun(@(a) size(a,1), rawImages);
            W = cellfun(@(a) size(a,2), rawImages);
            Hmax = max(H);  Wmax = max(W);

            rawImagesCat = zeros(Hmax, Wmax, T, 'like', rawImages{1});   % 0パディング

            for t = 1:T
                h = H(t); w = W(t);
                r0 = floor((Hmax - h)/2);
                c0 = floor((Wmax - w)/2);
                % r0 = 0;
                % c0 = 0;
                rawImagesCat( (1:h)+r0, (1:w)+c0, t ) = rawImages{t};
            end
        end

        labelImageCat = nucData.labelImageCat;
        maskedImages = uint16(labelImageCat) .* rawImagesCat;

        for i=1:N
            % ----------- 1. robust normalisation -------------------------------------
            img  = double(maskedImages(:,:,i));
            medI = median(img(:));
            madI = mad(img(:),1);                     % robust σ̂
            imgZ = (img - medI) / max(madI, eps);     % (≈ z-score)
            imgN = rescale(imgZ, 0, 1);               % to [0,1] for thresholding

            % Global contrast gate (skip unmistakably empty frames fast)
            %
            if madI/ max(img(:)) < stdT
                BW = false(size(img)); 
                continue;
            end
            if prctile(imgN(:),99) / max(prctile(imgN(:),33),eps) < ratT
                BW = false(size(img));  
                continue;
            end

            %}
            % ----------- 2. background suppression (white top-hat) -------------------
            se   = strel('disk', R);
            toph = imtophat(imgN, se);        % emphasise objects of ~size R

            % ----------- 3. Otsu threshold & area filtering --------------------------
            level  = graythresh(toph);
            %level  = multithresh(toph,3);
            BWc    = imbinarize(toph, level);
            BW     = bwareaopen(BWc, minA);   % remove crumbs < minA pixels
            
            if isempty(BW)
                %level  = graythresh(toph);
                level  = multithresh(toph,3);
                BWc    = imbinarize(toph, level(1));
                BW     = bwareaopen(BWc, minA);   % remove crumbs < minA pixels
            end

            % ----------- 4. enforce “1–3 objects” rule -------------------------------
            CC = bwconncomp(BW, 8);
            if CC.NumObjects==0 || CC.NumObjects>maxN
                BW = false(size(img));        % treat as “no valid bright object”
            else
                % (Optional) keep only the largest ‘maxN’ components
                sizes  = cellfun(@numel, CC.PixelIdxList);
                [~,ix] = sort(sizes,'descend');
                keep   = ix(1:min(maxN, numel(ix)));
                BW = false(size(img));
                BW(cat(1,CC.PixelIdxList{keep})) = true;
            end

            % 2) calculation of detected object stats
            stats = regionprops(BW, 'Area', 'Eccentricity', 'Solidity', 'PixelIdxList');

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
                touchesEdge = any( regionDilated(:) & ~img(:) );

                if badShape && touchesEdge
                    toRemove(h) = true;
                end

                % detected reagion is too large ignore the area (0.2 based on
                % maximum area nucleus is ~10% of droplet area)
                d_area = sum(img(:)>0);
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
            tmpNucMask = nucMask(:,:,i);
            binaryDnaMask = dnaMask(:,:,i);
            NucDNAMask(:,:,i) = tmpNucMask | binaryDnaMask;

        end
        positiveHoechst = logical(dnaMask);
        hoechstArea = squeeze(sum(positiveHoechst,[1 2]));
        fprintf(" -- DNA quantification completed.\n");

    else
        fprintf(' -- No nuclei detection in FP channel -- END.\n');
        hoechstArea = NaN*ones(1,size(nucMask,3));
        idxToFrame = NaN*ones(1,size(nucMask,3));
        return
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

    %% (Optional) Save overlay images.
    %{
    for i = 1:N
        overlay = imoverlay(rawImagesCat(:,:,i), bwperim(dnaMask(:,:,i)), [0, 1, 0]);
        name = sprintf("DnaMask_Overlay_%03d.tif",i);
        
        % Uncomment the line below to save overlay images.
        imwrite(overlay, fullfile(imgroot, name));
    end
    %}
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
