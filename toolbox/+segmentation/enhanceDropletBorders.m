function enhancedImage = enhanceDropletBorders(preProcessedBrightFieldImage, segmentationParameters)
    %   enhanceDropletBorders Creates a binary mask that enhances droplet borders via gaussian, hessian, and region filtering
    %
    %   enhancedImage = enhanceDropletBorders(preProcessedBrightFieldImage, segmentationParameters)
    %
    %   Inputs:
    %       preProcessedBrightFieldImage (MxN) double - Bright field image that has been pre-processed with `preProcessBrightFieldImage`
    %       segmentationParameters (1,1) segmentation.Parameters - Parameters for segmentation
    %
    %   Output:
    %       enhancedImage (MxN) logical - Binary mask of droplet borders. Borders are set to 1 and background is set to 0.
    arguments
        preProcessedBrightFieldImage (:,:) double
        segmentationParameters (1,1) segmentation.Parameters
    end
    % Difference of gaussian filters
    largeFilter = imgaussfilt(preProcessedBrightFieldImage, segmentationParameters.gaussianFilterSizeLarge);
    smallFilter = imgaussfilt(preProcessedBrightFieldImage, segmentationParameters.gaussianFilterSizeSmall);
    dropletBorders = largeFilter - smallFilter;
    dropletBorders(dropletBorders < 0) = 0;
    % Hessian filtering
    [Dxx, Dxy, Dyy] = segmentation.hessian2D(dropletBorders, segmentationParameters);
    lambda = segmentation.largestHessianEigenvalue(Dxx, Dxy, Dyy);
    reconstruction = (ones(size(dropletBorders)) - exp(-(lambda.^2)/segmentationParameters.hessianEigenvalueScale));
    hessianFiltered = reconstruction .* (lambda < 0) .* dropletBorders;
    % Thresholding
    threshold = graythresh(hessianFiltered) * segmentationParameters.thresholdFactor;
    binaryMask = hessianFiltered > threshold;
    % Region filtering
    p = regionprops(logical(binaryMask), 'Area', 'Eccentricity', 'PixelIdxList');
    tooSmall = [p.Area] < (segmentationParameters.borderSmallAreaThreshold * segmentationParameters.resolution ^ 2);
    mediumSize = [p.Area] < segmentationParameters.borderMediumAreaThreshold * segmentationParameters.resolution ^ 2;
    tooSymmetric = [p.Eccentricity] < segmentationParameters.borderEccentricityCutoff;
    toRemove = tooSmall | (mediumSize & tooSymmetric); 
    labels = bwlabel(binaryMask);
    filteredLabels = labels;
    filteredLabels(vertcat(p(toRemove).PixelIdxList)) = 0;
    selectedRegions = logical(filteredLabels);
    enhancedImage = dropletBorders & selectedRegions;
end