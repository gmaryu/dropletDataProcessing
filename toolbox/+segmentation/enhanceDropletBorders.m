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
    hessianFiltered = dropletBorders;
    % Thresholding
    threshold = graythresh(hessianFiltered);
    binaryMask = hessianFiltered > threshold;
    % Region filtering
    p = regionprops(logical(binaryMask), 'Area', 'Solidity', 'PixelIdxList');
    tooSmall = [p.Area] < (segmentationParameters.borderAreaThreshold * segmentationParameters.resolution ^ 2);
    tooSolid = [p.Solidity] > (segmentationParameters.borderSolidityThreshold);
    toRemove = tooSmall & tooSolid; 
    labels = bwlabel(binaryMask);
    filteredLabels = labels;
    filteredLabels(vertcat(p(toRemove).PixelIdxList)) = 0;
    selectedRegions = logical(filteredLabels);
    enhancedImage = dropletBorders & selectedRegions;
end