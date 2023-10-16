function segmentationResult = segmentBrightFieldImage(brightFieldImage, segmentationParameters)
    % TODO: Add description
    arguments (Input)
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters segmentation.Parameters
    end

    arguments (Output)
        segmentationResult segmentation.Result
    end
    % Pre-process image to make it suitable for thresholding and segmentation
    preProcessedImage = segmentation.preProcessBrightFieldImage(brightFieldImage, segmentationParameters);
    % Convert image into binary mask
    level = graythresh(preProcessedImage);
    binaryMask = imbinarize(preProcessedImage, level);
    % Apply watershed algorithm (see https://www.mathworks.com/help/images/ref/watershed.html)
    distanceTransform = -bwdist(~binaryMask);
    distanceTransform(~binaryMask) = -Inf; % Set background pixels to -Inf so they become a catchment basin
    labeledImage = watershed(distanceTransform);
    % Calculate the properties of each region (see https://www.mathworks.com/help/images/ref/regionprops.html)
    regionProperties = regionprops(labeledImage, brightFieldImage, ...
                                   'Centroid','Area', ...
                                   'PixelIdxList','PixelValues', ...
                                   'Eccentricity','EulerNumber');
    % Filter segmentation result
    tooSmall = [regionProperties.Area] < segmentationParameters.minSegmentedArea;
    tooLarge = [regionProperties.Area] > segmentationParameters.maxSegmentedArea;
    notRound = [regionProperties.Eccentricity] > segmentationParameters.maxEccentricity;
    tooManyHoles = [regionProperties.EulerNumber] < 1 - segmentationParameters.maxNumberOfHoles;
    filter = tooSmall | tooLarge | notRound | tooManyHoles;
    % Remove regions that do not pass the filter from the labeled image and region properties
    filteredLabeledImage = labeledImage;
    filteredLabeledImage(vertcat(regionProperties(filter).PixelIdxList)) = 0;
    filteredRegionProperties = regionProperties(~filter); 
    % Create segmentation result
    segmentationResult = segmentation.Result(filteredLabeledImage, filteredRegionProperties);
end