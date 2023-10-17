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
    binaryMask = segmentation.binarizeBrightFieldImage(preProcessedImage, segmentationParameters);
    % Apply watershed algorithm (see https://www.mathworks.com/help/images/ref/watershed.html)
    distanceTransform = -bwdist(~binaryMask);
    distanceTransform(~binaryMask) = -Inf; % Set background pixels to -Inf so they become a catchment basin
    labeledImage = watershed(distanceTransform);
    % Filter segmentation result
    [filteredLabeledImage, filteredRegionProperties] = segmentation.filterRegions(brightFieldImage, ...
                                                                                  labeledImage, ... 
                                                                                  segmentationParameters);
    % Create segmentation result
    segmentationResult = segmentation.Result(filteredLabeledImage, filteredRegionProperties);
end