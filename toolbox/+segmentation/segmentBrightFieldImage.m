function segmentationResult = segmentBrightFieldImage(brightFieldImage, segmentationParameters)
    % segmentBrightFieldImage Segments a bright field image into droplets via the watershed algorithm
    arguments (Input)
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters (1,1) segmentation.Parameters
    end

    arguments (Output)
        segmentationResult segmentation.Result
    end
    preProcessed = segmentation.preProcessBrightFieldImage(brightFieldImage, segmentationParameters);
    enhancedBordersImage = segmentation.enhanceDropletBorders(preProcessed, segmentationParameters);
    watershedMask = segmentation.getMaskForWatershed(brightFieldImage, enhancedBordersImage, segmentationParameters);
    labeledImage = segmentation.labelObjects(watershedMask, segmentationParameters);
    [filteredLabeledImage, filteredRegionProperties] = segmentation.filterLabels(labeledImage, segmentationParameters);
    segmentationResult = segmentation.Result(filteredLabeledImage, filteredRegionProperties);
end