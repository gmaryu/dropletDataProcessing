function watershedMask = getMaskForWatershed(brightFieldImage, enhancedBordersImage, segmentationParameters)
    % getMaskForWatershed Create a mask with minima at droplet locations (and background) to perform watershed on
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        enhancedBordersImage (:,:) logical
        segmentationParameters (1,1) segmentation.Parameters
    end
    dropletMarkers = segmentation.getDropletMarkers(enhancedBordersImage, segmentationParameters);
    maxCorrectedBrightField = segmentation.getMaxCorrectedBrightField(brightFieldImage);
    backgroundMarkers = segmentation.getBackgroundMarkers(maxCorrectedBrightField, dropletMarkers);
    watershedMask = imimposemin(maxCorrectedBrightField, dropletMarkers | backgroundMarkers);
end