function watershedMask = getMaskForWatershed(preProcessedBrightField, enhancedBordersImage, segmentationParameters)
    %   getMaskForWatershed Create a mask with minima at droplet locations (and background locations) to perform watershed on
    %
    %   watershedMask = getMaskForWatershed(preProcessedBrightField, enhancedBordersImage, segmentationParameters)
    %
    %   Inputs:
    %       preProcessedBrightField (MxN) double - Pre-processed bright field image via `preProcessBrightField`
    %       enhancedBordersImage (MxN) logical - Image with enhanced borders via `enhanceDropletBorders`
    %       segmentationParameters (1x1) segmentation.Parameters - Parameters for segmentation
    %
    %   Outputs:
    %       watershedMask (MxN) logical - Mask with minima at droplet locations (and background locations) to perform watershed on
    arguments
        preProcessedBrightField (:,:) double
        enhancedBordersImage (:,:) logical
        segmentationParameters (1,1) segmentation.Parameters
    end
    dropletMarkers = segmentation.getDropletMarkers(enhancedBordersImage, segmentationParameters);
    maxCorrectedBrightField = segmentation.getMaxCorrectedBrightField(preProcessedBrightField);
    backgroundMarkers = segmentation.getBackgroundMarkers(maxCorrectedBrightField, dropletMarkers);
    watershedMask = imimposemin(maxCorrectedBrightField, dropletMarkers | backgroundMarkers);
end