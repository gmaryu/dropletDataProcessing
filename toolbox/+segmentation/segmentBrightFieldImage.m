function segmentationResult = segmentBrightFieldImage(brightFieldImage, segmentationParameters)
    %   segmentBrightFieldImage Segments a bright field image into droplets via the watershed algorithm
    %
    %   segmentationResult = segmentBrightFieldImage(brightFieldImage, segmentationParameters) 
    %
    %   Inputs:
    %       brightFieldImage (MxN) numeric - Bright field image to segment
    %       segmentationParameters (1,1) segmentation.Parameters - Parameters for the segmentation
    %
    %   Output:
    %       segmentationResult (1,1) segmentation.Result - Result object containing the labeled image and region properties
    arguments (Input)
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters (1,1) segmentation.Parameters
    end

    arguments (Output)
        segmentationResult segmentation.Result
    end
    preProcessed = segmentation.preProcessBrightFieldImage(brightFieldImage, segmentationParameters);
    enhancedBordersImage = segmentation.enhanceDropletBorders(preProcessed, segmentationParameters);
    watershedMask = segmentation.getMaskForWatershed(preProcessed, enhancedBordersImage, segmentationParameters);
    labeledImage = segmentation.labelObjects(watershedMask, segmentationParameters);
    [filteredLabeledImage, filteredRegionProperties] = segmentation.filterLabels(labeledImage, segmentationParameters);
    segmentationResult = segmentation.Result(filteredLabeledImage, filteredRegionProperties);
end