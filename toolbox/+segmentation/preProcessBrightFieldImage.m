function processedImage = preProcessBrightFieldImage(brightFieldImage, segmentationParameters)
    % TODO: Add description
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters segmentation.Parameters
    end

    % TODO: Implement pre-processing steps
    gaussFilteredImage = imgaussfilt(brightFieldImage, segmentationParameters.gaussianFilterSigma);
    processedImage = gaussFilteredImage;
end