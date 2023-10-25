function preProcessed = preProcessBrightFieldImage(brightFieldImage, segmentationParameters)
    % preProcessBrightFieldImage Corrects the bright field image for uneven illumination and normalizes it.
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters (1,1) segmentation.Parameters
    end
    corrected = imflatfield(brightFieldImage, segmentationParameters.illuminationCorrectionFilterSize);
    corrected = double(corrected);
    imgMin = min(corrected(:));
    imgMax = max(corrected(:));
    preProcessed = (corrected - imgMin) / (imgMax - imgMin);
end