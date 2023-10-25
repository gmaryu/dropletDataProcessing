function preProcessed = preProcessBrightFieldImage(brightFieldImage, segmentationParameters)
    %   preProcessBrightFieldImage Corrects the bright field image for uneven illumination and normalizes it.
    %
    %   preProcessed = preProcessBrightFieldImage(brightFieldImage, segmentationParameters)
    %
    %   Inputs:
    %       brightFieldImage (MxN) numeric - The bright field image to be corrected.
    %       segmentationParameters (1x1) segmentation.Parameters - The parameters for the segmentation.
    %
    %   Output:
    %       preProcessed (MxN) double - The corrected and normalized bright field image.
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