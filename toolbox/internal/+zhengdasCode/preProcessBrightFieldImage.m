function processedImage = preProcessBrightFieldImage(brightFieldImage, segmentationParameters)
    % Pre-processing steps applied in the original algorithm developed by Zhengda Li
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters zhengdasCode.Parameters
    end
    
    brightFieldImage = double(brightFieldImage);
    % Normalize
    imgMin = min(brightFieldImage(:));
    imgMax = max(brightFieldImage(:));
    processedImage = (brightFieldImage - imgMin) / (imgMax - imgMin);
end