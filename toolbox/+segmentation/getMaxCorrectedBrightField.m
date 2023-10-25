function maxCorrectedBrightField = getMaxCorrectedBrightField(preProcessedBrightField)
    %   getMaxCorrectedBrightField Calculate the difference between the maximum value of the image and the image itself
    %
    %   maxCorrectedBrightField = getMaxCorrectedBrightField(preProcessedBrightField)
    %   
    %   Input:
    %       preProcessedBrightField - The image to be processed
    %
    %   Output:
    %       maxCorrectedBrightField - Same image where all the values have been subtracted from the maximum value of the image
    arguments
        preProcessedBrightField (:,:) double
    end
    maxCorrectedBrightField = max(preProcessedBrightField(:)) - preProcessedBrightField;
end