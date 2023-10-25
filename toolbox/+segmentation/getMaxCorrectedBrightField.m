function maxCorrectedBrightField = getMaxCorrectedBrightField(brightFieldImage)
    % getMaxCorrectedBrightField Calculate the difference between the maximum value of the image and the image itself
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
    end
    image = double(brightFieldImage);
    maxCorrectedBrightField = max(image(:)) - image;
end