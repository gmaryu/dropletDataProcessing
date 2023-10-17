function binaryMask = binarizeBrightFieldImage(brightFieldImage, segmentationParameters)
    % TODO: Add description
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        segmentationParameters segmentation.Parameters
    end

    level = graythresh(brightFieldImage);
    binaryMask = imbinarize(brightFieldImage, level);
    % Enhance round shapes via morphological opening with a disk
    binaryMask = imopen(binaryMask, ...
                        strel('disk', round(segmentationParameters.morphologicalOpeningRadius)));
    % Remove small objects and holes
    binaryMask = bwareaopen(binaryMask, segmentationParameters.minSegmentedArea);
    binaryMask = imfill(binaryMask, "holes");
end