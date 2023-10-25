function labeledImage = labelObjects(watershedMask, segmentationParameters)
    % getDropletMarkers Label all objects in `watershedMask` using watershed
    arguments
        watershedMask (:,:) {mustBeNumeric} % Output from getWatershedMask
        segmentationParameters segmentation.Parameters
    end
    labels = watershed(watershedMask);
    binaryLabels = logical(labels);
    diskSize = segmentationParameters.openingDiskSize * segmentationParameters.resolution;
    disk = strel('disk', diskSize);
    binaryLabels = imopen(binaryLabels, disk); % Make labels more circular
    labeledImage = double(labels) .* double(binaryLabels);
end