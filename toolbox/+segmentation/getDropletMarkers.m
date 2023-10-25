function dropletMarkers = getDropletMarkers(enhancedBordersImage, segmentationParameters)
    % getDropletMarkers Calculate droplet center markers using the distance transform of the enhanced borders image
    arguments
        enhancedBordersImage (:,:) logical
        segmentationParameters (1,1) segmentation.Parameters
    end
    distanceTransform = bwdist(enhancedBordersImage);
    dropletMarkers = imextendedmax(distanceTransform, segmentationParameters.HMaximaSuppressionThreshold);
end