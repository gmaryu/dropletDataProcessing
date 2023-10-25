function dropletMarkers = getDropletMarkers(enhancedBordersImage, segmentationParameters)
    %   getDropletMarkers Calculate droplet center markers using the distance transform of the enhanced borders image
    %
    %   dropletMarkers = getDropletMarkers(enhancedBordersImage, segmentationParameters)
    % 
    %   Inputs:
    %       enhancedBordersImage (MxN) logical - Image of enhanced borders obtained via `enhanceDropletBorders`
    %       segmentationParameters (1x1) segmentation.Parameters - Parameters for the segmentation
    %
    %   Output:
    %       dropletMarkers (MxN) logical - Image with important regions marked as 1
    arguments
        enhancedBordersImage (:,:) logical
        segmentationParameters (1,1) segmentation.Parameters
    end
    distanceTransform = bwdist(enhancedBordersImage);
    dropletMarkers = imextendedmax(distanceTransform, segmentationParameters.HMaximaSuppressionThreshold);
end