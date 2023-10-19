function segmentationCentroids = getSegmentationCentroids(labeledImage, regionProperties)
    % TODO: Add description
    arguments
        labeledImage (:,:) {mustBeNumeric}
        regionProperties (:,:) struct
    end
    % Get the label of each region
    regionLabels = NaN(1,length(regionProperties));
    for i = 1:length(regionProperties)
        pixel = regionProperties(i).PixelIdxList(1,1);
        label = labeledImage(pixel);
        regionLabels(i) = label;
    end
    % Get the centroids of each region
    segmentationXCentroids = NaN(length(regionProperties),1);
    segmentationYCentroids = NaN(length(regionProperties),1);
    for i = 1:length(regionProperties)
        centroid = regionProperties(i).Centroid;
        segmentationXCentroids(i) = centroid(1);
        segmentationYCentroids(i) = centroid(2);
    end
    % Create a struct with the region labels and centroids
    segmentationCentroids = struct('Label', regionLabels, 'XCentroid', segmentationXCentroids, 'YCentroid', segmentationYCentroids);
end