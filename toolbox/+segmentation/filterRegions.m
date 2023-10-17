function [filteredLabeledImage, filteredRegionProperties] = filterRegions(brightFieldImage, labeledImage, segmentationParameters)
    % TODO: Add description
    arguments
        brightFieldImage (:,:) {mustBeNumeric}
        labeledImage (:,:) {mustBeNumeric}
        segmentationParameters segmentation.Parameters
    end

    % Calculate the properties of each region. For an explanation on each property, see
    % https://www.mathworks.com/help/images/ref/regionprops.html
    regionProperties = regionprops(labeledImage, brightFieldImage, ...
                                   'Centroid','Area', 'Circularity',...
                                   'PixelIdxList','PixelValues', ...
                                   'Eccentricity','EulerNumber');
    % Filter segmentation result
    tooSmall = [regionProperties.Area] < segmentationParameters.minSegmentedArea;
    tooLarge = [regionProperties.Area] > segmentationParameters.maxSegmentedArea;
    notSymmetrical = [regionProperties.Eccentricity] > segmentationParameters.maxEccentricity;
    notRound = [regionProperties.Circularity] < segmentationParameters.minCircularity;
    tooManyHoles = [regionProperties.EulerNumber] < 1 - segmentationParameters.maxNumberOfHoles;
    filter = tooSmall | tooLarge | notSymmetrical | notRound | tooManyHoles;
    % Remove regions that do not pass the filter from the labeled image and region properties
    filteredLabeledImage = labeledImage;
    filteredLabeledImage(vertcat(regionProperties(filter).PixelIdxList)) = 0;
    filteredRegionProperties = regionProperties(~filter); 
end