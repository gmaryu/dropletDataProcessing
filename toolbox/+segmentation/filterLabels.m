function [filteredLabeledImage, filteredRegionProperties] = filterLabels(labeledImage, segmentationParameters)
    %   filterLabels Extract only droplet labels by filtering out regions that do not meet criteria for area, eccentricity, circularity and number of holes.
    % 
    %   [filteredLabeledImage, filteredRegionProperties] = filterLabels(labeledImage, segmentationParameters)
    %
    %   Inputs:
    %       labeledImage (MxN) double - labeled image of droplets created by `getDropletLabels`
    %       segmentationParameters (1x1) segmentation.Parameters - Parameters for segmentation
    %   
    %   Outputs:
    %       filteredLabeledImage (MxN) double - New labeled image with only those labels that pass the filter
    %       filteredRegionProperties (:,:) struct - Properties of the filtered labels (Label, Area, Centroid, Circularity, Eccentricity)
    arguments
        labeledImage (:,:) {mustBeNumeric}
        segmentationParameters segmentation.Parameters
    end
    regionProperties = regionprops(labeledImage, ...
                                   'Centroid','Area', 'Circularity',...
                                   'PixelIdxList', 'Eccentricity','EulerNumber');
    % Define filter conditions
    minRadiusInPixels = segmentationParameters.minDropletRadius * segmentationParameters.resolution * MicroscopeInfo.micronToPixel;
    minArea = pi * minRadiusInPixels^2;
    maxRadiusInPixels = segmentationParameters.maxDropletRadius * segmentationParameters.resolution * MicroscopeInfo.micronToPixel;
    maxArea = pi * maxRadiusInPixels^2;
    tooSmall = [regionProperties.Area] < minArea;
    tooLarge = [regionProperties.Area] > maxArea;
    notSymmetrical = [regionProperties.Eccentricity] > segmentationParameters.maxEccentricity;
    notRound = [regionProperties.Circularity] < segmentationParameters.minCircularity;
    tooManyHoles = [regionProperties.EulerNumber] < 1 - segmentationParameters.maxNumberOfHoles;
    filter = tooSmall | tooLarge | notSymmetrical | notRound | tooManyHoles;
    % Remove regions that do not pass the filter from the labeled image and region properties
    filteredLabeledImage = labeledImage;
    filteredLabeledImage(vertcat(regionProperties(filter).PixelIdxList)) = 0;
    filteredRegionProperties = regionProperties(~filter); 
    % Re-label the filtered image to avoid gaps in label numbers
    filteredLabeledImage = bwlabel(filteredLabeledImage);
    % Assign new label numbers to region properties
    for i = 1:length(filteredRegionProperties)
        pixelOfNewLabel = filteredRegionProperties(i).PixelIdxList(1);
        newLabel = filteredLabeledImage(pixelOfNewLabel);
        filteredRegionProperties(i).Label = newLabel;
    end
    % Remove PixelIdxList and EulerNumber field from region properties
    filteredRegionProperties = rmfield(filteredRegionProperties, {'PixelIdxList', 'EulerNumber'});
    % Sort region properties by Label
    filteredRegionProperties = struct2table(filteredRegionProperties);
    filteredRegionProperties = sortrows(filteredRegionProperties, 'Label');
    filteredRegionProperties = table2struct(filteredRegionProperties);
    % Reorder fields
    columnOrder = {'Label', 'Area', 'Centroid', 'Circularity', 'Eccentricity'};
    filteredRegionProperties = orderfields(filteredRegionProperties, columnOrder);
end