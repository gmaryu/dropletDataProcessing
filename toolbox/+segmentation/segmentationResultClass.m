classdef segmentationResultClass
    % TODO: Add description
    properties
        labeledImage (:,:) {mustBeNumeric}
        regionProperties struct
    end

    methods
        function obj = segmentationResultClass(labeledImage, regionProperties)
            obj.labeledImage = labeledImage;
            obj.regionProperties = regionProperties;
        end
    end
end