classdef Result
    % TODO: Add description
    properties
        labeledImage (:,:) {mustBeNumeric}
        regionProperties struct
    end

    methods
        function obj = Result(labeledImage, regionProperties)
            obj.labeledImage = labeledImage;
            obj.regionProperties = regionProperties;
        end
    end
end