classdef Result < handle
    % TODO: Add description
    properties
        labeledImage (:,:) {mustBeNumeric}
        regionProperties struct
    end

    methods
        function obj = Result(labeledImage, regionProperties)
            % TODO: Add description
            obj.labeledImage = labeledImage;
            obj.regionProperties = regionProperties;
        end

        function plot(obj, originalImage)
            % TODO: Add description
            % Adjust contrast for better visualization
            image = histeq(originalImage); 
            overlay = labeloverlay(image, obj.labeledImage);
            imshow(overlay);
        end
    end
end