classdef Result < handle
    % TODO: Add description
    properties
        labeledImage (:,:) {mustBeNumeric}
        segmentationCentroids (:,:) struct
    end

    methods
        function obj = Result(labeledImage, segmentationCentroids)
            % TODO: Add description
            obj.labeledImage = labeledImage;
            obj.segmentationCentroids = segmentationCentroids;
        end

        function save(obj, folderPath, fileName)
            % TODO: Add description
            % See advantages of using v7.3 here: https://www.mathworks.com/help/matlab/import_export/load-parts-of-variables-from-mat-files.html
            arguments
                obj segmentation.Result
                folderPath string
                fileName string
            end
            % check if fileName contains extension
            if ~endsWith(fileName, '.mat')
                fileName = fileName + '.mat';
            end
            fullPath = folderPath + fileName;
            save (fullPath, 'obj', '-v7.3');
        end

        function plot(obj, originalImage)
            % TODO: Add description
            arguments
                obj segmentation.Result
                originalImage (:,:) {mustBeNumeric}
            end
            % Adjust contrast for better visualization
            image = histeq(originalImage); 
            overlay = labeloverlay(image, obj.labeledImage);
            imshow(overlay);

        end
    end
end