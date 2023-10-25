classdef Result < handle
    %   Result Class for storing segmentation results as a labeled image and region properties
    %
    %   Constructor:
    %       segmentationResult = segmentation.Result(labeledImage, regionProperties)
    %
    %   Methods:
    %       segmentationResult.save(folderPath, fileName) - Saves the segmentation result to a .mat file at folderPath/fileName
    %       segmentationResult.plot(originalImage) - Plots the labeled image on top of the provided image as an overlay
    properties
        labeledImage (:,:) {mustBeNumeric}
        regionProperties (:,:) struct
    end

    methods
        function obj = Result(labeledImage, regionProperties)
            %  Creates a segmentation result object
            %
            %   Inputs:
            %       labeledImage - A labeled image where each pixel is assigned a label
            %       regionProperties - A struct array containing the properties of each region
            arguments
                labeledImage (:,:) {mustBeNumeric};
                regionProperties (:,:) struct;
            end
            obj.labeledImage = labeledImage;
            obj.regionProperties = regionProperties;
        end

        function save(obj, folderPath, fileName)
            %   save Saves the segmentation result to a .mat file at folderPath/fileName.m
            %
            %   Inputs:
            %       folderPath - Path to the folder where the file should be saved
            %       fileName - Name of the file to be saved
            %
            %   Note:
            %       If folderPath does not end with a slash, one will be added automatically
            %       Similarlym if fileName does not end with .m, it will be added automatically
            arguments
                obj segmentation.Result
                folderPath string
                fileName string
            end
            % check if folderPath contains trailing slash
            if ~endsWith(folderPath, '/')
                folderPath = folderPath + '/';
            end
            % check if fileName contains extension
            if ~endsWith(fileName, '.mat')
                fileName = fileName + '.mat';
            end
            fullPath = folderPath + fileName;
            save (fullPath, 'obj', '-v7.3');
        end

        function plot(obj, originalImage)
            %   plot Plots the segmentation result on top of the original image as an overlay
            %
            %   Inputs:
            %       originalImage - The original image on top of which the segmentation result should be plotted
            %
            %   Notes:
            %       - The original image is adjusted for better visualization with `histeq`
            %       - If the plot does not show up, try calling `figure` before calling this method or use set(gcf, 'Visible', 'on')
            arguments
                obj segmentation.Result
                originalImage (:,:) {mustBeNumeric}
            end
            image = histeq(originalImage); % Adjust contrast for better visualization
            overlay = labeloverlay(image, obj.labeledImage);
            imshow(overlay);
        end
    end
end