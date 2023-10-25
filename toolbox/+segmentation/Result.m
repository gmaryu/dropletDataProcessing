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
            arguments
                labeledImage (:,:) {mustBeNumeric};
                regionProperties (:,:) struct;
            end
            % Creates a segmentation result object
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
            % plot Plots the segmentation result on top of the original image as an overlay
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