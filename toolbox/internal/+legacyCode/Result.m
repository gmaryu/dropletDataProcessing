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

        function save(obj, folderPath, fileName)
            %   save Saves the segmentation result to a .mat file at folderPath/fileName.m
            %
            %   Inputs:
            %       folderPath string - Path to the folder where the file should be saved
            %       fileName string - Name of the file to be saved
            %
            %   Note:
            %       If folderPath does not end with a slash, one will be added automatically
            %       Similarlym if fileName does not end with .m, it will be added automatically
            arguments
                obj legacyCode.Result
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
    end
end