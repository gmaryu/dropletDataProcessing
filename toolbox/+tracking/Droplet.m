classdef Droplet < handle
    %   Droplet Class for storing information about tracked droplets
    %
    %   Constructor:
    %       tracking.Droplet(id, segmentationLabels)
    %
    %   Methods:
    %       
    properties
        id (1,1) 
        segmentationLabels (1,:) 
    end

    methods
        function obj = Droplet(id, segmentationLabels)
            %   Creates a Droplet object
            %
            %   Inputs:
            %       id (1,1) positive integer - Droplet's id
            %       segmentationLabels (1,:) array - Droplet's segmentation labels for each frame
            arguments
                id (1,1) {mustBeInteger, mustBePositive} = 1;
                segmentationLabels (1,:) = [];
            end
            obj.id = id;
            obj.segmentationLabels = segmentationLabels;
        end
    end
end