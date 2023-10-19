classdef Droplet < handle
    % TODO: Add description
    properties
        id (1,1) 
        segmentationLabels (1,:) 
    end

    methods
        function obj = Droplet(id, segmentationLabels)
            % TODO: Add description
            arguments
                id (1,1) {mustBeInteger, mustBePositive}
                segmentationLabels (1,:) 
            end
            obj.id = id;
            obj.segmentationLabels = segmentationLabels;
        end

        function setAvgChannelValue(obj, channel)
            % TODO: Implement
            % Creates an attribute with the chosen channel's average value
        end

        function setStdChannelValue(obj, channel)
            % TODO: Implement
            % Creates an attribute with the chosen channel's standard deviation value
        end
    end

end