classdef Result < handle
    %   Result Tracking result for storing droplet and segmentation link results
    %
    %   Constructor:
    %       trackingResult = tracking.Result(linkResultArray, dropletArray);
    %
    %   Properties:
    %       linkResultArray (1,:) tracking.LinkResult - Array of results linking segmented frames together
    %       dropletArray (1,:) tracking.Droplet - Array of tracked droplets.
    properties
        % Array containing the results of linking segmented frames together
        linkResultArray (1,:) tracking.LinkResult
        % Array containing the tracked droplets
        dropletArray (1,:) tracking.Droplet
    end

    methods
        function obj = Result(linkResultArray, dropletArray)
            %   Creates a new Result object
            %
            %   Inputs:
            %       linkResultArray (1,:) tracking.LinkResult - Array of results linking segmented frames together
            %       dropletArray (1,:) tracking.Droplet - Array of tracked droplets.
            arguments
                linkResultArray (1,:) tracking.LinkResult = tracking.LinkResult.empty
                dropletArray (1,:) tracking.Droplet = tracking.Droplet.empty
            end
            obj.linkResultArray = linkResultArray;
            obj.dropletArray = dropletArray;
        end
    end
end