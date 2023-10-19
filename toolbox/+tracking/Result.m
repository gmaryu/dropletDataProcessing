classdef Result
    % TODO: Add description
    properties
        linkResultArray (1, :) tracking.LinkResult % Array of link results
        dropletArray (1, :) tracking.Droplet % Array of tracked droplets
    end
    methods
        function obj = Result(linkResultArray, dropletArray)
            arguments
                linkResultArray (1, :) tracking.LinkResult
                dropletArray (1, :) tracking.Droplet
            end
            obj.linkResultArray = linkResultArray;
            obj.dropletArray = dropletArray;
        end
    end
end