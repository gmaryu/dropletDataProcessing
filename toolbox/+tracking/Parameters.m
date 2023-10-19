classdef Parameters
    % TODO: Add description
    properties
        % linkSegmentattionResults
        maxCost = 100;
        costOfNonAssignment = 100;
        % createDropletTracks
        minTrackLength = 3;
        mustStartOnFirstFrame = false;
    end

    methods
        function obj = Parameters()
            % TODO: Implement the constructor
        end
    end
end