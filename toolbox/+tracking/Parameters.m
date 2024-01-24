classdef Parameters
    %   Parameters Tracking parameters used in the analysis of segmented droplets
    %
    %   Constructor syntax:
    %       % Default values
    %       trackingParameters = tracking.Parameters();
    %       % Custom values
    %       trackingParameters = tracking.Parameters(maxCost, costOfNonAssignment, minTrackLength, mustStartOnFirstFrame);
    %
    %   Properties:
    %       maxCost - Maximum cost for linking two segments. Values higher than this in the cost matrix will be set to Inf.
    %       costOfNonAssignment - Cost of not linking a segment to any other segment.
    %       minTrackLength - Minimum length of a track (in number of frames). Tracks shorter than this will be discarded.
    %       mustStartOnFirstFrame - If true, only tracks that start on the first frame will be kept.
    properties
        % Max cost for linking two segments. Values higher than this in the cost matrix will be set to Inf.
        maxCost = 100; % linkSegmentationResults
        % Cost of not linking a segment to any other segment.
        costOfNonAssignment = 100; % linkSegmentationResults
        % Minimum length of a track (in number of frames). Tracks shorter than this will be discarded.
        minTrackLength = 3; % createDropletTracks
        % If true, only tracks that start on the first frame will be kept.
        mustStartOnFirstFrame = false; % createDropletTracks
    end

    methods
        function obj = Parameters(maxCost, costOfNonAssignment, minTrackLength, mustStartOnFirstFrame)
            %   Create a Parameters object. If no arguments are given, the default values are used.
            %
            %   Optional Inputs:
            %       maxCost (1,1) double = 100 - Maximum cost for linking two segments. Values higher than this in the cost matrix will be set to Inf.
            %       costOfNonAssignment (1,1) double = 100 - Cost of not linking a segment to any other segment.
            %       minTrackLength (1,1) integer = 3 - Minimum length of a track (in number of frames). Tracks shorter than this will be discarded.
            %       mustStartOnFirstFrame (1,1) logical = false - If true, only tracks that start on the first frame will be kept.
            arguments
                maxCost (1,1) double = 100
                costOfNonAssignment (1,1) double = 100
                minTrackLength (1,1) {mustBeInteger} = 3
                mustStartOnFirstFrame (1,1) logical = false
            end
            obj.maxCost = maxCost;
            obj.costOfNonAssignment = costOfNonAssignment;
            obj.minTrackLength = minTrackLength;
            obj.mustStartOnFirstFrame = mustStartOnFirstFrame;
        end
    end
end