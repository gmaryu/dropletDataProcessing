classdef Parameters
    %   Parameters Tracking parameters used in the analysis of segmented droplets
    %
    %   Constructor syntax:
    %       % Default values
    %       trackingParameters = tracking.Parameters(minTrackLength, mustStartOnFirstFrame);
    %       % Custom values
    %       trackingParameters = tracking.Parameters(minTrackLength, mustStartOnFirstFrame, maxCost, costOfNonAssignment);
    %
    %   Properties:
    %       maxCost - Maximum cost for linking two segments. Values higher than this in the cost matrix will be set to Inf.
    %       costOfNonAssignment - Cost of not linking a segment to any other segment.
    %       minTrackLength - Minimum length of a track (in number of frames). Tracks shorter than this will be discarded.
    %       mustStartOnFirstFrame - If true, only tracks that start on the first frame will be kept.
    properties
        % Minimum length of a track (in number of frames). Tracks shorter than this will be discarded.
        minTrackLength; % createDropletTracks
        % If true, only tracks that start on the first frame will be kept.
        mustStartOnFirstFrame; % createDropletTracks
        % Max cost for linking two segments. Values higher than this in the cost matrix will be set to Inf.
        maxCost = 100; % linkSegmentationResults
        % Cost of not linking a segment to any other segment.
        costOfNonAssignment = 100; % linkSegmentationResults
    end

    methods
        function obj = Parameters(minTrackLength, mustStartOnFirstFrame, maxCost, costOfNonAssignment)
            %   Create a Parameters object. If no arguments are given, the default values are used.
            %
            %   Required Inputs:
            %       minTrackLength (1,1) integer - Minimum length of a track (in number of frames). Tracks shorter than this will be discarded.
            %       mustStartOnFirstFrame (1,1) logical - If true, only tracks that start on the first frame will be kept.
            %
            %   Optional Inputs:
            %       maxCost (1,1) double = 100 - Maximum cost for linking two segments. Values higher than this in the cost matrix will be set to Inf.
            %       costOfNonAssignment (1,1) double = 100 - Cost of not linking a segment to any other segment.
            arguments
                minTrackLength (1,1) {mustBeInteger}
                mustStartOnFirstFrame (1,1) logical
                maxCost (1,1) double = 100
                costOfNonAssignment (1,1) double = 100
            end
            obj.maxCost = maxCost;
            obj.costOfNonAssignment = costOfNonAssignment;
            obj.minTrackLength = minTrackLength;
            obj.mustStartOnFirstFrame = mustStartOnFirstFrame;
        end
    end
end