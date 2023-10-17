classdef Parameters
    % TODO: Add description
    properties
        % Pre-processing
        gaussianFilterSigma = 1;
        % Binary mask processing
        morphologicalOpeningRadius = 6;
        % Filtering of segmented image
        minSegmentedArea = 1000; % unit: number of pixels in the region
        maxSegmentedArea = 20000; % unit: number of pixels in the region
        maxEccentricity = 0.75; % value between 0 (perfect circle) and 1 (line segment)
        maxNumberOfHoles = 0; % unit: number of holes in the region
        minCircularity = 0.9; % value between 0 (line segment) and 1 (perfect circle)
    end

    methods
        function obj = Parameters()
            % TODO: Implement the constructor
        end
    end
end