classdef segmentationParametersClass
    % TODO: Add description
    properties
        minSegmentedArea = 100; % unit: number of pixels in the region
        maxSegmentedArea = 10000; % unit: number of pixels in the region
        maxEccentricity = 0.8; % value between 0 (perfect circle) and 1 (line segment)
        maxNumberOfHoles = 0; % unit: number of holes in the region
    end

    methods
        function obj = segmentationParametersClass()
            % TODO: Implement the constructor
        end
    end
end