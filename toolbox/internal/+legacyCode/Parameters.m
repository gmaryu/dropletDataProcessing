classdef Parameters
    % TODO: Add description
    properties
        % Binary mask processing
        c = 0.0005;
        gaussianFilterSizeLarge = 6;
        gaussianFilterSizeSmall = 2;
        thresholdFactor = 0.25;
        hessianSigma = 3;
        resAmp = 1.0;
        maxEccentricity = 0.95;
        minSegmentedArea = 1000;
        prthresh = 0.35;
        % Filtering of segmented image
        maxSegmentedArea = 10000; % unit: number of pixels in the region
        maxNumberOfHoles = 0; % unit: number of holes in the region
        minCircularity = 0.9; % value between 0 (line segment) and 1 (perfect circle)
    end

    methods
        function obj = Parameters()
            % TODO: Implement the constructor
        end
    end
end