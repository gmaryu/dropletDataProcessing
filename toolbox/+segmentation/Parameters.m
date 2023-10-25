classdef Parameters
    %   Parameters Segmentation parameters for droplet detection spanning from image pre-processing to label filtering
    % 
    %   Constructor syntax:
    %       % Default values
    %       segmentationParameters = segmentation.Parameters()
    %       % Custom values for resolution, minDropletRadius and maxDropletRadius
    %       segmentationParameters = segmentation.Parameters(resolution, minDropletRadius, maxDropletRadius)
    %       % To change other parameters, use dot notation after creating the object
    %       segmentationParameters = segmentation.Parameters()
    %       segmentationParameters.maxEccentricity = 0.8;
    %
    %   Parameters properties:
    %       resolution - Resolution relative to 1x1 binning with 20x objective. Use 1 for 1x1 binning, 0.5 for 2x2 binning, etc.
    %       illuminationCorrectionFilterSize - Size of the Gaussian filter used for illumination correction.
    %       gaussianFilterSizeLarge - Size of the large Gaussian filter used for border enhancement.
    %       gaussianFilterSizeSmall - Size of the small Gaussian filter used for border enhancement.
    %       hessianSigma - Sigma of the Gaussian filter used for Hessian matrix calculation.
    %       hessianEigenvalueScale - Scale of the eigenvalues of the Hessian matrix used for border enhancement.
    %       thresholdFactor - Factor used to scale Otsu's threshold for border enhancement.
    %       borderSmallAreaThreshold - Threshold used to remove small regions after border enhancement.
    %       borderMediumAreaThreshold - Threshold used to characterize medium-sized regions after border enhancement.
    %       borderEccentricityCutoff - Eccentricity cutoff used to characterize regions after border enhancement.
    %       HMaximaSuppressionThreshold - Threshold used to suppress local maxima in the H-maxima transform.
    %       openingDiskSize - Size of the disk used for morphological opening after applying the watershed algorithm.
    %       minDropletRadius - Minimum droplet radius in um.
    %       maxDropletRadius - Maximum droplet radius in um.
    %       maxEccentricity - Maximum eccentricity of droplets. Value between 0 (perfect circle) and 1 (line segment).
    %       maxNumberOfHoles - Maximum number of holes in droplets. Related to EulerNumber region property.
    %       minCircularity - Minimum circularity of droplets. Value between 0 (line segment) and 1 (perfect circle).
    properties
        % Resolution relative to 1x1 binning with 20x objective. Use 1 for 1x1 binning, 0.5 for 2x2 binning, etc.
        resolution (1,1) double = 1; % Image property 
        % Size of the Gaussian filter used for illumination correction
        illuminationCorrectionFilterSize (1,1) double = 30; % preProcessBrightFieldImage
        % Size of the large Gaussian filter used for border enhancement
        gaussianFilterSizeLarge (1,1) double = 6; % enhanceDropletBorders
        % Size of the small Gaussian filter used for border enhancement
        gaussianFilterSizeSmall (1,1) double = 2; % enhanceDropletBorders
        % Sigma of the Gaussian filter used for Hessian matrix calculation
        hessianSigma (1,1) double = 3; % enhanceDropletBorders
        % Scale of the eigenvalues of the Hessian matrix used for border enhancement
        hessianEigenvalueScale (1,1) double = 0.0005; % enhanceDropletBorders
        % Factor used to scale Otsu's threshold for border enhancement
        thresholdFactor (1,1) double = 0.25; % enhanceDropletBorders
        % Threshold used to remove small regions after border enhancement
        borderSmallAreaThreshold (1,1) double = 50; % enhanceDropletBorders
        % Threshold used to characterize medium-sized regions after border enhancement
        borderMediumAreaThreshold (1,1) double = 1000; % enhanceDropletBorders
        % Eccentricity cutoff used to characterize regions after border enhancement
        borderEccentricityCutoff (1,1) double = 0.95; % enhanceDropletBorders
        % Threshold used to suppress local maxima in the H-maxima transform
        HMaximaSuppressionThreshold (1,1) double = 2; % getDropletMarkers
        % Size of the disk used for morphological opening after applying the watershed algorithm
        openingDiskSize (1,1) double = 2; % labelObjects
        % Minimum droplet radius in um
        minDropletRadius (1,1) double = 20; % filterLabels
        % Maximum droplet radius in um
        maxDropletRadius (1,1) double = 70; % filterLabels
        % Maximum eccentricity of droplets. Value between 0 (perfect circle) and 1 (line segment)
        maxEccentricity (1,1) double = 0.95; % filterLabels
        % Maximum number of holes in droplets. Related to EulerNumber region property
        maxNumberOfHoles (1,1) double = 0; % filterLabels
        % Minimum circularity of droplets. Value between 0 (line segment) and 1 (perfect circle)
        minCircularity (1,1) double = 0.75; % filterLabels
    end

    methods
        function obj = Parameters(resolution, minDropletRadius, maxDropletRadius)
            % Creates a new Parameters object. If no arguments are provided, default values are used
            % but resolution, minDropletRadius and maxDropletRadius can be provided as optional arguments
            arguments
                resolution (1,1) double = 1 % 1 for 1x1 binning, 0.5 for 2x2 binning, etc. Assumes 20x objective was used
                minDropletRadius (1,1) double = 20 % Minimum droplet radius in um
                maxDropletRadius (1,1) double = 70 % Maximum droplet radius in um
            end
            obj.resolution = resolution;
            obj.minDropletRadius = minDropletRadius;
            obj.maxDropletRadius = maxDropletRadius;
        end
    end
end