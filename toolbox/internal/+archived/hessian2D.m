function [Dxx, Dxy, Dyy] = hessian2D(image, segmentationParameters)
    %   hessian2D Filters the `image` using the second derivatives of a Gaussian with standard deviation `sigma`
    % 
    %   [Dxx, Dxy, Dyy] = hessian2D(image, sigma)
    %
    %   Inputs:
    %       image (MxN) double - Input image to be processed
    %       segmentation.Parameters (1x1) segmentation.Parameters - Parameters of the segmentation
    %
    %   Outputs:
    %    Dxx (MxN) double - Second derivative in x, scale corrected
    %    Dxy (MxN) double - Second derivative in x and y, scale corrected
    %    Dyy (MxN) double - Second derivative in y, scale corrected
    %
    %   Note:
    %       Function adapted from https://github.com/timjerman/JermanEnhancementFilter
    %       first by Zhengda Li and then by Franco Tavella
    arguments
        image (:,:) {mustBeNumeric}
        segmentationParameters (1,1) segmentation.Parameters
    end
    sigma = segmentationParameters.hessianSigma;
    if nargin < 2, sigma = 1; end
    % Make kernel coordinates
    [X,Y] = ndgrid(-round(3*sigma):round(3*sigma));
    % Gaussian 2nd derivatives filters
    DGaussxx = 1/(2*pi*sigma^4) * (X.^2/sigma^2 - 1) .* exp(-(X.^2 + Y.^2)/(2*sigma^2));
    DGaussxy = 1/(2*pi*sigma^6) * (X .* Y) .* exp(-(X.^2 + Y.^2)/(2*sigma^2));
    DGaussyy = DGaussxx';
    Dxx = imfilter(image, DGaussxx, 'conv');
    Dxy = imfilter(image, DGaussxy, 'conv');
    Dyy = imfilter(image, DGaussyy, 'conv');
    % Correct for scale
    Dxx = (segmentationParameters.hessianSigma^2) * Dxx; 
    Dxy = (segmentationParameters.hessianSigma^2) * Dxy;
    Dyy = (segmentationParameters.hessianSigma^2) * Dyy;
end