function [Dxx, Dxy, Dyy] = hessian2D(image, segmentationParameters)
    %  hessian2D Filters the `image` using the second derivatives of a Gaussian with standard deviation `sigma`
    % 
    % outputs,
    %   Dxx, Dxy, Dyy: The 2nd derivatives, scale-corrected
    %
    % Function adapted from https://github.com/timjerman/JermanEnhancementFilter
    % by Zhengda Li first and then by Franco Tavella
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