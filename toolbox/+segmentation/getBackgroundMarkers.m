function backgroundMarkers = getBackgroundMarkers(maxCorrectedBrightField, dropletMarkers)
    %   getBackgroundMarkers Calculate the location of the background using ridge lines of `dropletMarkers`
    % 
    %   backgroundMarkers = getBackgroundMarkers(maxCorrectedBrightField, dropletMarkers)
    %   
    %   Inputs:
    %       maxCorrectedBrightField (MxN) double - Brightfield image obtained via `getMaxCorrectedBrightField`
    %       dropletMarkers (MxN) logical - Droplet markers obtained via `getDropletMarkers`
    %
    %   Output:
    %       backgroundMarkers (MxN) logical - Binary image with background markers set to 1
    arguments
        maxCorrectedBrightField (:,:) {mustBeNumeric} % Brightfield image obtained via `getMaxCorrectedBrightField`
        dropletMarkers (:,:) logical % Droplet markers obtained via `getDropletMarkers`
    end
    minImposedImage = imimposemin(maxCorrectedBrightField, dropletMarkers);
    labels = watershed(minImposedImage);
    backgroundMarkers = labels == 0;
end