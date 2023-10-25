function backgroundMarkers = getBackgroundMarkers(maxCorrectedBrightField, dropletMarkers)
    % getBackgroundMarkers Calculate location of background using ridge lines of `dropletMarkers`
    arguments
        maxCorrectedBrightField (:,:) {mustBeNumeric} % Brightfield image obtained via `getMaxCorrectedBrightField`
        dropletMarkers (:,:) logical % Droplet markers obtained via `getDropletMarkers`
    end
    minImposedImage = imimposemin(maxCorrectedBrightField, dropletMarkers);
    labels = watershed(minImposedImage);
    backgroundMarkers = labels == 0;
end