% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function yNorm = amplitudeNormalize(y)

    yMin = min(y, [], "omitnan");
    yMax = max(y, [], "omitnan");

    if yMax > yMin
        yNorm = (y - yMin) ./ (yMax - yMin);
    else
        yNorm = nan(size(y));
    end
end

