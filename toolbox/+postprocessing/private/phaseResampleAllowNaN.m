% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function yInterp = phaseResampleAllowNaN(y, targetPhase)

    y = y(:)';
    n = numel(y);

    if n < 2 || all(isnan(y))
        yInterp = nan(size(targetPhase));
        return;
    end

    oldPhase = linspace(0, 1, n);
    validIdx = ~isnan(y);

    if sum(validIdx) < 2
        yInterp = nan(size(targetPhase));
        return;
    end

    yInterp = interp1(oldPhase(validIdx), y(validIdx), targetPhase, "linear", NaN);
end

