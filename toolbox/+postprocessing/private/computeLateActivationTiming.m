% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function [timeToHalf, peakPhase, maxSlopeLate] = computeLateActivationTiming(y, searchStart, lateWindow, nPhasePoints)

    timeToHalf = NaN;
    peakPhase = NaN;
    maxSlopeLate = NaN;

    if all(isnan(y))
        return;
    end

    [~, peakIdx] = max(y);
    peakPhase = peakIdx / nPhasePoints;

    halfMax = min(y, [], "omitnan") + 0.5 * (max(y, [], "omitnan") - min(y, [], "omitnan"));
    idxHalf = find(y(searchStart:end) >= halfMax, 1, "first");

    if ~isempty(idxHalf)
        timeToHalf = (searchStart + idxHalf - 1) / nPhasePoints;
    end

    if numel(lateWindow) >= 2
        dyLate = diff(y(lateWindow));
        maxSlopeLate = max(dyLate, [], "omitnan");
    end
end

