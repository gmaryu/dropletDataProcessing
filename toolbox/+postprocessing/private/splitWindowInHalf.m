% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function [earlyW, lateW] = splitWindowInHalf(w)
    nW = numel(w);
    splitPoint = max(1, floor(nW / 2));
    earlyW = w(1:splitPoint);
    lateW = w(min(splitPoint + 1, nW):end);
end

