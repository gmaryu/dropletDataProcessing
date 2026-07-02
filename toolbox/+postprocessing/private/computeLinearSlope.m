% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function slope = computeLinearSlope(x, y)
    valid = ~isnan(x) & ~isnan(y);
    if sum(valid) >= 3
        p = polyfit(x(valid), y(valid), 1);
        slope = p(1);
    else
        slope = NaN;
    end
end

