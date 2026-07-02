% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function q = prctileOmitNaN(y, p)

    y = y(:);
    y = y(~isnan(y));

    if isempty(y)
        q = NaN;
    else
        q = prctile(y, p);
    end
end

