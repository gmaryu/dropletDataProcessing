% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function tf = isUsableBoundary(row, idx1, idx2, nDroplets, nTime)
    tf = ~(isnan(row) || isnan(idx1) || isnan(idx2) || idx2 <= idx1);
    if ~tf
        return;
    end
    row = round(row);
    idx1 = round(idx1);
    idx2 = round(idx2);
    tf = row >= 1 && row <= nDroplets && idx1 >= 1 && idx2 <= nTime;
end

