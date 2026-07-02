% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function firstIdx = firstConsecutiveTrue(logicalVec, consecutiveN)

    logicalVec = logicalVec(:)';
    firstIdx = NaN;

    if isempty(logicalVec) || consecutiveN < 1 || numel(logicalVec) < consecutiveN
        return;
    end

    for i = 1:(numel(logicalVec) - consecutiveN + 1)
        if all(logicalVec(i:(i + consecutiveN - 1)))
            firstIdx = i;
            return;
        end
    end
end

