% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function idxNorm = frameToPhaseIndex(frame, peak1Frame, peak2Frame, nPhasePoints)

    if isnan(frame) || peak2Frame <= peak1Frame
        idxNorm = NaN;
        return;
    end

    phaseValue = (frame - peak1Frame) / (peak2Frame - peak1Frame);
    idxNorm = round(phaseValue * (nPhasePoints - 1)) + 1;
    idxNorm = max(1, min(nPhasePoints, idxNorm));
end

