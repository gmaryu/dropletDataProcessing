% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function frame = phaseIndexToFrame(idxNorm, peak1Frame, peak2Frame, nPhasePoints)

    if isnan(idxNorm) || nPhasePoints <= 1 || peak2Frame <= peak1Frame
        frame = NaN;
        return;
    end

    phaseValue = (idxNorm - 1) / (nPhasePoints - 1);
    frame = round(peak1Frame + phaseValue * (peak2Frame - peak1Frame));
    frame = max(peak1Frame, min(peak2Frame, frame));
end

