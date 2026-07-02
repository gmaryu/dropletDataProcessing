% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function reason = getInterphaseEndReason(traceRow, idx1, idx2, troughFrame, XnucVol, P)

    nucVolCycleRaw = XnucVol(traceRow, idx1:idx2);
    frameIdxCycle = idx1:idx2;
    searchStartFrame = max(idx1, min(idx2, troughFrame + P.minFramesAfterTroughForLossSearch));
    searchMask = frameIdxCycle >= searchStartFrame;
    nucVolSearch = nucVolCycleRaw(searchMask);

    if P.useRelativeNucVolLossThreshold
        nucVolCycleMax = max(nucVolCycleRaw, [], "omitnan");
        threshold = P.relativeNucVolLossThreshold * nucVolCycleMax;
    else
        threshold = P.nucVolLossThreshold;
    end

    lostNucleus = isnan(nucVolSearch) | nucVolSearch <= threshold;
    if ~isnan(firstConsecutiveTrue(lostNucleus, P.nucVolLossConsecutiveN))
        reason = "nucVolLoss";
    else
        reason = "fallbackFixedOffset";
    end
end

