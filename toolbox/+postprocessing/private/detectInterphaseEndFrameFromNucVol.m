% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function endFrame = detectInterphaseEndFrameFromNucVol(traceRow, idx1, idx2, troughFrame, troughIdxNorm, XnucVol, nPhasePoints, P)

    fallbackEndIdx = min(troughIdxNorm + P.interphaseOffsetAfterTrough, ...
        nPhasePoints - P.excludeLastNPhasePointsFromInterphase);
    fallbackEndIdx = max(troughIdxNorm, fallbackEndIdx);
    fallbackEndFrame = phaseIndexToFrame(fallbackEndIdx, idx1, idx2, nPhasePoints);

    if ~P.useNucVolForInterphaseEnd
        endFrame = fallbackEndFrame;
        return;
    end

    nucVolCycleRaw = XnucVol(traceRow, idx1:idx2);
    frameIdxCycle = idx1:idx2;
    searchStartFrame = max(idx1, min(idx2, troughFrame + P.minFramesAfterTroughForLossSearch));
    searchMask = frameIdxCycle >= searchStartFrame;

    nucVolSearch = nucVolCycleRaw(searchMask);
    frameSearch = frameIdxCycle(searchMask);

    if P.useRelativeNucVolLossThreshold
        nucVolCycleMax = max(nucVolCycleRaw, [], "omitnan");
        nucVolLossThresholdThisCycle = P.relativeNucVolLossThreshold * nucVolCycleMax;
    else
        nucVolLossThresholdThisCycle = P.nucVolLossThreshold;
    end

    lostNucleus = isnan(nucVolSearch) | nucVolSearch <= nucVolLossThresholdThisCycle;
    lossLocalIdx = firstConsecutiveTrue(lostNucleus, P.nucVolLossConsecutiveN);

    if ~isnan(lossLocalIdx)
        firstLossFrame = frameSearch(lossLocalIdx);
        endFrame = max(troughFrame, min(idx2, firstLossFrame - 1));
    else
        endFrame = fallbackEndFrame;
    end
end

