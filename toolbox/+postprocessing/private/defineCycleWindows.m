% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function [windows, F] = defineCycleWindows(F, c, nucCycleInfo, XnucVol, nPhasePoints, P)

    windows = struct("skip", true);

    idx1 = nucCycleInfo.Peak1Frame(c);
    idx2 = nucCycleInfo.Peak2Frame(c);
    troughFrame = nucCycleInfo.TroughFrame(c);
    traceRow = nucCycleInfo.TraceRow(c);

    if isnan(idx1) || isnan(idx2) || isnan(troughFrame) || isnan(traceRow)
        return;
    end

    idx1 = round(idx1);
    idx2 = min(size(XnucVol, 2), round(idx2));
    troughFrame = max(idx1, min(idx2, round(troughFrame)));
    traceRow = round(traceRow);

    troughIdxNorm = frameToPhaseIndex(troughFrame, idx1, idx2, nPhasePoints);
    if isnan(troughIdxNorm)
        return;
    end

    interphaseEndFrame = detectInterphaseEndFrameFromNucVol(traceRow, idx1, idx2, troughFrame, troughIdxNorm, XnucVol, nPhasePoints, P);
    interphaseStart = frameToPhaseIndex(troughFrame, idx1, idx2, nPhasePoints);
    interphaseEnd = frameToPhaseIndex(interphaseEndFrame, idx1, idx2, nPhasePoints);

    interphaseStart = max(1, min(nPhasePoints, interphaseStart));
    interphaseEnd = max(interphaseStart, min(nPhasePoints, interphaseEnd));
    interphaseEnd = min(interphaseEnd, nPhasePoints - P.excludeLastNPhasePointsFromInterphase);
    interphaseEnd = max(interphaseStart, interphaseEnd);

    if interphaseEnd - interphaseStart + 1 < P.minWindowPoints
        interphaseStart = max(1, interphaseEnd - P.minWindowPoints + 1);
    end

    lateStart = min(interphaseEnd + 1, nPhasePoints);
    lateEnd = nPhasePoints;
    if lateEnd - lateStart + 1 < P.minWindowPoints
        lateStart = max(1, nPhasePoints - P.minWindowPoints + 1);
    end

    windows.skip = false;
    windows.troughIdxNorm = troughIdxNorm;
    windows.interphase = interphaseStart:interphaseEnd;
    windows.late = lateStart:lateEnd;

    F.TroughIdxNorm(c) = troughIdxNorm;
    F.InterphaseStartIdx(c) = interphaseStart;
    F.InterphaseEndIdx(c) = interphaseEnd;
    F.LateStartIdx(c) = lateStart;
    F.LateEndIdx(c) = lateEnd;
    F.InterphaseStartFrame(c) = troughFrame;
    F.InterphaseEndFrame(c) = interphaseEndFrame;

    if P.useNucVolForInterphaseEnd
        F.InterphaseEndReason(c) = getInterphaseEndReason(traceRow, idx1, idx2, troughFrame, XnucVol, P);
    else
        F.InterphaseEndReason(c) = "fixedOffset";
    end
end

