% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function F = computeCdk1WindowFeatures(F, c, normCycles, shapeCycles, phase, windows, P)

    iw = windows.interphase;
    lw = windows.late;

    F.InterphaseNucMean(c) = mean(normCycles.nuc(c, iw), "omitnan");
    F.InterphaseCytoMean(c) = mean(normCycles.cyto(c, iw), "omitnan");
    F.InterphaseNucCytoRatioMean(c) = mean(normCycles.nucCytoRatio(c, iw), "omitnan");
    F.InterphaseNucMinusCytoMean(c) = mean(normCycles.nucMinusCyto(c, iw), "omitnan");

    ratioTrace = normCycles.nucCytoRatio(c, iw);
    F.InterphaseNucCytoRatioMax(c) = max(ratioTrace, [], "omitnan");
    F.InterphaseNucCytoRatioP90(c) = prctileOmitNaN(ratioTrace, P.interphasePercentileForState2);

    [earlyW, lateInterphaseW] = splitWindowInHalf(iw);
    F.EarlyInterphaseNucCytoRatioMean(c) = mean(normCycles.nucCytoRatio(c, earlyW), "omitnan");
    F.LateInterphaseNucCytoRatioMean(c) = mean(normCycles.nucCytoRatio(c, lateInterphaseW), "omitnan");
    F.InterphaseNucCytoRatioIncrease(c) = ...
        F.LateInterphaseNucCytoRatioMean(c) - F.EarlyInterphaseNucCytoRatioMean(c);

    F.InterphaseNucCytoRatioSlope(c) = computeLinearSlope(phase(iw)', ratioTrace(:));

    F.LateNucMean(c) = mean(normCycles.nuc(c, lw), "omitnan");
    F.LateCytoMean(c) = mean(normCycles.cyto(c, lw), "omitnan");
    F.LateActivationScore_nucShape(c) = ...
        mean(shapeCycles.nuc(c, lw), "omitnan") - mean(shapeCycles.nuc(c, iw), "omitnan");

    yNucShape = shapeCycles.nuc(c, :);
    [F.TimeToHalfActivation_nuc(c), F.PeakPhase_nuc(c), F.MaxSlopeLate_nucShape(c)] = ...
        computeLateActivationTiming(yNucShape, windows.troughIdxNorm, lw, numel(yNucShape));
end

