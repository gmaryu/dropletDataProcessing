% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function [nucCycleInfo, ctx] = extractCycleFeatures(nucCycleInfo, normCycles, shapeCycles, X, cytoRef, phase, nPhasePoints, P)

    n = height(nucCycleInfo);
    F = initializeFeatureStruct(n);

    for c = 1:n
        F = computeCytoplasmicDistance(F, c, nucCycleInfo.CycleInDroplet(c), shapeCycles.total(c, :), cytoRef, P);

        [windows, F] = defineCycleWindows(F, c, nucCycleInfo, X.nucVol, nPhasePoints, P);
        if windows.skip
            continue;
        end

        F = computeCdk1WindowFeatures(F, c, normCycles, shapeCycles, phase, windows, P);
        F = computeVolumeWindowFeatures(F, c, normCycles, windows);
    end

    nucCycleInfo = addStructFieldsToTable(nucCycleInfo, F);

    ctx.interphaseWindow = [];
    ctx.lateWindow = [];
end

