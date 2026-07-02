% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function F = computeVolumeWindowFeatures(F, c, normCycles, windows)

    iw = windows.interphase;

    F.InterphaseNucVolMean(c) = mean(normCycles.nucVol(c, iw), "omitnan");
    F.InterphaseDropletVolMean(c) = mean(normCycles.dropletVol(c, iw), "omitnan");
    F.InterphaseNCvolRatioMean(c) = mean(normCycles.NCvolRatio(c, iw), "omitnan");

    F.CycleMeanNucVol(c) = mean(normCycles.nucVol(c, :), "omitnan");
    F.CycleMeanNCvolRatio(c) = mean(normCycles.NCvolRatio(c, :), "omitnan");

    F.NucVolValidFraction(c) = mean(~isnan(normCycles.nucVol(c, :)));
    F.NCvolRatioValidFraction(c) = mean(~isnan(normCycles.NCvolRatio(c, :)));
    F.InterphaseNucVolValidFraction(c) = mean(~isnan(normCycles.nucVol(c, iw)));
    F.InterphaseNCvolRatioValidFraction(c) = mean(~isnan(normCycles.NCvolRatio(c, iw)));
end

%% =========================
%  Local functions: classification
%  =========================

