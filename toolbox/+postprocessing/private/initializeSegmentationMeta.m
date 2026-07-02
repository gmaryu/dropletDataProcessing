% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function meta = initializeSegmentationMeta(n)
    meta.CycleLength = nan(n, 1);
    meta.Peak1Frame = nan(n, 1);
    meta.Peak2Frame = nan(n, 1);
    meta.TroughFrame = nan(n, 1);
    meta.TroughFrameRelative = nan(n, 1);
    meta.TroughPhase = nan(n, 1);
    meta.TotalPeak1Value = nan(n, 1);
    meta.TotalPeak2Value = nan(n, 1);
    meta.TotalTroughValue = nan(n, 1);
    meta.TotalAmplitudeFromTroughToPeak2 = nan(n, 1);
    meta.TotalAmplitudeMaxMinusMin = nan(n, 1);
    meta.NaNFractionTotal = nan(n, 1);
    meta.NaNFractionNucVol = nan(n, 1);
    meta.NaNFractionNCvolRatio = nan(n, 1);
end

