% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function featuresRaw = smoothAndDeriveSignals(raw, P)

    featuresRaw.totalDetect = smoothTraceMatrixNoBridge(raw.total, P.detectSmoothWindow);

    featuresRaw.total = smoothTraceMatrixNoBridge(raw.total, P.featureSmoothWindow);
    featuresRaw.nuc = smoothTraceMatrixNoBridge(raw.nuc, P.featureSmoothWindow);
    featuresRaw.cyto = smoothTraceMatrixNoBridge(raw.cyto, P.featureSmoothWindow);
    featuresRaw.nucVol = smoothTraceMatrixNoBridge(raw.nucVol, P.featureSmoothWindow);
    featuresRaw.dropletVol = smoothTraceMatrixNoBridge(raw.dropletVol, P.featureSmoothWindow);
    featuresRaw.NCvolRatio = smoothTraceMatrixNoBridge(raw.NCvolRatio, P.featureSmoothWindow);

    epsilon = 1e-6;
    featuresRaw.nucCytoRatio = featuresRaw.nuc ./ (featuresRaw.cyto + epsilon);
    featuresRaw.nucMinusCyto = featuresRaw.nuc - featuresRaw.cyto;
end

