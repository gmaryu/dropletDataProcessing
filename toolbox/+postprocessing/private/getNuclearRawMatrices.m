% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function raw = getNuclearRawMatrices(nucData)

    raw.total = nucData.WholeMean_RATIO;
    raw.nuc = nucData.nucMean_RATIO;
    raw.cyto = nucData.cytMean_RATIO;
    raw.nucVol = nucData.nucVol;
    raw.dropletVol = nucData.dropletVol;
    raw.NCvolRatio = nucData.NCVolRatio;
end

