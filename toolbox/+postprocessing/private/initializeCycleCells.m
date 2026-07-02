% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function cycles = initializeCycleCells(n)
    names = ["total", "nuc", "cyto", "nucVol", "dropletVol", "NCvolRatio", "nucCytoRatio", "nucMinusCyto"];
    for i = 1:numel(names)
        cycles.(names(i)) = cell(n, 1);
    end
end

