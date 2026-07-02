% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function [nucCycleInfo, cycles] = segmentPeakToPeakCycles(cycleBoundaryTable, X, nDroplets, nTime)

    nCyclesTotal = height(cycleBoundaryTable);
    cycles = initializeCycleCells(nCyclesTotal);
    nucCycleInfo = cycleBoundaryTable;

    meta = initializeSegmentationMeta(nCyclesTotal);

    for c = 1:nCyclesTotal
        row = nucCycleInfo.TraceRow(c);
        idx1 = nucCycleInfo.Peak1Frame(c);
        idx2 = nucCycleInfo.Peak2Frame(c);
        troughIdx = nucCycleInfo.TroughFrame(c);

        if ~isUsableBoundary(row, idx1, idx2, nDroplets, nTime)
            continue;
        end

        row = round(row);
        idx1 = round(idx1);
        idx2 = round(idx2);
        frameIdx = idx1:idx2;

        cycles.total{c} = X.total(row, frameIdx);
        cycles.nuc{c} = X.nuc(row, frameIdx);
        cycles.cyto{c} = X.cyto(row, frameIdx);
        cycles.nucVol{c} = X.nucVol(row, frameIdx);
        cycles.dropletVol{c} = X.dropletVol(row, frameIdx);
        cycles.NCvolRatio{c} = X.NCvolRatio(row, frameIdx);
        cycles.nucCytoRatio{c} = X.nucCytoRatio(row, frameIdx);
        cycles.nucMinusCyto{c} = X.nucMinusCyto(row, frameIdx);

        meta = updateSegmentationMeta(meta, c, row, idx1, idx2, troughIdx, X.total, frameIdx);
        meta.NaNFractionNucVol(c) = mean(isnan(cycles.nucVol{c}));
        meta.NaNFractionNCvolRatio(c) = mean(isnan(cycles.NCvolRatio{c}));
    end

    nucCycleInfo = addStructFieldsToTable(nucCycleInfo, meta);
end

