% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function meta = updateSegmentationMeta(meta, c, row, idx1, idx2, troughIdx, Xtotal, frameIdx)

    thisTotal = Xtotal(row, frameIdx);

    meta.CycleLength(c) = numel(frameIdx);
    meta.Peak1Frame(c) = idx1;
    meta.Peak2Frame(c) = idx2;
    meta.TotalPeak1Value(c) = Xtotal(row, idx1);
    meta.TotalPeak2Value(c) = Xtotal(row, idx2);

    if ~isnan(troughIdx)
        troughIdx = round(troughIdx);
        meta.TroughFrame(c) = troughIdx;

        if troughIdx >= idx1 && troughIdx <= idx2
            troughRel = troughIdx - idx1 + 1;
            meta.TroughFrameRelative(c) = troughRel;
            meta.TroughPhase(c) = troughRel / numel(frameIdx);
            meta.TotalTroughValue(c) = Xtotal(row, troughIdx);
        end
    end

    if ~isnan(meta.TotalTroughValue(c))
        meta.TotalAmplitudeFromTroughToPeak2(c) = meta.TotalPeak2Value(c) - meta.TotalTroughValue(c);
    end

    meta.TotalAmplitudeMaxMinusMin(c) = max(thisTotal, [], "omitnan") - min(thisTotal, [], "omitnan");
    meta.NaNFractionTotal(c) = mean(isnan(thisTotal));
end

