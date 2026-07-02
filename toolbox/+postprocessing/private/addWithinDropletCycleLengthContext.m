% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function nucCycleInfo = addWithinDropletCycleLengthContext(nucCycleInfo, P)

    n = height(nucCycleInfo);
    CycleLengthFoldChangePrev = nan(n, 1);
    CycleLengthDeltaPrev = nan(n, 1);
    CycleLengthZWithinDroplet = nan(n, 1);
    IsWithinDropletLongCycle = false(n, 1);

    uniqueDropletIDs = unique(string(nucCycleInfo.DropletID));

    for i = 1:numel(uniqueDropletIDs)
        thisID = uniqueDropletIDs(i);
        idx = find(string(nucCycleInfo.DropletID) == thisID & nucCycleInfo.IsValidCycle);

        if numel(idx) < 2
            continue;
        end

        [~, order] = sort(nucCycleInfo.CycleInDroplet(idx));
        idx = idx(order);

        L = nucCycleInfo.CycleLength(idx);
        prevL = [NaN; L(1:end-1)];

        CycleLengthFoldChangePrev(idx) = L ./ prevL;
        CycleLengthDeltaPrev(idx) = L - prevL;

        medL = median(L, "omitnan");
        madL = mad(L, 1);

        if madL > 0
            CycleLengthZWithinDroplet(idx) = (L - medL) ./ (1.4826 * madL);
        end

        IsWithinDropletLongCycle(idx) = ...
            CycleLengthFoldChangePrev(idx) >= P.withinDropletFoldChangeThreshold | ...
            CycleLengthDeltaPrev(idx) >= P.withinDropletDeltaThreshold | ...
            CycleLengthZWithinDroplet(idx) >= P.withinDropletZThreshold;
    end

    nucCycleInfo.CycleLengthFoldChangePrev = CycleLengthFoldChangePrev;
    nucCycleInfo.CycleLengthDeltaPrev = CycleLengthDeltaPrev;
    nucCycleInfo.CycleLengthZWithinDroplet = CycleLengthZWithinDroplet;
    nucCycleInfo.IsWithinDropletLongCycle = IsWithinDropletLongCycle;
end

