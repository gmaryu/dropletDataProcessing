% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function [existingCycleTable, cycleBoundaryTable, traceInfo] = prepareCycleBoundaryTable(dataSet, traceInfo, nDroplets)

    existingCycleTable = dataSet.cycle;
    existingCycleTable.OriginalCycleRow = (1:height(existingCycleTable))';
    existingCycleTable.DropletID = makeDropletID(existingCycleTable.POS_ID, existingCycleTable.TRACK_ID);
    existingCycleTable.Diameter = NaN*ones(size(existingCycleTable,1),1);
    for i = 1:size(existingCycleTable,1)
        existingCycleTable.Diameter(i) = dataSet.info.ORIGINAL_MED_DIAMETER(dataSet.info.POS_ID == existingCycleTable.POS_ID(i) & dataSet.info.TRACK_ID == existingCycleTable.TRACK_ID(i));
    end
    %

    if ~ismember("DropletID", traceInfo.Properties.VariableNames)
        if all(ismember(["PositionNumber", "DropletNumber"], string(traceInfo.Properties.VariableNames)))
            traceInfo.DropletID = makeDropletID(traceInfo.PositionNumber, traceInfo.DropletNumber);
        else
            error("traceInfo must contain DropletID or PositionNumber/DropletNumber.");
        end
    end

    if ~ismember("TraceRow", traceInfo.Properties.VariableNames)
        error("traceInfo must contain TraceRow, the row index in the time-series matrices.");
    end

    checkTraceInfo(traceInfo, nDroplets);

    cycleBoundaryTable = table;
    cycleBoundaryTable.OriginalCycleRow = existingCycleTable.OriginalCycleRow;
    cycleBoundaryTable.Position = existingCycleTable.POS_ID;
    cycleBoundaryTable.Droplet = existingCycleTable.TRACK_ID;
    cycleBoundaryTable.CycleInDroplet = existingCycleTable.CYCLE_ID;
    cycleBoundaryTable.DropletID = existingCycleTable.DropletID;
    cycleBoundaryTable.Peak1Frame = existingCycleTable.START_INDEX;
    cycleBoundaryTable.Peak2Frame = existingCycleTable.END_INDEX;
    cycleBoundaryTable.TroughFrame = existingCycleTable.TROUGH_INDEX;

    optionalColumns = ["Peak1Value", "Peak2Value", "TroughValue"];
    for i = 1:numel(optionalColumns)
        col = optionalColumns(i);
        if ismember(col, string(existingCycleTable.Properties.VariableNames))
            cycleBoundaryTable.(col) = existingCycleTable.(col);
        end
    end

    cycleBoundaryTable = outerjoin( ...
        cycleBoundaryTable, ...
        traceInfo(:, ["DropletID", "TraceRow"]), ...
        "Keys", "DropletID", ...
        "MergeKeys", true);

    cycleBoundaryTable = sortrows(cycleBoundaryTable, "OriginalCycleRow");

    missingTrace = isnan(cycleBoundaryTable.TraceRow);
    fprintf("Cycles without matched TraceRow: %d / %d\n", ...
        sum(missingTrace), height(cycleBoundaryTable));

    if any(missingTrace)
        disp(cycleBoundaryTable(missingTrace, ...
            ["Position", "Droplet", "DropletID", "CycleInDroplet"]));
    end
end

