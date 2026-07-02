% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function checkTraceInfo(traceInfo, nDroplets)

    Gtrace = groupcounts(traceInfo, "DropletID");
    if any(Gtrace.GroupCount > 1)
        warning("Some DropletIDs map to multiple TraceRows. Inspect traceInfo before trusting classification.");
        disp(Gtrace(Gtrace.GroupCount > 1, :));
    end

    badRows = traceInfo.TraceRow < 1 | traceInfo.TraceRow > nDroplets | isnan(traceInfo.TraceRow);
    if any(badRows)
        warning("Some traceInfo.TraceRow values are outside the time-series matrix row range.");
        disp(traceInfo(badRows, :));
    end
end

%% =========================
%  Local functions: preprocessing
%  =========================

