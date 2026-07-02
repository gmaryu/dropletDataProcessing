% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function Xsmooth = smoothTraceMatrixNoBridge(Xraw, smoothWindow)

    Xsmooth = nan(size(Xraw));

    for i = 1:size(Xraw, 1)
        y = Xraw(i, :);

        if smoothWindow <= 1
            Xsmooth(i, :) = y;
            continue;
        end

        validIdx = ~isnan(y);
        d = diff([false, validIdx, false]);
        segStart = find(d == 1);
        segEnd = find(d == -1) - 1;

        for s = 1:numel(segStart)
            idx = segStart(s):segEnd(s);
            if numel(idx) < smoothWindow
                Xsmooth(i, idx) = y(idx);
            else
                Xsmooth(i, idx) = smoothdata(y(idx), "movmedian", smoothWindow);
            end
        end
    end
end

