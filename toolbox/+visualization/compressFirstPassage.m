function [compressed, complementary] = compressFirstPassage(osciPeakTable, thresFrame)
        pids = unique(osciPeakTable.POS_ID);
        dids = unique(osciPeakTable.TRACK_ID);

        compressed = table();
        complementary = table();

        for i = 1:length(pids)
            for j = 1:length(dids)
                dat = osciPeakTable(osciPeakTable.POS_ID == pids(i), :);
                dat = dat(dat.TRACK_ID == dids(j), :);

                minp = min(dat.END_FRAME - dat.START_FRAME, [], "omitnan");
                % unit in frame

                dat = sortrows(dat(dat.CYCLE_ID > 1, :), "CYCLE_ID");

                % initial period
                passflag = false;
                for k = 1:size(dat, 1)
                    if dat.END_FRAME(k) - dat.START_FRAME(k) > minp + thresFrame
                        passflag = true;
                        break;
                    end
                end
                if passflag == true
                    compressed(end + 1, :) = dat(k, :);
                    for m = 1:size(dat, 1)
                        if m ~= k
                            complementary(end + 1, :) = dat(m, :);
                        end
                    end
                end
            end
        end
end